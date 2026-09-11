#!/usr/bin/env python3
"""
List or extract members of an uncompressed .tar served over HTTP, without
downloading the whole archive.

Two phases, deliberately different:

  Listing   -- header-hopping. A tar is 512-byte headers each followed by
               padded data, so tarfile can read a header, consult its size
               field, and Range-seek past the payload. Reads a few hundred KB
               regardless of archive size.

  Extract   -- ONE streaming request. Once listing has located the member's
               offset_data, we issue a single open-ended Range GET and read
               the body continuously. This is the important difference from
               v1, which issued a Range request per 1 MB read and produced
               thousands of connections per member -- enough for a server to
               start refusing them.

If the connection drops mid-stream, extraction resumes from the last byte
received rather than starting over.

Does NOT work on .tar.gz / .tar.bz2 / .tar.xz -- single compressed streams
have no seekable member boundaries. Stream those: curl -s URL | tar -tzf -

Usage:
    python3 tarls2.py URL                    # list members
    python3 tarls2.py URL --long             # list with mode/owner/mtime
    python3 tarls2.py URL -x PATH            # stream one member to stdout
    python3 tarls2.py URL -x PATH -o FILE    # ...to a file, resumable
"""

import argparse
import datetime
import difflib
import http.client
import io
import os
import sys
import tarfile
import time
import urllib.parse

HEADER_CHUNK = 1 << 20   # cap on a single listing-phase Range request
STREAM_CHUNK = 4 << 20   # read granularity while streaming a member body


def _split(url):
    p = urllib.parse.urlsplit(url)
    if p.scheme not in ("http", "https"):
        raise ValueError(f"unsupported scheme: {p.scheme!r}")
    return p.scheme, p.netloc, p.path + (f"?{p.query}" if p.query else "")


def _connect(scheme, host, timeout):
    cls = (http.client.HTTPSConnection if scheme == "https"
           else http.client.HTTPConnection)
    return cls(host, timeout=timeout)


class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only file object over HTTP Range. Used for LISTING only."""

    def __init__(self, url, timeout=60, retries=5):
        self.url, self.timeout, self.retries = url, timeout, retries
        self._pos = 0
        self.bytes_read = 0
        self.requests = 0
        self._scheme, self._host, self._path = _split(url)
        self._conn = None
        self.size, self.accepts_ranges = self._head()
        if self.size is None:
            raise RuntimeError("server did not report Content-Length; cannot seek")

    def _request(self, method, headers):
        """Issue a request, retrying with exponential backoff on transport errors."""
        last = None
        for attempt in range(self.retries):
            try:
                if self._conn is None:
                    self._conn = _connect(self._scheme, self._host, self.timeout)
                self._conn.request(method, self._path, headers=headers)
                return self._conn.getresponse()
            except (http.client.HTTPException, OSError) as e:
                last = e
                if self._conn is not None:
                    try:
                        self._conn.close()
                    except Exception:
                        pass
                    self._conn = None
                if attempt < self.retries - 1:
                    delay = min(2 ** attempt, 30)
                    print(f"warning: {type(e).__name__}: {e}; retrying in {delay}s "
                          f"({attempt + 1}/{self.retries})", file=sys.stderr)
                    time.sleep(delay)
        raise RuntimeError(f"request failed after {self.retries} attempts: {last}")

    def _head(self):
        r = self._request("HEAD", {"Accept-Encoding": "identity"})
        r.read()
        if r.status >= 400:
            raise RuntimeError(f"HEAD {self.url} -> HTTP {r.status}")
        length = r.getheader("Content-Length")
        accepts = (r.getheader("Accept-Ranges") or "").lower() == "bytes"
        return (int(length) if length is not None else None), accepts

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self._pos

    def seek(self, offset, whence=io.SEEK_SET):
        if whence == io.SEEK_SET:
            new = offset
        elif whence == io.SEEK_CUR:
            new = self._pos + offset
        elif whence == io.SEEK_END:
            new = self.size + offset
        else:
            raise ValueError(f"bad whence: {whence}")
        self._pos = max(0, new)
        return self._pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self._pos
        n = min(n, self.size - self._pos, HEADER_CHUNK)
        if n <= 0:
            return b""
        start, end = self._pos, self._pos + n - 1
        r = self._request("GET", {"Range": f"bytes={start}-{end}",
                                  "Accept-Encoding": "identity"})
        data = r.read()
        if r.status == 200:
            print("warning: server ignored Range (HTTP 200); listing will be slow",
                  file=sys.stderr)
            data = data[start:end + 1]
        elif r.status != 206:
            raise RuntimeError(f"GET range {start}-{end} -> HTTP {r.status}")
        self.requests += 1
        self.bytes_read += len(data)
        self._pos += len(data)
        return data

    def close(self):
        if self._conn is not None:
            try:
                self._conn.close()
            except Exception:
                pass
            self._conn = None
        super().close()


def stream_range(url, start, length, out, timeout=60, retries=8, progress=True):
    """
    Stream bytes [start, start+length) to `out` using ONE request, resuming
    from the last received byte if the connection drops.
    """
    scheme, host, path = _split(url)
    done = 0
    attempt = 0
    t0 = time.time()
    last_report = 0.0

    while done < length:
        conn = None
        try:
            conn = _connect(scheme, host, timeout)
            lo = start + done
            hi = start + length - 1
            conn.request("GET", path, headers={"Range": f"bytes={lo}-{hi}",
                                               "Accept-Encoding": "identity"})
            r = conn.getresponse()
            if r.status not in (200, 206):
                raise RuntimeError(f"GET range {lo}-{hi} -> HTTP {r.status}")
            if r.status == 200 and done > 0:
                raise RuntimeError("server ignored Range; cannot resume")

            while done < length:
                buf = r.read(min(STREAM_CHUNK, length - done))
                if not buf:
                    break
                out.write(buf)
                done += len(buf)
                attempt = 0  # progress made; reset the failure counter
                now = time.time()
                if progress and now - last_report > 5:
                    rate = done / max(now - t0, 1e-9) / (1 << 20)
                    print(f"\r  {human(done)} / {human(length)} "
                          f"({100 * done / length:.1f}%)  {rate:.1f} MB/s",
                          end="", file=sys.stderr, flush=True)
                    last_report = now

            if done < length:
                raise RuntimeError(f"stream ended early at {done}/{length} bytes")

        except (http.client.HTTPException, OSError, RuntimeError) as e:
            attempt += 1
            if attempt >= retries:
                raise RuntimeError(
                    f"failed at byte {start + done} after {retries} attempts: {e}")
            delay = min(2 ** attempt, 60)
            print(f"\nwarning: {type(e).__name__}: {e}\n"
                  f"  resuming from byte {done:,} in {delay}s "
                  f"({attempt}/{retries})", file=sys.stderr)
            time.sleep(delay)
        finally:
            if conn is not None:
                try:
                    conn.close()
                except Exception:
                    pass

    if progress:
        print(f"\r  {human(done)} complete{' ' * 30}", file=sys.stderr)
    return done


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def main():
    ap = argparse.ArgumentParser(
        description="List or extract from a remote uncompressed .tar via HTTP Range.",
        epilog="For .tar.gz: curl -s URL | tar -tzf -",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("url")
    ap.add_argument("-l", "--long", action="store_true",
                    help="show mode, owner, and mtime")
    ap.add_argument("-x", "--extract", metavar="MEMBER",
                    help="stream one member (to stdout unless -o given)")
    ap.add_argument("-o", "--output", metavar="FILE",
                    help="write extracted member here; enables resume if interrupted")
    ap.add_argument("--timeout", type=float, default=60)
    ap.add_argument("--retries", type=int, default=8)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()

    fobj = HTTPRangeFile(args.url, timeout=args.timeout, retries=args.retries)
    if not args.quiet:
        print(f"# object: {human(fobj.size)} ({fobj.size:,} bytes), "
              f"Accept-Ranges: {'yes' if fobj.accepts_ranges else 'not advertised'}",
              file=sys.stderr)

    try:
        with tarfile.open(fileobj=fobj, mode="r:") as tf:
            if args.extract:
                try:
                    m = tf.getmember(args.extract)
                except KeyError:
                    names = [n for n in tf.getnames()]
                    near = difflib.get_close_matches(args.extract, names, n=8, cutoff=0.5)
                    if not near:
                        tail = args.extract.rsplit("/", 1)[-1]
                        near = difflib.get_close_matches(
                            tail, [n.rsplit("/", 1)[-1] for n in names], n=8, cutoff=0.4)
                    msg = f"member not found: {args.extract!r}"
                    if near:
                        msg += "\ndid you mean one of:\n  " + "\n  ".join(near)
                    sys.exit(msg)
                if not m.isfile():
                    sys.exit(f"{args.extract!r} is not a regular file")
                offset, size = m.offset_data, m.size
            else:
                count = total = 0
                for m in tf:
                    count += 1
                    total += m.size
                    if args.long:
                        kind = "d" if m.isdir() else ("l" if m.issym() else "-")
                        ts = datetime.datetime.fromtimestamp(m.mtime)
                        print(f"{kind}{oct(m.mode)[2:]:>4}  "
                              f"{m.uname or m.uid}/{m.gname or m.gid:<8}  "
                              f"{m.size:>14,}  {ts:%Y-%m-%d %H:%M}  {m.name}")
                    else:
                        print(f"{m.size:>14,}  {m.name}")
                if not args.quiet:
                    print(f"# {count} members, {human(total)} uncompressed",
                          file=sys.stderr)
                    print(f"# transferred {human(fobj.bytes_read)} in "
                          f"{fobj.requests} requests "
                          f"({100 * fobj.bytes_read / fobj.size:.2f}% of object)",
                          file=sys.stderr)
                return
    except tarfile.ReadError as e:
        sys.exit(f"not a readable uncompressed tar: {e}")
    finally:
        fobj.close()

    # --- extraction: single streaming request, outside the tarfile context ---
    if not args.quiet:
        print(f"# streaming {args.extract} ({human(size)}) "
              f"from offset {offset:,}", file=sys.stderr)

    if args.output:
        # Resume support: if a partial file exists, continue from its length.
        existing = os.path.getsize(args.output) if os.path.exists(args.output) else 0
        if existing > size:
            sys.exit(f"{args.output} is larger than the member ({existing} > {size}); "
                     "remove it and retry")
        if existing == size:
            if not args.quiet:
                print("# already complete", file=sys.stderr)
            return
        if existing and not args.quiet:
            print(f"# resuming from {human(existing)}", file=sys.stderr)
        with open(args.output, "ab") as fh:
            stream_range(args.url, offset + existing, size - existing, fh,
                         timeout=args.timeout, retries=args.retries,
                         progress=not args.quiet)
    else:
        stream_range(args.url, offset, size, sys.stdout.buffer,
                     timeout=args.timeout, retries=args.retries,
                     progress=not args.quiet)
        sys.stdout.buffer.flush()


if __name__ == "__main__":
    main()
