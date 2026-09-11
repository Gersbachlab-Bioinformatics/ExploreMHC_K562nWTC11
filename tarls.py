#!/usr/bin/env python3
"""
List (or extract from) an uncompressed .tar served over HTTP, without
downloading the whole thing.

A tar file is a flat sequence of 512-byte headers, each followed by that
member's data padded to a 512-byte boundary. So if the transport can seek,
tarfile can read a header, consult its size field, and jump straight over
the payload to the next header. Over HTTP, "seek" means a Range request.

Result: listing a 50 GB archive reads a few hundred KB.

Requires the server to honor Range requests. Does NOT work on .tar.gz /
.tar.bz2 / .tar.xz -- those are single compressed streams with no seekable
member boundaries; stream those instead (see --help epilog).

Usage:
    python3 tarls.py URL                     # list members
    python3 tarls.py URL --long              # list with mode/owner/mtime
    python3 tarls.py URL --extract PATH      # write one member to stdout
"""

import argparse
import http.client
import io
import sys
import tarfile
import urllib.parse

CHUNK_CAP = 8 << 20  # never ask for more than 8 MB in a single Range request


class HTTPRangeFile(io.RawIOBase):
    """A seekable, read-only file object backed by HTTP Range requests."""

    def __init__(self, url, timeout=30):
        self.url = url
        self.timeout = timeout
        self._pos = 0
        self.bytes_read = 0
        self.requests = 0

        p = urllib.parse.urlsplit(url)
        if p.scheme not in ("http", "https"):
            raise ValueError(f"unsupported scheme: {p.scheme!r}")
        self._scheme, self._host = p.scheme, p.netloc
        self._path = p.path + (f"?{p.query}" if p.query else "")
        self._conn = None

        self.size, self.accepts_ranges = self._head()
        if self.size is None:
            raise RuntimeError("server did not report Content-Length; cannot seek")

    # -- connection handling -------------------------------------------------

    def _connect(self):
        cls = (http.client.HTTPSConnection if self._scheme == "https"
               else http.client.HTTPConnection)
        self._conn = cls(self._host, timeout=self.timeout)

    def _request(self, method, headers):
        """Issue a request, reconnecting once if the pooled socket went stale."""
        for attempt in (1, 2):
            try:
                if self._conn is None:
                    self._connect()
                self._conn.request(method, self._path, headers=headers)
                return self._conn.getresponse()
            except (http.client.HTTPException, OSError):
                try:
                    self._conn.close()
                except Exception:
                    pass
                self._conn = None
                if attempt == 2:
                    raise

    def _head(self):
        r = self._request("HEAD", {"Accept-Encoding": "identity"})
        r.read()
        if r.status >= 400:
            raise RuntimeError(f"HEAD {self.url} -> HTTP {r.status}")
        length = r.getheader("Content-Length")
        accepts = (r.getheader("Accept-Ranges") or "").lower() == "bytes"
        return (int(length) if length is not None else None), accepts

    # -- io.RawIOBase interface ---------------------------------------------

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
        n = min(n, self.size - self._pos, CHUNK_CAP)
        if n <= 0:
            return b""

        start = self._pos
        end = start + n - 1
        r = self._request("GET", {
            "Range": f"bytes={start}-{end}",
            "Accept-Encoding": "identity",
        })
        data = r.read()

        if r.status == 200:
            # Server ignored Range and sent the whole object. Salvage the slice
            # we wanted so the caller still gets correct bytes, but warn loudly:
            # every subsequent read will also transfer the entire file.
            print(f"warning: server ignored Range (HTTP 200, {len(data)} bytes); "
                  "listing will be slow", file=sys.stderr)
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


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def main():
    ap = argparse.ArgumentParser(
        description="List an uncompressed remote .tar via HTTP Range requests.",
        epilog=("For .tar.gz and friends there is no seeking; stream instead:\n"
                "  curl -s URL | tar -tzf -"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("url")
    ap.add_argument("-l", "--long", action="store_true",
                    help="show mode, owner, and mtime")
    ap.add_argument("-x", "--extract", metavar="MEMBER",
                    help="write one member's bytes to stdout")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="suppress the transfer summary")
    args = ap.parse_args()

    fobj = HTTPRangeFile(args.url)
    if not args.quiet:
        print(f"# object: {human(fobj.size)} ({fobj.size:,} bytes), "
              f"Accept-Ranges: {'yes' if fobj.accepts_ranges else 'not advertised'}",
              file=sys.stderr)

    count = total = 0
    try:
        # mode="r:" forces uncompressed; it stops tarfile from trying gzip/bz2
        # decoders, which would read the whole stream.
        with tarfile.open(fileobj=fobj, mode="r:") as tf:
            if args.extract:
                member = tf.getmember(args.extract)
                src = tf.extractfile(member)
                if src is None:
                    sys.exit(f"{args.extract!r} is not a regular file")
                out = sys.stdout.buffer
                while True:
                    buf = src.read(1 << 20)
                    if not buf:
                        break
                    out.write(buf)
                out.flush()
            else:
                for m in tf:
                    count += 1
                    total += m.size
                    if args.long:
                        kind = "d" if m.isdir() else ("l" if m.issym() else "-")
                        print(f"{kind}{oct(m.mode)[2:]:>4}  "
                              f"{m.uname or m.uid}/{m.gname or m.gid:<8}  "
                              f"{m.size:>12,}  "
                              f"{m.mtime and __import__('datetime').datetime.fromtimestamp(m.mtime):%Y-%m-%d %H:%M}  "
                              f"{m.name}")
                    else:
                        print(f"{m.size:>12,}  {m.name}")
    except tarfile.ReadError as e:
        sys.exit(f"not a readable uncompressed tar: {e}\n"
                 "(if the URL ends in .tar.gz/.tgz, stream it instead: "
                 "curl -s URL | tar -tzf -)")
    finally:
        pct = 100 * fobj.bytes_read / fobj.size if fobj.size else 0
        if not args.quiet:
            if not args.extract:
                print(f"# {count} members, {human(total)} uncompressed",
                      file=sys.stderr)
            print(f"# transferred {human(fobj.bytes_read)} in {fobj.requests} "
                  f"requests ({pct:.2f}% of object)", file=sys.stderr)
        fobj.close()


if __name__ == "__main__":
    main()
