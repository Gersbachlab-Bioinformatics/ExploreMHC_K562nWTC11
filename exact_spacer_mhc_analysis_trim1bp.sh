#!/usr/bin/env bash
set -euo pipefail
# Usage: bash exact_spacer_mhc_analysis_trim1bp.sh IGVFFI2404DYFG.tsv.gz GRCh38.fa output_dir
INPUT_TSV_GZ="${1:?Input IGVFFI TSV.GZ required}"
GRCH38_FA="${2:?GRCh38 FASTA required}"
OUTDIR="${3:-spacer_mhc_exact_trim1bp}"
K562_URL='https://raw.githubusercontent.com/Haozhe-Yuan/MHC_Haplotype_Assemble/refs/heads/main/K562.hap2.MHC.fa'
THREADS="${THREADS:-8}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUTDIR"/{reference,index,sam,results,plots,tmp}
for cmd in bowtie2 bowtie2-build python3 wget; do command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing $cmd" >&2; exit 1; }; done
K562_FA="$OUTDIR/reference/K562.hap2.MHC.fa"
[[ -s "$K562_FA" ]] || wget -O "$K562_FA" "$K562_URL"
python3 "$SCRIPT_DIR/extract_spacers_trim1bp.py" "$INPUT_TSV_GZ" "$OUTDIR/tmp/unique_trimmed_19nt_spacers.fa" "$OUTDIR/tmp/trimmed_spacer_row_map.tsv.gz"
[[ -e "$OUTDIR/index/grch38.1.bt2" || -e "$OUTDIR/index/grch38.1.bt2l" ]] || bowtie2-build --threads "$THREADS" "$GRCH38_FA" "$OUTDIR/index/grch38"
[[ -e "$OUTDIR/index/k562_mhc.1.bt2" || -e "$OUTDIR/index/k562_mhc.1.bt2l" ]] || bowtie2-build --threads "$THREADS" "$K562_FA" "$OUTDIR/index/k562_mhc"
bowtie2 --end-to-end -a --no-unal -f -p "$THREADS" -x "$OUTDIR/index/grch38" -U "$OUTDIR/tmp/unique_trimmed_19nt_spacers.fa" -S "$OUTDIR/sam/grch38_trim1bp.sam"
bowtie2 --end-to-end -a --no-unal -f -p "$THREADS" -x "$OUTDIR/index/k562_mhc" -U "$OUTDIR/tmp/unique_trimmed_19nt_spacers.fa" -S "$OUTDIR/sam/k562_mhc_trim1bp.sam"
python3 "$SCRIPT_DIR/parse_exact_alignments_trim1bp.py" "$OUTDIR/sam/grch38_trim1bp.sam" GRCh38 "$OUTDIR/results/grch38_exact_alignments_trim1bp.tsv.gz"
python3 "$SCRIPT_DIR/parse_exact_alignments_trim1bp.py" "$OUTDIR/sam/k562_mhc_trim1bp.sam" K562_hap2_MHC "$OUTDIR/results/k562_mhc_exact_alignments_trim1bp.tsv.gz"
python3 "$SCRIPT_DIR/summarize_and_cluster_exact_trim1bp.py" "$INPUT_TSV_GZ" "$OUTDIR/tmp/trimmed_spacer_row_map.tsv.gz" "$OUTDIR/results/grch38_exact_alignments_trim1bp.tsv.gz" "$OUTDIR/results/k562_mhc_exact_alignments_trim1bp.tsv.gz" "$OUTDIR/results" "$OUTDIR/plots"
if command -v mafft >/dev/null 2>&1; then mafft --auto "$OUTDIR/results/k562_matching_unique_trimmed_spacers.fa" > "$OUTDIR/results/k562_matching_unique_trimmed_spacers.mafft.fa"; else echo 'NOTE: MAFFT not installed; MSA skipped.' >&2; fi
echo 'Completed.'
echo "$OUTDIR/results/guide_exact_match_summary_trim1bp.tsv.gz"
echo "$OUTDIR/results/exact_match_category_counts_trim1bp.tsv"
echo "$OUTDIR/plots/k562_trim1bp_spacer_hamming_dendrogram.png"
echo "$OUTDIR/plots/k562_mhc_trim1bp_spacer_positions.png"
