# Capsella RNA-seq Expression Pipeline — Status & Notes

Status (as of 2025-12-19)
- Data generation & preprocessing: ✅ Completed
- Reference genome prep: ✅ Completed
- Alignment to reference: ✅ Completed
- Gene-level counting: ✅ Completed
- Data transfer for local analysis: ✅ Completed
- Exploratory QC (PCA): ✅ Completed
- Next major analysis: Introgression → Differential expression / allele-specific analyses

Project goal
-----------
Generate high-quality, tissue-resolved RNA-seq expression datasets in Capsella and use them to study how introgression may affect gene expression. This repo documents completed preprocessing pipelines, QC, and outlines the next analysis steps.

High-level summary of completed pipelines
----------------------------------------
- 1) SRR download → FASTQ (paired-end)  
  - Tool: SRA Toolkit (`fasterq-dump`) run as a SLURM array  
  - Output: gzipped FASTQ files organized by tissue:  
    - `flower/raw_fastq/*.fastq.gz`  
    - `leaf/raw_fastq/*.fastq.gz`  
    - `root/raw_fastq/*.fastq.gz`
  - Status: ✔ All SRRs downloaded and compressed

- 2) Reference genome preparation  
  - Files used: `CBP_JLv4.fa` (genome), `CBP_JLv4_augustus.gff3` (annotation)  
  - Indexes created for alignment (HISAT2) and FASTA indexed (`.fai`)  
  - Status: ✔ Reference indexed and ready

- 3) RNA-seq alignment → coordinate-sorted BAMs  
  - Tool(s): HISAT2 (index), `samtools` for sorting/indexing/validation  
  - Output: `alignments/flower/*.bam`, `alignments/leaf/*.bam`, `alignments/root/*.bam`  
  - Status: ✔ All samples aligned; BAMs validated and readable; consistent depth

- 4) Gene-level read counting (featureCounts)  
  - Tool: `featureCounts` (Subread v2.0.6)  
  - Key params: paired-end mode (`-p`), exon features (`-t exon`), gene grouping by Parent attribute in GFF3  
  - Outputs: `flower_counts.txt`, `leaf_counts.txt`, `root_counts.txt` + `.summary` files  
  - Typical assignment rates: ~50–65% (within expected range for plant RNA-seq)  
  - Status: ✔ Counts generated

- 5) Data transfer to local machine for R analysis  
  - Files copied: counts and annotation (`CBP_JLv4_augustus.gff3`)  
  - Purpose: downstream statistics and visualization in R / RStudio  
  - Status: ✔ Files available locally

Exploratory QC (completed)
--------------------------
- PCA performed on normalized expression (e.g., logCPM / vst / rlog)  
- Observations:
  - Strong tissue separation
  - PC1: root vs aerial tissues
  - PC2: leaf vs flower
  - Replicates cluster tightly
- Conclusion: Expression data quality is high and biologically meaningful

Exact script used for SRR download (SLURM array)
------------------------------------------------
The SLURM array script used to download SRRs and compress FASTQs exactly as run:

```bash
# Load SRA Toolkit
module load SRA-Toolkit/3.0.10-gompi-2023a
BASE_DIR=/mnt/scratch/f0111256/capsella_project
cd "$BASE_DIR"

# Build a combined list for indexing this array task
combined_list=$(mktemp)
cat flower_srr.txt leaf_srr.txt root_srr.txt > "$combined_list"
total_lines=$(wc -l < "$combined_list")

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > total_lines )); then
  echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} out of range (1..${total_lines})"
  rm -f "$combined_list"
  exit 1
fi

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$combined_list")
echo "Task $SLURM_ARRAY_TASK_ID processing SRR: $SRR"

# Determine tissue
tissue=""
for t in flower leaf root; do
  if grep -qxF "$SRR" "${t}_srr.txt"; then
    tissue="$t"
    break
  fi
done

[[ -z "$tissue" ]] && tissue="unknown"
outdir="$BASE_DIR/$tissue/raw_fastq"
mkdir -p "$outdir"

# Skip if already downloaded
if [[ -f "$outdir/${SRR}_1.fastq.gz" && -f "$outdir/${SRR}_2.fastq.gz" ]]; then
  echo "Files for $SRR already exist, skipping."
  rm -f "$combined_list"
  exit 0
fi

# Download
fasterq-dump "$SRR" --split-files --threads "$SLURM_CPUS_PER_TASK" --outdir "$outdir"

# Compress
if command -v pigz >/dev/null 2>&1; then
  pigz -p "$SLURM_CPUS_PER_TASK" "$outdir/${SRR}"_*.fastq
else
  gzip "$outdir/${SRR}"_*.fastq
fi

rm -f "$combined_list"
echo "Finished $SRR"
```

Reproducible notes & locations
------------------------------
- Raw FASTQs: `<BASE_DIR>/{flower,leaf,root}/raw_fastq/*.fastq.gz`  
- Reference genome & annotation: `<BASE_DIR>/references/CBP_JLv4.fa`, `CBP_JLv4_augustus.gff3`  
- Alignment outputs: `<BASE_DIR>/alignments/{flower,leaf,root}/*.bam`  
- Count matrices: `<BASE_DIR>/counts/{flower_counts.txt, leaf_counts.txt, root_counts.txt}`  
- Local copies: these count files + `CBP_JLv4_augustus.gff3` were transferred to the analyst's local machine for R

Tools / environment (used)
--------------------------
- SRA Toolkit (fasterq-dump) — for SRR → FASTQ  
- pigz / gzip — multithreaded / single-threaded compression  
- HISAT2 — alignment (indexes prepared from CBP_JLv4.fa)  
- samtools — BAM sorting, indexing, validation  
- Subread / featureCounts v2.0.6 — gene-level read counting  
- SLURM — job array for parallel SRR downloads  
- R / Bioconductor — downstream normalization, PCA, visualization

Key QC numbers & interpretation
-------------------------------
- Read assignment rate to genes: ~50–65%  
  - Interpretation: typical for plant RNA-seq when counting to exon features from a GFF3 annotation; not unusually low
- PCA: strong tissue separation and tight replicates indicate robust biological signal and low technical noise.<img width="1260" height="778" alt="image" src="https://github.com/user-attachments/assets/a1f58da7-19f5-468e-9ecc-aef08c5f9938" />
-Heatmaps. <img width="1260" height="778" alt="image" src="https://github.com/user-attachments/assets/dbdecc94-f2bd-428f-b4b0-ef7cc803530e" />


Next steps — introgression & expression analyses (recommended)
------------------------------------------------------------
Planned high-level analyses to study introgression effects on expression:
1. Sample metadata & covariates
   - Ensure complete sample metadata (population/line, tissue, replicate, batch, sequencing lane, SRR).
2. Normalize & filter
   - Use e.g. DESeq2 (rlog / vst) or edgeR (TMM + voom) for normalization and variance-stabilizing transforms.
   - Filter lowly expressed genes (e.g., keep genes with CPM > 1 in at least X samples).
3. Differential expression (DE)
   - Model: expression ~ introgression_status + tissue + batch + (other covariates)
   - Perform tissue-specific DE and cross-tissue comparisons.
4. Allele-specific expression (ASE)
   - If phased variants / read-backed phasing available, quantify ASE to detect expression bias from introgressed haplotypes.
   - Tools: WASP, GATK ASEReadCounter, alleleCounts pipelines.
5. Local ancestry / introgression mapping
   - Use available genotype / ancestry calls to map local ancestry and correlate local ancestry with gene expression.
   - Methods: eQTL-like tests using ancestry as predictor; linear mixed models to account for relatedness.
6. Integrative visualizations
   - Heatmaps, PCA per tissue, volcano plots, gene set enrichment for introgression-associated DE genes.
7. Validation / follow-up
   - qPCR validation for top candidate genes, and cross-checks across independent populations if available.

How to reproduce the preprocessing quickly
-----------------------------------------
- Requirements: SLURM cluster, SRA Toolkit, HISAT2, samtools, Subread/featureCounts, pigz (recommended)
- Typical commands (high-level):
  - fasterq-dump SRR --split-files --threads N --outdir PATH
  - hisat2 -x ref_index -1 sample_1.fastq.gz -2 sample_2.fastq.gz | samtools sort -o sample.bam
  - samtools index sample.bam
  - featureCounts -T N -p -t exon -g Parent -a CBP_JLv4_augustus.gff3 -o counts.txt *.bam. 

