# Capsella bursa‑pastoris — Self‑synteny & homeolog analysis
This study shows the synteny analysis performed on the Capsella bursa‑pastoris genome (an allopolyploid). The analysis focuses on detecting orthologs, paralogs, and homeologs and on testing whether chromosomes in the two subgenomes retain expected homeologous relationships (1 ↔ 9, 2 ↔ 10, …, 8 ↔ 16).

Capsella bursa‑pastoris is an allopolyploid formed from two diploid parental lineages, producing two subgenomes commonly referred to as A and B. Here we treat:

- Subgenome A: chromosomes jlSCF_1 … jlSCF_8  
- Subgenome B: chromosomes jlSCF_9 … jlSCF_16

The main goals:
- Identify gene matches and syntenic blocks within the genome
- Test chromosome pairings between A and B
- Visualize cross‑group homology with a circos plot and a chromosome–chromosome heatmap

---

## Data & tools

- Input: Capsella bursa‑pastoris genome FASTA and GFF3 annotation
- GENESPACE (R package) — used to run OrthoFinder, BLAST, and MCScanX and to produce downstream R objects and tables
- R and the following packages:
  - dplyr, circlize, scales, RColorBrewer, ggplot2
- GENESPACE output files used in this analysis:
  - `Capsella_vs_Capsella.synHits.txt.gz`
  - `Capsella_vs_Capsella.allBlast.txt.gz`
  - `syntenicBlock_coordinates.csv`
  - `Capsella_geneOrder_rSourceData.rda`
  - `Capsella_bp_rSourceData.rda`

---

## Workflow (high level)

1. Run GENESPACE on the HPC cluster with a config pointing to the Capsella genome and annotation. This runs OrthoFinder, BLAST, and MCScanX internally and writes synteny and BLAST outputs.
2. Load GENESPACE outputs into R for downstream filtering, counting, and visualization.
3. Remove contigs/unplaced scaffolds and focus only on canonical chromosomes.
4. Define subgenome groups (A = 1–8, B = 9–16).
5. Subset self‑synteny hits to cross‑group matches (A vs B) and count matches per chromosome pair.
6. Summarize block coordinates to look for large syntenic blocks connecting expected homeologous chromosomes.
7. Visualize cross‑group homology using a circos plot and a chromosome–chromosome heatmap.

---

## Key R snippets used in the analysis

Below is the R code used to create the circos plot and the chromosome–chromosome heatmap (heatmap uses counts of all BLAST/synteny hits between chromosomes). Update `path` to where your GENESPACE outputs are stored (e.g., Downloads or cluster directory). The code includes a simple filter to remove contigs and scaffolds and organizes chromosomes numerically.

```r
library(dplyr)
library(circlize)
library(scales)
library(RColorBrewer)
library(ggplot2)

# ============================================================
# Load BLAST data
# ============================================================
path <- "~/Downloads"

blast <- read.table(
  file.path(path, "Capsella_vs_Capsella.allBlast.txt.gz"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ============================================================
# Clean data (remove scaffolds / contigs)
# ============================================================
blast <- blast %>%
  filter(
    !grepl("contig|scaffold", chr1, ignore.case = TRUE),
    !grepl("contig|scaffold", chr2, ignore.case = TRUE)
  )

# Assign chromosome groups (1-8 = A, 9-16 = B)

get_group <- function(chr) {
  num <- as.numeric(gsub("[^0-9]", "", chr))
  if (num >= 1 & num <= 8) return("A")
  if (num >= 9 & num <= 16) return("B")
  return(NA)
}

blast <- blast %>%
  mutate(
    group1 = sapply(chr1, get_group),
    group2 = sapply(chr2, get_group)
  ) %>%
  # Keep only cross-group matches
  filter(!is.na(group1) & !is.na(group2) & group1 != group2)


# Robust chromosome ordering (1 → 16)

chroms <- unique(c(blast$chr1, blast$chr2))
chrom_nums <- as.numeric(gsub("[^0-9]", "", chroms))
chroms <- chroms[order(chrom_nums)]

blast$chr1 <- factor(blast$chr1, levels = chroms)
blast$chr2 <- factor(blast$chr2, levels = chroms)

# ============================================================
# ---------------- CIRCOS PLOT -----------------------------
# ============================================================
# Subsample for visibility
set.seed(123)
blast_small <- blast %>%
  sample_n(min(15000, nrow(.)))

n <- length(chroms)
chrom_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n),
  chroms
)

circos.clear()
circos.par(start.degree = 90, gap.degree = 3, track.margin = c(0.01, 0.01))
circos.initialize(factors = chroms, xlim = cbind(rep(0, n), rep(1, n)))

# Chromosome ideograms
circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                col = chrom_colors[chr], border = NA)
    circos.text(CELL_META$xcenter, 0.5, chr, facing = "inside", cex = 0.6)
  }
)

# Links
for (i in seq_len(nrow(blast_small))) {
  circos.link(
    blast_small$chr1[i], runif(1),
    blast_small$chr2[i], runif(1),
    col = adjustcolor(chrom_colors[blast_small$chr1[i]], alpha.f = 0.1),
    border = NA
  )
}

# ============================================================
# -------------- CHROMOSOME-CHROMOSOME HEATMAP -------------
# ============================================================
chr_matrix <- blast %>%
  count(chr1, chr2)

ggplot(chr_matrix, aes(x = chr1, y = chr2, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "#fde0dd", high = "#a50f15", trans = "log10", guide = "none") +
  coord_fixed() +
  labs(
    title = "Cross-group homology (1–8 vs 9–16)",
    x = "Chromosome",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )
```

---

## Notes on data cleaning and subgenome definition

- Remove contigs and unplaced scaffolds before counting chromosome pairs to avoid spurious matches.
- Define subgenomes by chromosome number ranges (1–8 A, 9–16 B) and restrict tests to cross‑group matches. This makes the homeologous signal clearer.
- When summarizing syntenic blocks, filter for large block sizes and/or high gene density to prioritize biologically meaningful matches.

---

## Block‑level analysis

Load `syntenicBlock_coordinates.csv` (output from GENESPACE/MCScanX) and summarize block lengths and their connecting chromosomes. Example (in R):

```r
blocks <- read.csv("~/Downloads/syntenicBlock_coordinates.csv")
# Filter to cross-group blocks and summarize length / gene counts
blocks_AB <- blocks %>%
  filter((chr1 %in% paste0("jlSCF_", 1:8) & chr2 %in% paste0("jlSCF_", 9:16)) |
         (chr1 %in% paste0("jlSCF_", 9:16) & chr2 %in% paste0("jlSCF_", 1:8)))

blocks_AB %>%
  group_by(chr1, chr2) %>%
  summarize(
    n_blocks = n(),
    total_len = sum(block_length, na.rm = TRUE),
    mean_len = mean(block_length, na.rm = TRUE)
  ) %>%
  arrange(desc(total_len))
```

Note:(Adjust column names as needed to match your `syntenicBlock_coordinates.csv`.)
<img width="1259" height="778" alt="image" src="https://github.com/user-attachments/assets/ceeaae7d-4950-4c0b-be43-4242cac48a8e" /> 
<img width="1259" height="778" alt="image" src="https://github.com/user-attachments/assets/b6af650f-4640-4db2-a6a5-1e686b09d8c2" />


---

## Interpretation & expectations

Given the allopolyploid origin:
- Expect strong syntenic connections between homeologous chromosome pairs (1 ↔ 9, 2 ↔ 10, …).
- The heatmap should show strong off‑diagonal tiles corresponding to these pairings.
- The circos plot should show abundant links between corresponding A and B chromosomes.
- Block‑level summaries should show larger/denser blocks linking homeologous pairs compared to other cross‑group combinations.

---

## Recommended next steps

1. Commit this Markdown file to your repo (suggested filename `synteny-analysis.md`).
2. Run the R script on the full BLAST/synHits files on the cluster to generate high‑resolution figures (save as PNG/PDF).
3. Add the generated figures to the repo (e.g., `figures/circos.png`, `figures/chr_heatmap.png`) and reference them in this Markdown.
4. Summarize any deviations from expectation (e.g., translocations, fusions, weak homeology) and add a discussion section with figures.
5. Optional: add a small table of chr1 ↔ chr9 counts using the counted results (from `syn_AB %>% count(chr1, chr2)`).

