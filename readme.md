# Gene Ontology Enrichment Analysis Pipeline: From Annotations to Dotplot

## Overview

This pipeline extracts transcription start sites (TSS) from human gene annotation data, identifies promoter regions by extending these coordinates upstream, extracts the DNA sequences for these promoters, and searches for a specific GC-rich transcription factor binding motif (`GCGC..GCGC`). Finally, it extracts the genes containing this motif and performs a Gene Ontology (GO) enrichment analysis to visualize their biological functions.

---

## Prerequisites

- **Python 3.x**
- **bedtools**
- **EMBOSS suite (`dreg`)**
- **Conda** with Bioconductor packages:

```
conda install -c bioconda bioconductor-clusterprofiler bioconductor-org.hs.eg.db
```

- **Required Reference Files:**
  - `human_gene_annotation.tsv` (Ensembl human gene annotations)
  - `hg38.chrom.sizes` (UCSC chromosome sizes for hg38)
  - `hg38.fa` (Reference human genome sequence)

---

## Step 1: Extracting Transcription Start Sites (TSS)

### 1A: Filtering for Valid Chromosomes

```
cat human_gene_annotation.tsv | python filter.py > filtered.tsv
```

### 1B: Transforming to BED Format

```
cat filtered.tsv | python transform.py > hg38_tss.bed
```

---

## Step 2: Defining Promoter Regions

### 2A: Harmonizing Chromosome Names

```
awk 'BEGIN{FS="\t"; OFS="\t"} {$1 = "chr"$1; print $0}' hg38_tss.bed > hg38_tss_chr.bed
```

### 2B: Extending Coordinates with `bedtools slop`

```
bedtools slop -i hg38_tss_chr.bed -g hg38.chrom.sizes -l 500 -r 0 -s > hg38_promoters_500bp.bed
```

---

## Step 3: Extracting Promoter Sequences

```
bedtools getfasta -fi hg38.fa -bed hg38_promoters_500bp.bed -fo hg38_promoters_500bp.fa -name
```

---

## Step 4: Motif Searching

```
dreg -sequence hg38_promoters_500bp.fa -pattern "GCGC..GCGC" -outfile hg38_promoters_500bp_GCGC_hits.txt
```

---

## Step 5: Extracting the Target Gene List

### 5A: Extract Positive Hit Coordinates

```
awk '/Sequence:/ {seq=$3} /HitCount:/ {if ($3 > 0) print seq}' hg38_promoters_500bp_GCGC_hits.txt > hit_coords.txt
```

### 5B: Map Back to BED and Extract Gene Symbols

```
awk 'NR==FNR {hits[$1]=1; next} {coord=$2"-"$3; if (coord in hits) print $4}' hit_coords.txt hg38_promoters_500bp.bed | cut -d '|' -f 2 | sort -u > hg38_GCGC_target_genes.tsv
```

---

## Step 6: Gene Ontology (GO) Enrichment Analysis

```
Rscript go_analysis.R
```

**Final Outputs:**

| File | Description |
|---|---|
| `hg38_GCGC_GO_results.csv` | Raw tabular data of enriched biological processes |
| `hg38_GCGC_GO_dotplot.pdf` | Top 20 enriched biological processes dotplot |
