# Promoter Region Extraction Pipeline
### hg38 | 500 bp upstream | Strand-aware | BEDTools

---

## What This Does

Takes raw human gene annotation data and produces **strand-correct promoter windows** (500 bp upstream of each TSS) ready for motif scanning or GO enrichment.

```
human_gene_annotation.tsv.gz
        │
        ▼
  [Parse TSS positions]
        │
        ▼
  [Filter to valid chromosomes]
        │
        ▼
  [Extend 500 bp upstream]
        │
        ▼
  promoters_500bp.bed  ✓
```

---

## Quick Start

```bash
# 1. Create environment
mamba create -n go_enrichment python=3.12 && mamba activate go_enrichment
mamba install -c bioconda bedtools samtools emboss

# 2. Index genome
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.genome

# 3. Parse TSS coordinates
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"}
NR>1{
    if($8==-1 || $8=="" || $7=="") next
    chrom="chr"$5
    if(chrom=="chrMT") chrom="chrM"
    tss=$8 + 0
    strand = ($6==-1) ? "-" : "+"
    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand
}' > genes_tss_clean.bed

# 4. Drop unrecognized chromosomes
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed

# 5. Extend to promoter windows
bedtools slop -i genes_tss_final.bed -g hg38.genome -l 500 -r 0 -s > promoters_500bp.bed
```

---

## Inputs

| Resource | Source |
|---|---|
| `hg38.fa` | [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) |
| `human_gene_annotation.tsv.gz` | Provided annotation table |

---

## Output Files

| File | Description |
|---|---|
| `genes_tss_final.bed` | Cleaned TSS positions (chromosome-filtered) |
| `promoters_500bp.bed` | Final promoter intervals — use this for downstream analysis |

**BED format (6 columns):**

```
chr1    11868    12369    chr1@11868-11869|DDX11L1    .    +
 [1]     [2]      [3]              [4]               [5]  [6]
```

| # | Field | Value |
|---|---|---|
| 1 | Chromosome | e.g. `chr1`, `chrM` |
| 2 | Start | 500 bp upstream of TSS |
| 3 | End | TSS position |
| 4 | ID | `chr@tss_start-tss_end\|gene_name` |
| 5 | Score | `.` (unused) |
| 6 | Strand | `+` or `-` |

---

## Strand Logic

The `-s` flag in `bedtools slop` makes upstream mean *biologically* upstream:

```
(+) strand gene:    <=====[500bp]=====[TSS]>>>>>>>>>>>>>
(-) strand gene:    <<<<<<<<<<<[TSS]=====[500bp]=====>
```

So for a minus-strand gene, coordinates extend toward *larger* numbers:

```text
# Before slop:
chrM    4400    4401    chrM@4400-4401|MT-TQ    .    -

# After slop -l 500 -s:
chrM    4400    4901    chrM@4400-4401|MT-TQ    .    -
```

---

## Use Cases

The output BED file plugs directly into:

- **Motif analysis** — MEME, HOMER, JASPAR scans
- **GO enrichment** — link promoter features back to gene function
- **Peak overlap** — intersect with ChIP-seq or ATAC-seq data via `bedtools intersect`
