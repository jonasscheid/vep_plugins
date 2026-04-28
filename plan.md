# AlphaGenome VEP Plugin — Implementation Plan

## Overview

Build `AlphaGenome.pm`, a VEP plugin that annotates variants with pre-computed AlphaGenome variant-effect scores across all 11 modalities. The plugin follows the established `BaseVepTabixPlugin` pattern (like SpliceAI.pm, Enformer.pm) and reads from a tabix-indexed TSV file containing scores exported from the AlphaGenome Python SDK via `tidy_scores()`.

Peptide extraction from cryptic splice junctions is **not** part of this plugin — that is a downstream process (AGFUSION, GTF/FASTA translation) that consumes the splice junction annotations this plugin produces.

---

## AlphaGenome Modalities (verified from API)

The AlphaGenome SDK defines 11 `OutputType` values, each corresponding to a modality:

| OutputType | Description | Scorer aggregation |
|------------|-------------|--------------------|
| `RNA_SEQ` | RNA-seq gene expression tracks | log2 fold change of summed signals |
| `CAGE` | CAGE expression at TSS | log2 fold change of summed signals |
| `PROCAP` | PRO-cap transcription initiation | log2 fold change of summed signals |
| `DNASE` | DNase I chromatin accessibility | log2 fold change of summed signals |
| `ATAC` | ATAC-seq chromatin accessibility | log2 fold change of summed signals |
| `CHIP_HISTONE` | ChIP-seq histone modifications | log2 fold change of summed signals |
| `CHIP_TF` | ChIP-seq transcription factor binding | log2 fold change of summed signals |
| `SPLICE_SITES` | Splice site donor/acceptor probabilities | max abs difference (0–1 scale) |
| `SPLICE_SITE_USAGE` | Splice site usage fractions | max abs difference (0–1 scale) |
| `SPLICE_JUNCTIONS` | Splice junction counts (donor–acceptor pairs) | max abs log fold change |
| `CONTACT_MAPS` | 3D chromatin interaction maps | *(not gene-centric)* |

All 11 should be supported. `CONTACT_MAPS` is the one modality that is not gene-centric (returns a single global score per variant, not per-gene), so it needs slightly different handling.

---

## Data Format

The plugin reads a **bgzipped, tabix-indexed TSV** exported from AlphaGenome's `tidy_scores()` output. This aligns directly with the SDK's recommended export format.

### TSV columns

```
#CHROM  POS  REF  ALT  GENE_ID  GENE_NAME  OUTPUT_TYPE  TRACK_NAME  RAW_SCORE  QUANTILE_SCORE
```

| Column | Description |
|--------|-------------|
| `CHROM` | Chromosome (matching VEP input, no `chr` prefix by default) |
| `POS` | 1-based variant position |
| `REF` | Reference allele |
| `ALT` | Alternate allele |
| `GENE_ID` | Ensembl gene ID (e.g. `ENSG00000100342`), or `.` for non-gene-centric scorers |
| `GENE_NAME` | HGNC gene symbol (e.g. `APOL1`), or `.` |
| `OUTPUT_TYPE` | One of the 11 modality names above |
| `TRACK_NAME` | Specific track/tissue (e.g. `UBERON:0036149 total RNA-seq`) |
| `RAW_SCORE` | Raw variant effect score from the scorer |
| `QUANTILE_SCORE` | Quantile-normalized score (0–1 unsigned, or -1 to 1 signed) against common variant background |

This is a **one-row-per-(variant, gene, modality, track)** long format, matching `tidy_scores()` output.

### Why this format

- Directly exportable from `tidy_scores()` → pandas → TSV → bgzip → tabix
- No custom post-processing needed beyond sort + compress + index
- Keeps the VEP plugin as a pure lookup — all ML inference happens upstream
- Both `raw_score` and `quantile_score` are included so users can filter on either

### Data preparation (documented in plugin POD)

```bash
# 1. Run AlphaGenome batch variant scoring (Python SDK)
# 2. Export via tidy_scores() to TSV
# 3. Sort, compress, index:
sort -k1,1 -k2,2n alphageome_scores.tsv | bgzip -c > alphageome_scores.tsv.gz
tabix -s 1 -b 2 -e 2 alphageome_scores.tsv.gz
```

---

## Implementation Steps

### Step 1: Create `AlphaGenome.pm` scaffold

```perl
package AlphaGenome;
use strict;
use warnings;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);
```

### Step 2: Implement `new()` constructor

Parameters:
- `file` (required) — path to tabix-indexed AlphaGenome TSV
- `cutoff` (optional, default 0) — minimum absolute raw score to report (0 = no filtering)
- `quantile_cutoff` (optional, default 0) — minimum quantile score to report
- `modalities` (optional, default all) — colon-separated list of modalities to include (e.g. `SPLICE_JUNCTIONS:RNA_SEQ:SPLICE_SITES`)

```
--plugin AlphaGenome,file=/path/to/scores.tsv.gz
--plugin AlphaGenome,file=/path/to/scores.tsv.gz,cutoff=0.5,modalities=SPLICE_JUNCTIONS:SPLICE_SITES:SPLICE_SITE_USAGE
```

Logic:
- Parse params via `params_to_hash()`
- Validate file exists, add via `$self->add_file()`
- `expand_left(0)`, `expand_right(0)` (exact position lookup)
- Parse modalities into a hash set for O(1) lookup
- Store cutoff values on `$self`

### Step 3: Implement `feature_types()`

Return `['Feature', 'Intergenic']` — like Enformer, this operates at the variant-feature level (not transcript-specific), because AlphaGenome scores are per-gene not per-transcript. The plugin matches results to genes via `GENE_NAME`/`GENE_ID`.

### Step 4: Implement `get_header_info()`

```
AlphaGenome_GENE_ID        — Ensembl gene ID
AlphaGenome_GENE_NAME      — Gene symbol
AlphaGenome_OUTPUT_TYPE    — AlphaGenome modality
AlphaGenome_TRACK_NAME     — Specific track/tissue
AlphaGenome_RAW_SCORE      — Raw variant effect score
AlphaGenome_QUANTILE_SCORE — Quantile-normalized score
```

For JSON/REST output, nest under `{AlphaGenome => [...]}` (array of result hashes, one per passing modality/track).

For tab/VCF output, return the **top-scoring** result across requested modalities as a pipe-delimited string: `GENE_NAME|OUTPUT_TYPE|TRACK_NAME|RAW_SCORE|QUANTILE_SCORE`.

### Step 5: Implement `parse_data()`

```perl
sub parse_data {
  my ($self, $line) = @_;
  chomp $line;
  my @f = split /\t/, $line;
  return {
    chr    => $f[0],
    start  => $f[1],
    ref    => $f[2],
    alt    => $f[3],
    result => {
      gene_id        => $f[4],
      gene_name      => $f[5],
      output_type    => $f[6],
      track_name     => $f[7],
      raw_score      => $f[8],
      quantile_score => $f[9],
    }
  };
}
```

### Step 6: Implement `run()`

1. Extract variant info from `$tva` (VariationFeature, chr, start, end, alleles)
2. Query `$self->get_data($chr, $start, $end)`
3. For each row:
   a. Match alleles via `get_matched_variant_alleles()`
   b. Apply modality filter (skip if `output_type` not in requested set)
   c. Apply cutoff filters (skip if `abs(raw_score) < cutoff` or `abs(quantile_score) < quantile_cutoff`)
   d. Collect passing results
4. Output formatting:
   - **JSON**: return `{AlphaGenome => \@results}` (array of all passing result hashes)
   - **Tab/VCF**: select the result with the highest `abs(raw_score)` and return as pipe-delimited string

### Step 7: Implement `get_start()` / `get_end()`

```perl
sub get_start { return $_[1]->{start}; }
sub get_end   { return $_[1]->{end}; }
```

### Step 8: Write POD documentation

- NAME, SYNOPSIS, DESCRIPTION
- All parameters with defaults
- Output format per VEP mode (tab, JSON, VCF)
- Data preparation instructions (AlphaGenome SDK → tidy_scores → TSV → bgzip → tabix)
- Citation: Avsec et al., Nature 2026, https://www.nature.com/articles/s41586-025-10014-0
- Links to AlphaGenome SDK: https://github.com/google-deepmind/alphagenome

---

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `AlphaGenome.pm` | **Create** | New VEP plugin (~200–250 lines) |

---

## Design Decisions

1. **Pre-computed scores, not live inference**: AlphaGenome requires 1 Mb context and GPU. The plugin reads pre-computed tabix-indexed scores, matching the convention of every other deep-learning VEP plugin (SpliceAI, Enformer, AlphaMissense).

2. **No peptide extraction in the plugin**: Peptide extraction from cryptic splice junctions is a downstream process. This plugin annotates variants with AlphaGenome's splice junction scores (and all other modalities). Downstream tools (AGFUSION, custom GTF/FASTA translators) consume these annotations to extract cryptic peptides. This avoids duplicating translation logic and keeps the plugin minimal.

3. **TSV mirrors `tidy_scores()` output**: The data format maps 1:1 to the AlphaGenome SDK's recommended export. Users run batch scoring → `tidy_scores()` → write TSV → bgzip + tabix. No custom data reshaping needed.

4. **Both raw and quantile scores**: The SDK produces both. Raw scores are useful for modality-specific thresholds (e.g. splice site probabilities are 0–1). Quantile scores enable cross-modality comparison against a common variant background.

5. **All 11 modalities in one plugin**: Single plugin with modality filtering via the `modalities` parameter. Users doing splice analysis pass `modalities=SPLICE_JUNCTIONS:SPLICE_SITES:SPLICE_SITE_USAGE`. Users wanting everything omit the parameter.

6. **`CONTACT_MAPS` handling**: This modality is not gene-centric (single global score per variant). The plugin handles this by allowing `GENE_ID` and `GENE_NAME` to be `.` — non-gene-centric results are returned regardless of gene context.

7. **Feature-level, not Transcript-level**: AlphaGenome scores are per-gene, not per-transcript. Using `feature_types => ['Feature', 'Intergenic']` matches variants at the position level (like Enformer), then results are matched by gene symbol where applicable.
