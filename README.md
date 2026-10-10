# Toxo MLST Allele Caller

A free, browser-based tool for calling *Toxoplasma gondii* MLST alleles from Sanger reads. It replaces the manual Geneious Prime workflow and needs no software installation or computational skills.

**Status:** proof of concept. PCR primers are built in for 11 markers. Reference alleles are built in for 7 of them, taken from 15 *T. gondii* genomes in ToxoDB (ARI, BR9, CAST, COUG, DOM2, FOU, GT1, MAS, ME49, P89, PRC2, RH88, RUB, VAND, VEG). PK1 also includes the lab's own Type I, II, III, X and COUG references.

## Reference panel

| Marker | Expected product | Compared region | Distinct alleles among the 15 genomes |
|---|---|---|---|
| PK1 | 902–903 bp | 845 bp | 10 |
| SAG1 | 476 bp | 435 bp | 6 |
| alt.SAG2 | 726–729 bp | 686 bp | 7 |
| 3′-SAG2 | 326–327 bp | 287 bp | 5 |
| SAG3 | 311 bp | 271 bp | 5 |
| BTUB | 411 bp | 371 bp | 5 |
| GRA6 | 344 bp | 306 bp | 9 |
| c22-8, c29-2, L358, Apico | primers only | — | references not yet added |

Every primer pair was checked against all 15 genomes (both sites found, at most one mismatch). Notes:

- **SAG1:** the forward primer sits in a 27 bp tandem repeat upstream of the start codon and can bind more than one copy. The tool uses the innermost pair of sites, which gives the same 476 bp product in every strain.
- **alt.SAG2:** the built-in pair is the outer PCR. The inner (nested) pair from Su et al. 2006 gives 546 bp and is listed as a note. If the nested product is what gets sequenced, switch to the inner pair.
- **PK1 and BTUB** references are ToxoDB transcript or gene sequences without flanks; both primer sites lie inside them, and the products match the published sizes.

When several references are identical over the compared region, the call shows the first name in panel order, and the others are listed as identical ("Type I = 5 more").

## How to use it

1. Open the tool's web page (or download `index.html` and double-click it).
2. Drop in your `.ab1` files, a folder, or a `.zip`. Forward and reverse reads are paired by sample name.
3. Clean samples are called automatically. Samples marked amber or red come with an explanation and the chromatogram at each flagged position.
4. Download the allele table (CSV, opens in Excel), the details and the consensus sequences.

Your sequence files are read inside your own browser and are never uploaded anywhere.

## What it does

| Step | Method |
|---|---|
| Match each read to a marker and orient it | Shared 12-mers with the reference panel |
| Quality trimming (.ab1 only) | Mott's algorithm, error limit 0.05 (the Geneious default) |
| Mixed / heterozygous positions (.ab1 only) | Second peak ≥ 0.33 × main peak, written as an IUPAC code |
| Forward + reverse assembly | Overlap alignment with affine gaps; disagreements are excluded and reported |
| Comparison with references | Alignment to every reference over the region between the marker's PCR primers (or, without primers or with the full-length setting, the region all references cover). References can be amplicons or whole genes |
| Allele call | Zero differences required. Anything else is flagged for review, with an NCBI BLAST link |

## Reference panel format

Load your own references from the **Reference panel** section of the page. The file is FASTA, with headers in this form:

```
>MARKER|Allele name|optional note
>SAG3|Type II|GenBank AB123456
```

Primers go on a header line with no sequence below it:

```
>SAG3|primers|FORWARD_PRIMER|REVERSE_PRIMER
```

## Files

- `index.html`: the complete tool in one file
- `src/core.js`: analysis code (also used by the tests)
- `src/app.html`: page template; `build.py` combines it with `core.js`, the panel and the primers to produce `index.html`
- `src/panel.json`: built-in reference alleles
- `src/primers.json`: built-in PCR primers and their sources
- `add_marker.py`, `check_marker.js`: add a marker's ToxoDB sequences to the panel, then check primer sites, alleles and simulated calls (`python3 add_marker.py GRA6 gra6.fasta 'genomic sequence'`, then `node check_marker.js GRA6` and `python3 build.py`)
- `test/`: Node test runner and a synthetic `.ab1` generator

Developed by VEuPathDB for the *Toxoplasma* research community.
