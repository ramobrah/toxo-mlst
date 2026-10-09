# Toxo MLST Allele Caller

A free, browser-based tool for calling *Toxoplasma gondii* MLST alleles from Sanger reads. It replaces the manual Geneious Prime workflow and needs no software installation or computational skills.

**Status:** proof of concept. PCR primers are built in for 11 markers (c22-8, c29-2, L358, SAG1, alt.SAG2, 3′-SAG2, PK1, BTUB, GRA6, SAG3, Apico). Reference alleles are so far included only for PK1 (Types I, II, III, X and COUG/TgCgCa1).

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
- `test/`: Node test runner and a synthetic `.ab1` generator

Developed by VEuPathDB for the *Toxoplasma* research community.
