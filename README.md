# Toxo MLST Allele Caller

A free, browser-based tool for calling *Toxoplasma gondii* MLST alleles from Sanger reads. It replaces the manual Geneious Prime workflow and needs no software installation or computational skills.

**Status:** proof of concept. The built-in reference panel covers only the PK1 marker (Types I, II, III, X and COUG/TgCgCa1) and is not yet complete.

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
| Comparison with references | Alignment to every reference over the region all references cover |
| Allele call | Zero differences required. Anything else is flagged for review, with an NCBI BLAST link |

## Reference panel format

Load your own references from the **Reference panel** section of the page. The file is FASTA, with headers in this form:

```
>MARKER|Allele name|optional note
>SAG3|Type II|GenBank AB123456
```

## Files

- `index.html`: the complete tool in one file
- `src/core.js`: analysis code (also used by the tests)
- `src/app.html`: page template; `build.py` combines it with `core.js` and the panel to produce `index.html`
- `src/panel.json`: built-in reference panel
- `test/`: Node test runner and a synthetic `.ab1` generator

Developed by VEuPathDB for the *Toxoplasma* research community.
