---
output: 
  html_document: 
    keep_md: yes
---
# Genetic Analysis Project

## Binary Matrix Conversion and STRUCTURE Preparation

This repository contains two R Markdown notebooks that facilitate the processing and analysis of **multilocus** genetic data, from converting to a binary matrix to preparing data for the program **STRUCTURE**.

---

### 1. `binary.Rmd`

This notebook processes a dataset containing alleles per **locus** for multiple individuals and converts it into a **binary matrix**. Each column in this matrix represents the presence (1) or absence (0) of a specific allele at a given locus.

#### Main features:

- Loads and preprocesses data, replacing missing values coded as `-9` with `NA`.
- Converts data from wide to long format for easier manipulation.
- Creates a unique identifier for each locus-allele combination.
- Removes duplicates and handles potential issues caused by homozygous genotypes.
- Generates a binary matrix ready for downstream analysis.

> **Notebook link:** [binary.Rmd](binary.md)

---

### 2. `structure.Rmd`

This notebook transforms a typical fragment analysis data matrix into the input format required by the program [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html), a widely used software for investigating population structure using multilocus genotype data.

#### Main features:

- Loads multilocus data and checks file integrity, removing empty columns if necessary.
- Replaces missing values with the STRUCTURE-specific code `-9`.
- Separates loci data columns from metadata columns (e.g., sample ID, population).
- Properly reorganizes and encodes data, accounting for polyploid samples (e.g., tetraploids).
- Combines metadata with processed genetic data.
- Generates and saves a STRUCTURE-ready input file.

> **Notebook link:** [structure.Rmd](structure.md)

---

## Recommended workflow

1. Run `binary.Rmd` to convert your original genetic data into a binary presence/absence matrix.
2. Use `structure.Rmd` to prepare and generate the input file compatible with STRUCTURE.
3. Run STRUCTURE externally with the generated files to analyze population genetic structure.

---

## Requirements

- R and RStudio
- R packages: `polysat`, `dplyr`, `tidyr`, `stringr`
- STRUCTURE software installed to run population structure analyses (if applicable)

---

## Contact

For questions or suggestions, please open an issue or contact me directly.

---

Thank you for using this project for your genetic analyses!
