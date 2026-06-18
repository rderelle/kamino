[![Cargo Build & Test](https://github.com/rderelle/kamino/actions/workflows/ci.yml/badge.svg)](https://github.com/rderelle/kamino/actions/workflows/ci.yml)
[![Clippy check](https://github.com/rderelle/kamino/actions/workflows/clippy.yml/badge.svg)](https://github.com/rderelle/kamino/actions/workflows/clippy.yml)
[![codecov](https://codecov.io/github/rderelle/kamino/graph/badge.svg?token=6B8WIGZL2F)](https://codecov.io/github/rderelle/kamino)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/kamino/README.html)

<br><br>
<p align="center">
  <img src="logo_kamino.svg" alt="kamino logo" width="400">
</p>

<br><br>

Builds an amino acid alignment in a reference-free, alignment-free manner from a set of proteomes. Not ‘better’ than traditional marker-based pipelines, but simpler and much faster to run. Typical usages range from between-strains to within-family (prokaryotes) or within-phylum (eukaryotes) phylogenetic analyses.

## Installation
You can either compile the code locally using rustc, or install a precompiled binary from Bioconda:

```bash
conda install bioconda::kamino
```

## Usage
### input
Input consists of proteome files in FASTA format (gzipped or not), with one file per sample. Files can be placed in a single directory, specified with the `-i` argument. Files are then recognised by their extension (.fa, .fas, .fasta, .faa, .fna; gzipped ot not) and filenames minus extensions become isolate names in the final alignment. 
```bash
kamino -i <input_dir> -t 4
```
Alternatively , the path of input files can be provided in a headless tab-delimited file using the `-I` argument. This is useful when file names do not encode the sample name or when proteomes are located in multiple directories. 
```bash
kamino -I <tabular_file> -t 4
```
The first column of the tab-delimited file should contain the isolate names in the final alignment and the second column the path of proteome files.
```bash
species_soil    new_assemblies/proteins.fas.gz
ncbi_species_1    ncbi/GCA_00076868/proteome.fa.gz
ncbi_species_2    ncbi/GCA_01092864/proteome.fa.gz
```
For **bacterial** isolates, the phylogenomic alignment can also be generated directly from genome assemblies by selecting the `--genomes` argument, using either `-i` or `-I`. In this case, an ultra-fast but approximate protein prediction is performed, with predicted proteomes stored in a temporary directory.

```bash
kamino -i <input_dir> -t 4 --genomes
```
### output
kamino outputs a FASTA amino acid alignment, a file containing the percentage of missing data per isolate, and a file containing variant group coordinates in the FASTA alignment. By default, these output files use the prefix 'kamino_', but this can be changed with the `-o` argument.

Additionally, a Neighbour-joinging tree can be generated together with other output files by selecting the `--NJ` optional argument.

```bash
kamino -i <input_dir> -t 4 --NJ
```

## Documentation
More information are available at https://docs.rs/kamino-cli/latest/kamino_cli/

Please let me know if anything is unclear or missing, and I'll update the doc accordingly.

## Citation
If you use kamino, please cite:

> Romain Derelle, John A. Lees and Leonid Chindelevitch. 2026.
> kamino: proteome-wide variant calling for amino acid phylogenomics.
> https://www.biorxiv.org/content/10.64898/2026.05.21.726148v1

---
This codebase is provided under the MIT License. Some parts of the code were drafted using AI assistance.
