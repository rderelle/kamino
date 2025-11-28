[![Cargo Build & Test](https://github.com/rderelle/kamino/actions/workflows/ci.yml/badge.svg)](https://github.com/rderelle/kamino/actions/workflows/ci.yml)
[![Clippy check](https://github.com/rderelle/kamino/actions/workflows/clippy.yml/badge.svg)](https://github.com/rderelle/kamino/actions/workflows/clippy.yml)
[![codecov](https://codecov.io/gh/rderelle/kamino/branch/main/graph/badge.svg?token=YOUR_TOKEN_HERE)](https://codecov.io/gh/rderelle/kamino)

<br>
<h1 align="center">kamino</h1>
<br><br>

<h2> 🚧  work in progress </h2>

<br><br>

From the Spanish word for *path*.  

Builds an amino-acid alignment in a reference-free, alignment-free manner from a set of proteomes.  
Not ‘better’ than traditional marker-based pipelines, but simpler and faster to run.  
  
Typical usages range from between-strains to within-phylum phylogenetic analyses (bacteria, archaea and eukaryotes).

<br>

---
## under the hood
kamino performs the following successive steps:
- lists proteome files from the input directory (-i)
- recodes proteins with the Dayhoff6 recoding scheme
- simplifies proteomes by discarding within-proteome branching k-mers 
- builds a global assembly graph and identifies variant groups as described <a href="https://academic.oup.com/mbe/article/42/4/msaf077/8103706">here</a> (-d)
- converts variant group paths back to amino acids using a sliding window
- filters variant groups by missing data and middle-length thresholds (-m and -l)
- extracts middle positions and incorporate 'constant' positions (-c)
- outputs the final amino acid alignment (-o)

---

## installation

```bash
git clone https://github.com/rderelle/kamino.git
cd kamino
cargo build --release
```

This produces the executable at target/release/kamino

---

## running kamino

Input is a directory containing proteome files in FASTA format, gzipped or not.
The files are detected based on their extension (.faa, .fa, .fas or .fasta) and filenames minus extension become sequence names in the final amino acid alignment.  

Basic run using 4 threads:

```bash
kamino -i <input_dir> -t 4
```
---

## examples

All analyses were performed on a MacBook M4 pro using 4 threads (other parameters set to default):  

| dataset                     | taxonomic diversity  | runtime (min) | memory (GB) | alignment size (aa) |
|-----------------------------|----------------------|---------------|-------------|---------------------|
| 50 *Mycobacterium*          | within-genera        | 0.2           | 3           | 15.618              |
| 400 *Mycobacterium*         | within-genera        | 1.4           | 11          | 12.291              |
| 50 Polysporales (fungi)     | within-order         | 0.8           | 11          | 18.318              |
| 290 *A. fumigatus* (fungi)  | within-species       | 2.7           | 14          | 64.602              |
| 46 *Drosophila*             | within-genera        | 1             | 10          | 190.175             |
| 55 Mammalia                 | within-class         | 1.9           | 12          | 277.334             |  


And using the 400 *Mycobacterium* dataset to understand how parameters can impact the analyses (all using 4 threads):

| parameters                | runtime (min) | memory (GB) | alignment size (aa) |
|---------------------------|---------------|-------------|---------------------|
| [default: k=14, d=4]      | 1.4           | 11          | 12.291              |
| --k 15                    | 1.4           | 13          | 18.772              |
| --depth 5                 | 1.6           | 11          | 16.773              |
| --k 15 --depth 5          | 1.6           | 13          | 24.835              |


---

## FAQ

- **When not to use kamino?**
    * low diversity datasets (ie, within-strain), for which genome-based approaches will be more powerful 
    * very large datasets (eg, thousands of bacterial proteomes or hundreds of vertebrate proteomes)
    * very divergent datasets (eg, animal kingdom)
    * distant outgroup composed of a few isolates: these might have disproportionately more missing data
    * list to be completed ...

- **Is the output reproducible?**
<p>Yes, kamino is fully deterministic so will produce the exact same alignment for a given set of parameters and input proteomes.</p>

- **Where do the constant positions come from?**
<p>They are taken from the left flank of the end k-mer in each variant group, next to the middle positions. Because these positions are recoded under Dayhoff6, some may become polymorphic once converted back to amino acids. Using the default c=3, constant positions represent 60 to 75% of the alignment in the analyses mentioned above.</p>

- **How to get more phylogenetic positions?**
    * increase the k-mer size (-k), but can substantially raise memory usage
    * increase the maximum depth of the graph traversal (-d), but increases the runtime
    * lower the minimum proportion of isolates per position (-m), if that is acceptable for downstream analyses

- **What do '-' and 'X' correspond to in the alignment?**
    * '-' = *absent amino acid*. Created when an isolate is absent from a variant group
    * 'X' = *unknown amino acid*. Created by a strict consensus in two possible cases: when recoded paths are converted back to amino acids, and when an isolate has two or more paths in a given variant group
 

---

This codebase is provided under the MIT License. Some parts of the code were drafted using AI assistance.
