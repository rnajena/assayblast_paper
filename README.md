### Code repository for the AssayBLAST v2 publication
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19738487.svg)](https://doi.org/10.5281/zenodo.19738487)

This repository hosts results from the AssayBLAST v2 publication and code to reproduce these.

Accession numbers can be found in `data/organisms.json` and `assay/microarray_results.txt`.

Prerequisites: conda installation using conda-forge and bioconda channels.

First, please install the necessary packages:

```bash
conda create -n assayblast_paper "python=3.14" "blast=2.17.0" matplotlib numpy pandas seaborn scipy "biopython=1.85"
conda activate assayblast_paper
pip install "rnajena-sugar==1.0" "assay_blast==2.1"
```

#### Create the example result tables (tables 2, 3)

The overview and details tables are located in the example folder. To re-create them switch to the `example` folder and run the following commands.

```bash
assay_test -d .
assay_blast example_database.fasta -q example_queries.fasta --mismatch 2
assay_analyze blast_results.tsv --mismatch 2
```

#### Analyze score to bit score relationship (figure 1)

The first script downloads the genomes.
The second script calculates the bit scores and plots the relationship.

```bash
python prepare_genomes.py
python bitscore_vs_score.py
```

#### Comparison to AssayBLAST v1.0 (table 4)

This section requires the genomes downloaded in the previous section. We download AssayBLAST v1.0 from GitHub and create a modified version with a fixed E-value of 1000. The script `evaluation.py` runs the different AssayBLAST versions with various primers and organisms and generates the comparison table.

```bash py
git clone --branch v1.0 https://github.com/rnajena/AssayBLAST.git AssayBLASTv1
sed 's/evalue=100000/evalue=1000/g' AssayBLASTv1/assayBLAST.py > AssayBLASTv1/assayBLAST_e1k.py
python evaluation.py
```


#### Comparison with Microarray assay (figure 3)

The `microarray_results.txt` and `primer_and_probe.fasta` files were obtained from the supplementary material of Collatz et al. (2025, [doi](https://doi.org/10.3390/applbiosci4020018)).

The following code downloads the corresponding genomes, executes AssayBLAST v2, and generates the confusion matrix.

```
python assay_evaluation_prepare.py
cd assay
assay_blast -q primer_and_probe.fasta genomes.fasta
assay_analyze blast_results.tsv --distance 150
cd ..
python assay_evaluation_analyze.py
```
