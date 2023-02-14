# RNAseq

[![LICENSE](https://img.shields.io/github/license/akshayparopkari/RNAseq)](https://github.com/akshayparopkari/RNAseq/blob/master/LICENSE.md)

If you use this tool, please cite us -

```
Paropkari AD, Bapat PS, Sindi SS, Nobile CJ. A
Computational Workflow for Analysis of 3' Tag-Seq Data. Curr Protoc. 
2023 Feb;3(2):e664. doi: 10.1002/cpz1.664. PMID: 36779816.
```

**DESCRIPTION**

This bash file that stiches together the preprocessing steps for raw RNA-seq data in FASTQ files from 3′RNA-seq project as well as outputs differential gene expression results.

The [Wiki](https://github.com/akshayparopkari/RNAseq/wiki) page provides detailed information about the using scripts in this repository.

---

**MERCED SETUP**

This pipeline is implemented on [MERCED](https://github.com/ucmerced/merced-cluster/wiki) as a `RNA-seq` conda environment. MERCED users can log into MERCED and activate the environment to access the pipeline modules.

```sh
# activate RNA-seq conda evnvironment
source activate RNA-seq
```

---

**LOCAL INSTALLATION**

You can download this GitHub repository using the following command - 

```sh
# Navigate to desired directory to download this folder on your machine

git clone https://github.com/akshayparopkari/RNAseq.git
```

Read more about creating the local Conda environment, required for running the workflow in [section 3](https://github.com/akshayparopkari/RNAseq/wiki/3.-Loading-RNAseq-environment) of the Wiki.

```sh
# Making script files executable

cd RNAseq/

chmod u+x pipeline.sh

chmod u+x format_counts_table.py
```

---

**USAGE and EXAMPLE**

```sh
sh pipeline.sh [-h] inputs
```

*OR*

```sh
./pipeline.sh [-h] inputs
```
---

**OPTIONS**

inputs: Directory containing __all__ raw RNA-seq FASTQ files

---

**REFERENCES**

1. 3′RNA-seq - https://www.rna-seqblog.com/for-model-species-the-3-rna-seq-method-might-more-accurately-detect-differential-expression/
2. BBTools Suite - https://jgi.doe.gov/data-and-tools/bbtools/
3. FastQC - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
4. STAR - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

---

**SCRIPT INFORMATION**

- VERSION: 0.0.6
- AUTHOR: Akshay Paropkari
- LICENSE: [BSD 3-Clause License](LICENSE.md)
