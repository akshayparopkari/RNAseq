# RNAseq

[![LICENSE](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg?style=plastic)](https://github.com/akshayparopkari/RNAseq/blob/master/LICENSE.md)

**DESCRIPTION**

This bash file that stiches together the preprocessing steps for raw RNA-seq data in FASTQ files from 3′RNA-seq project as well as outputs differential gene expression results.

---

**INSTALLATION**

You can download this GitHub repository using the following command - 

```sh
# navigate to desired location on your machine and run -

git clone https://github.com/akshayparopkari/RNAseq.git
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

- VERSION: 0.0.4
- AUTHOR: Akshay Paropkari
- LICENSE: [BSD 3-Clause License](LICENSE.md)
