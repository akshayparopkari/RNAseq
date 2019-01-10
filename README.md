**DESCRIPTION**

This bash file that stiches together the preprocessing steps for RNA-seq sequences from 3′RNA-seq project.

**USAGE and EXAMPLE**

sh pipeline.sh [-h] inputs
./pipeline.sh [-h] inputs

**OPTIONS**

inputs: Folder containing all raw sequence FASTQ files, this will be updated in the
        next version NOTE: Care should be taken to supply full or relative path to the script file (this file) and input 
        folder.
                 
**REFERENCES**

1. 3′RNA-seq - https://www.rna-seqblog.com/for-model-species-the-3-rna-seq-method-might-more-accurately-detect-differential-expression/
2. BBTools Suite - https://jgi.doe.gov/data-and-tools/bbtools/
3. FastQC - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
4. STAR - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
                 
**SCRIPT INFORMATION**

- VERSION: 0.0.2
- AUTHOR: Akshay Paropkari
- LICENSE: [MIT](LICENSE.md)