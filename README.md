# Bulk RNA-seq Alignment and Analysis

### Pipeline Steps:
1. Runs an initial [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) quality assessment of the reads
2. Trims adapter sequences using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
3. Creates a second QC analysis of the trimmed reads
4. Maps the reads to the reference genome using [STAR](https://github.com/alexdobin/STAR)
   - Good tutorial on running STAR: [Introduction to RNA-seq using high-performance computing](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)
5. Filters out duplicate reads
6. Runs FPKM analysis on the aligned reads using [Cufflinks](https://cole-trapnell-lab.github.io/cufflinks/getting_started/)
7. Creates an integer count matrix using [featureCounts](https://subread.sourceforge.net/featureCounts.html)