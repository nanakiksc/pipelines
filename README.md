# Variant calling pipeline

The following pipeline for variant calling is based on the [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/). It starts from raw sequencing data in FASTQ format and produces a table containing the annotated variants.

Please note that the scripts contain local paths to the data resources files, so they must be edited to point to the files on your system.

## 1. Data pre-processing

The first part of the pipeline maps the raw reads to a reference genome, then marks duplicate reads, performs base quality score recalibration, sorts the BAM file and creates its index.

### Dependencies

Software:
- BWA
- GATK
- SAMtools

Data:
- FASTA sequence of the reference genome (hg38)
- BWA index files for the reference genome
- dbSNP database
- Mills and 1000G gold standard indel annotations
- Known indels annotation for reference genome

### Usage
For paired-end sequencing reads use:
```shell
$ ./data_pre_pe.sh <sample_name>
```

For single-end sequencing reads use:
```shell
$ ./data_pre_se.sh <sample_name>
```

Where `<sample_name>` is the base name of the FASTQ files (without the extension). In the case of paired-end files, the script assumes that the base name is suffixed with _1 and _2 for the files containing the reads 1 and reads 2, respectively, therefore **the suffixes _1 and _2 must not be included in the base name**.

This will produce the alignments file `<sample_name>.bam` and its corresponding BAI index.

## 2. Variant calling

The second part of the pipeline reads the alignment files, then calls variants, computes different metrics used to filter those variants and finally performs the filtering proper. Since most experiments done at the lab are at exome level, the analyses of this second part are limited to the exome coordinates to save time. This feature can be disabled, or modified to observe custom panel coordinates, by editing the scripts themeselves (`--intervals` option of the different GATK tools).

### Dependencies

Software:
- GATK

Data:
- FASTA sequence of the reference genome (hg38)
- dbSNP database
- Mills and 1000G gold standard indel annotations
- Known indels annotation for reference genome
- Exome coordinates
- Database of known germline variants
- Database of known germline variants (only biallelic sites)

### Usage

For somatic variant calling of a tumor sample with its matched normal sample:
```shell
$ ./vc_mn.sh <normal_sample_name> <tumor_sample_name>
```

For somatic variant calling of a tumor sample without a matched normal sample:
```shell
$ ./vc_to.sh <tumor_sample_name>
```

For germline variant calling:
```shell
$ conda activate gatk # To use a Python3 environment, if needed.
$ ./vc_germline.sh <sample_name>
```

Where `<(normal_|tumor_)?sample_name>` are the base names of the BAM files (without the extension) of the normal and tumor samples, respectively. In the case of germline variant calling, any type of sample (normal or tumor) can be used, but only germline variants will be identified.

This will produce the filtered variants file `<sample_name>.vcf.gz` and its corresponding tabix index.

## 3. Annotation

The last part of the pipeline reads in the variants and annotates them against a database of known variants, performs a functional consequence annotation of coding variants, and adds gene-level information. The output is an (Excel-compatible) tab-separated table suitable for manual curation.

### Dependencies

Software:
- BCFtools
- SAMtools
- ANNOVAR

Data:
- dbNSFP4.0a
- ANNOVAR humandb (with non-canonical transcripts excluded)
- COSMIC Cancer Gene Census

### Usage

To annotate somatic variants:
```shell
$ ./anno.sh <sample_name>
```

To annotate germline variants:
```shell
$ ./anno_germline.sh <sample_name>
```

Where `<sample_name>` is the base names of the VCF files (without the extension).

This will produce the annotated variants table `<sample_name>.xls` along with the annotated variants file `<sample_name>.anno.vcf.gz` and its corresponding tabix index.
