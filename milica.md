
# Required tools

- fastq-dump (sra-toolkit)
- minimap2 (at least version 2.17, I used 2.24)
- samtools
- seqkit
- bedtools
- bedToGenePred & genePredToGtf (UCSC https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads)

# Directories

- genes fasta files are stored in `${DATA}/${SPECIES}/${SPECIES}_gene_fasta_files` (named `${FASTADIR}`in the following)
- FASTQ files are stored in `${DATA}/${SPECIES}/${PLATFORM}/${FASTQ}` (named `${MAPPINGDIR}`in the following) and all files built from a gievn FASTQ file are stored in the same directory

You thus have to define some environment variables.

```bash
DATA="." 
SPECIES="zf" 
FASTQ="SRR8551562" 
PLATFORM="pacbio"
FASTADIR=${DATA}/${SPECIES}/${SPECIES}_gene_fasta_files
FASTA=all_${SPECIES}_genes
MAPPINGDIR=${DATA}/${SPECIES}/${PLATFORM}/${FASTQ}
MAPPING=${FASTA}-vs-${FASTQ}
```

# How to download genes

```bash
SPECIES="hs"; FASTADIR=${DATA}/${SPECIES}/${SPECIES}_gene_fasta_files; mkdir -p ${FASTADIR}; for gid in ENSG00000107643 ENSG00000050748 ENSG00000128641 ENSG00000010810 ENSG00000058404 ENSG00000140416 ENSG00000148737 ENSG00000197122 ENSG00000187122 ENSG00000169855 ENSG00000185008 ENSG00000092421 ENSG00000139220; do curl -H "Content-Type:application/text" "https://rest.ensembl.org/sequence/id/${gid}?content-type=text/x-fasta" > ${FASTADIR}/$gid.fasta; done
```  

```bash
SPECIES="mm"; FASTADIR=${DATA}/${SPECIES}/${SPECIES}_gene_fasta_files; mkdir -p ${FASTADIR}; for gid in ENSMUSG00000021936 ENSMUSG00000020366 ENSMUSG00000018417 ENSMUSG00000019843 ENSMUSG00000057897 ENSMUSG00000032366 ENSMUSG00000024985 ENSMUSG00000027646 ENSMUSG00000025020 ENSMUSG00000022883 ENSMUSG00000052516 ENSMUSG00000019647 ENSMUSG00000053825; do curl -H "Content-Type:application/text" "https://rest.ensembl.org/sequence/id/${gid}?content-type=text/x-fasta" > ${FASTADIR}/${gid}.fasta; done
``` 

```bash
SPECIES="zf"; FASTADIR=${DATA}/${SPECIES}/${SPECIES}_gene_fasta_files; mkdir -p ${FASTADIR}; mkdir -p ${FASTADIR}; for gid in ENSTGUG00000008103 ENSTGUG00000023050 ENSTGUG00000010724 ENSTGUG00000011998 ENSTGUG00000003464 ENSTGUG00000005002 ENSTGUG00000010805 ENSTGUG00000004937 ENSTGUG00000009734 ENSTGUG00000013537 ENSTGUG00000013538 ENSTGUG00000017547 ENSTGUG00000007778; do curl -H "Content-Type:application/text" "https://rest.ensembl.org/sequence/id/${gid}?content-type=text/x-fasta" > ${FASTADIR}/${gid}.fasta; done
``` 

Don't forget to redefine correctly `SPECIES` and `FASTADIR` afterwards.

# How to download FASTQ

```bash
mkdir -p ${MAPPINGDIR}; fastq-dump --gzip --skip-technical --readids --dumpbase --clip --outdir ${MAPPINGDIR} ${FASTQ}
```

# Mapping genes to FASTQ

We assume that all the 13 genes will be mapped at once. Then we concatenate all the genes in a single mutlifasta file, and we trun minimap of all the genes against the FASTQ file.

Here are the command lines for zebrafinch.

```bash
cat ${FASTADIR}/ENS*.fasta > ${FASTADIR}/${FASTA}.fasta
```

```bash
minimap2 -a -x splice:hq -u f --sam-hit-only ${FASTADIR}/${FASTA}.fasta ${MAPPINGDIR}/${FASTQ}.fastq.gz > ${MAPPINGDIR}/${FASTA}-vs-${FASTQ}.sam
```

Then we sort and index the sam file to have a sorted bam file.

```bash
cat ${MAPPINGDIR}/${MAPPING}.sam | samtools view -bS | samtools sort -o ${MAPPINGDIR}/${MAPPING}.bam ; samtools index ${MAPPINGDIR}/${MAPPING}.bam
```

# Getting stats

## Number of mapped reads

```bash
samtools flagstats ${MAPPINGDIR}/${MAPPING}.bam
```

## Number of mapped reads per gene

```bash
for i in `head -30 ${MAPPINGDIR}/${MAPPING}.sam | grep "@SQ" | cut -d":" -f2 | cut -d$'\t' -f1`; do echo -n $i; samtools view ${MAPPINGDIR}/${MAPPING}.bam $i | wc -l; done
```

# Extracting reads for every gene

We start by creating files containg read IDs.

```bash
for i in `head -30 ${MAPPINGDIR}/${MAPPING}.sam | grep "@SQ" | cut -d":" -f2 | cut -d$'\t' -f1`; do samtools view ${MAPPINGDIR}/${MAPPING}.bam $i | cut -f1 > ${MAPPINGDIR}/${FASTA}-vs-${FASTQ}-${i}.txt; done
```

Then we can extract reads (the method is not really efficient ...).

```bash
for i in `ls ${MAPPINGDIR}/${MAPPING}-*.txt`; do basename="${i%.*}"; seqkit grep -f $i -o ${basename}.fastq.gz ${MAPPINGDIR}/${FASTQ}.fastq.gz; done
```

Il you also want the fasta files of the reads.

```bash
for i in `ls ${MAPPINGDIR}/${MAPPING}-*.fastq.gz`; do basename="${i%.*}"; seqkit fq2fa $i -o ${basename}.fasta; done
```

We also extract a BAM file specific to each gene

```bash
for i in `head -30 ${MAPPINGDIR}/${MAPPING}.sam | grep "@SQ" | cut -d":" -f2 | cut -d$'\t' -f1`; do samtools view -b ${MAPPINGDIR}/${MAPPING}.bam $i > ${MAPPINGDIR}/${FASTA}-vs-${FASTQ}-${i}.bam; done
```

# Getting a BED or GTF file

First we build BED files from BAM files.

```bash
for i in `ls ${MAPPINGDIR}/${FASTA}-vs-${FASTQ}-*.bam`; do basename="${i%.*}"; bedtools bamtobed -bed12 -i $i > $basename.bed; done
```

Then we build genPred files from BED files.

```bash
for i in `ls ${MAPPINGDIR}/${FASTA}-vs-${FASTQ}-*.bed`; do basename="${i%.*}"; bedToGenePred $i $basename.pred; done
```

And finally, we build GTF files from genePred files.

```bash
for i in `ls ${MAPPINGDIR}/${FASTA}-vs-${FASTQ}-*.pred`; do basename="${i%.*}"; genePredToGtf "file" $i  $basename.gtf; done
```

# Example for zebrafinch

Here is an example for the dataset SRR8551562.

## Download of genes
After download of genes, the DATA directory contains 


```
.
└── zf
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        └── ENSTGUG00000023050.fasta

3 directories, 13 files
```


## Download of FASTQ
After download of the FASTQ file (for this dataset, 1.4GBytes, it took about 40 minutes)

```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       └── SRR8551562.fastq.gz
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        └── ENSTGUG00000023050.fasta

5 directories, 14 files
```

## Mapping

After making the multifasta file


```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       └── SRR8551562.fastq.gz
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 15 files
```

After running minimap2 (it took about 8 minutes)

```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       ├── SRR8551562.fastq.gz
    │       └── all_zf_genes-vs-SRR8551562.sam
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 16 files
```

After sort and index.


```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       ├── SRR8551562.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562.bam
    │       ├── all_zf_genes-vs-SRR8551562.bam.bai
    │       └── all_zf_genes-vs-SRR8551562.sam
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 18 files
```


## Stats

### Mapped reads
```
3746 + 0 in total (QC-passed reads + QC-failed reads)
3243 + 0 primary
266 + 0 secondary
237 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
3746 + 0 mapped (100.00% : N/A)
3243 + 0 primary mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Mapped reads per gene
```
ENSTGUG00000003464.2    1567
ENSTGUG00000004937.2     117
ENSTGUG00000005002.2      45
ENSTGUG00000007778.2     196
ENSTGUG00000008103.2      26
ENSTGUG00000009734.2     245
ENSTGUG00000010724.2     105
ENSTGUG00000010805.2       4
ENSTGUG00000011998.2      73
ENSTGUG00000013537.2     225
ENSTGUG00000013538.2    1008
ENSTGUG00000017547.2      69
ENSTGUG00000023050.1      66
```

## Extracting reads

### Reads IDs for each gene

```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       ├── SRR8551562.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.txt
    │       ├── all_zf_genes-vs-SRR8551562.bam
    │       ├── all_zf_genes-vs-SRR8551562.bam.bai
    │       └── all_zf_genes-vs-SRR8551562.sam
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 31 files
```

### FASTQ of the reads for each gene (it took about 11 minutes)

Output:
```
[INFO] 1528 patterns loaded from file
[INFO] 112 patterns loaded from file
[INFO] 44 patterns loaded from file
[INFO] 193 patterns loaded from file
[INFO] 25 patterns loaded from file
[INFO] 231 patterns loaded from file
[INFO] 102 patterns loaded from file
[INFO] 4 patterns loaded from file
[INFO] 73 patterns loaded from file
[INFO] 222 patterns loaded from file
[INFO] 773 patterns loaded from file
[INFO] 68 patterns loaded from file
[INFO] 66 patterns loaded from file
```


```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       ├── SRR8551562.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.txt
    │       ├── all_zf_genes-vs-SRR8551562.bam
    │       ├── all_zf_genes-vs-SRR8551562.bam.bai
    │       └── all_zf_genes-vs-SRR8551562.sam
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 44 files
```

### Transforming fastq.gz into fasta


```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       ├── SRR8551562.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.txt
    │       ├── all_zf_genes-vs-SRR8551562.bam
    │       ├── all_zf_genes-vs-SRR8551562.bam.bai
    │       └── all_zf_genes-vs-SRR8551562.sam
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 57 files
```

### BAM for each gene

```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       ├── SRR8551562.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.txt
    │       ├── all_zf_genes-vs-SRR8551562.bam
    │       ├── all_zf_genes-vs-SRR8551562.bam.bai
    │       └── all_zf_genes-vs-SRR8551562.sam
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 70 files
```

# BED and GTF files

```
.
└── zf
    ├── pacbio
    │   └── SRR8551562
    │       ├── SRR8551562.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000003464.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000004937.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000005002.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000007778.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000008103.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000009734.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010724.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000010805.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000011998.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013537.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000013538.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000017547.2.txt
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.bam
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.bed
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.fastq.fasta
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.fastq.gz
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.gtf
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.pred
    │       ├── all_zf_genes-vs-SRR8551562-ENSTGUG00000023050.1.txt
    │       ├── all_zf_genes-vs-SRR8551562.bam
    │       ├── all_zf_genes-vs-SRR8551562.bam.bai
    │       └── all_zf_genes-vs-SRR8551562.sam
    └── zf_gene_fasta_files
        ├── ENSTGUG00000003464.fasta
        ├── ENSTGUG00000004937.fasta
        ├── ENSTGUG00000005002.fasta
        ├── ENSTGUG00000007778.fasta
        ├── ENSTGUG00000008103.fasta
        ├── ENSTGUG00000009734.fasta
        ├── ENSTGUG00000010724.fasta
        ├── ENSTGUG00000010805.fasta
        ├── ENSTGUG00000011998.fasta
        ├── ENSTGUG00000013537.fasta
        ├── ENSTGUG00000013538.fasta
        ├── ENSTGUG00000017547.fasta
        ├── ENSTGUG00000023050.fasta
        └── all_zf_genes.fasta

5 directories, 109 files
```

