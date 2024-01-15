## Info


\
hg19.fa (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/):

```txt
>chr1
ATATATATATTATTTTATATATA...
TATTTTAAAGGGGGGGCCC...
>chr2
ATATATATATTATTTTATATATA...
...
```

\
Samples.vcf (https://figshare.com/articles/dataset/Whole_Exome_Data_VCF_files/13696750/1?file=26347630):

```txt
##fileformat=VCFv4.2
##FILTER=<ID=FSFilter,Description="FS > 60.0">
##FILTER=<ID=HaplotypeFilter,Description="HaplotypeScore > 13.0">
##FILTER=<ID=LowCovFilter,Description="DP <= 20">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=MQFilter,Description="MQRankSum < -12.5">
##FILTER=<ID=QDFilter,Description="QD < 2.0">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M46
chr1	13116	rs201725126	T	G	216.77	PASS	AC=1;AF=0.500;AN=2;BaseQRankSum=-3.141;ClippingRankSum=0.000;DB;DP=63;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=32.69;MQRankSum=-5.538;QD=3.44;ReadPosRankSum=0.802;SNPEFF_EFFECT=INTRON;SNPEFF_EXON_ID=2;SNPEFF_FUNCTIONAL_CLASS=NONE;SNPEFF_GENE_BIOTYPE=processed_transcript;SNPEFF_GENE_NAME=DDX11L1;SNPEFF_IMPACT=MODIFIER;SNPEFF_TRANSCRIPT_ID=ENST00000456328;SOR=0.312	GT:AB:AD:DP:GQ:MQ0:PL	0/1:0.780:49,14:62:99:0:245,0,3203
...
```

\
Output fasta:

```txt
...
>chr14:24657396-24657446
CCGAGAGAAGTGGGTAGGTACTTACTCTGTTTCTCTCTGCAGGGCCGTCA
...
>chrY:13311565-13311615
GACAGCCAGGAGAAGGTCGGAGAGGCTCTGCCAATGGCTGGGGGCCCTGC
...
```


## Usage

* Define paths for `fasta`, `vcf` input files and directory where result `fasta` files will be stored.

* Define `consensuses length` and `consensuses count` (integers).

* Define `allele frequency (AF)` (float).

* Run:

```sh
./main.py --help
```

* Example:

```sh
python3 ./main.py \
--reference=/Path/to/fasta \
--vcf=/Path/to/vcf \
--output=/Directory/to/save/result \
--length=100 \
--count=1000000 \
--allele_frequency=0.5
```

## Note

* With count=1000000 and length=100 program will take about 2 minutes and 30 seconds.
