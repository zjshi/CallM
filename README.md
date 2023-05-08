# CallM

CallM for efficient identification of microbial common SNPs  

## Important note

CallM now has been replaced by Maast (https://github.com/zjshi/Maast)  

## What CallM does

Recent spikes in available whole-genome sequences have greatly expanded intra-species diversity especially for prevalent species. However, the increasing availability also introduced high-level of redundancy which imposed computing burden for core-genome SNP calling. The issue will exacerbate as the trend is irreversible. CallM presented here is a tool free of read alignment and assembly, which detects subspecies structure in conspecific genomes and picks tag and reference genomes for rapid and sensitive core-genome SNP calling in both sequencing reads and whole genomes. CallM runs orders of magnitude faster with less RAM use and recovers more core-genome SNPs comparing to other the-state-of-art tools.

## How to cite

The publication of CallM is in preparation. Please cite this GitHub repo as alternative for now. 

## Installation

<b>Install python libraries</b>

* numpy
* Biopython
* pysam
* PyVCF
* ujson
* psutil
* pandas

<b>Install external programs</b>

* [Mash](https://github.com/marbl/Mash) (>= v2.2)
* [MUMmer4](https://github.com/mummer4/mummer) (>= v4.0.0)

## Tutorial to start

`./CallM -h`  

## Examples

### Call common SNPs from a set of whole genomes (default)

`CallM genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/CallM/output/`  

Note:  
By defaulty, CallM first purges redundancy in the input genomes and then call common SNPs from a subset of tag genomes. It also automatically identifies a centroid-genome and use it for the representative genome.

### Call common SNPs from a set of whole genomes without redundancy reduction

`CallM genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/CallM/output/ --skip-centroid --keep-redundancy`  

### Call common SNPs with customized Min. prevalence and MAF thresholds

`CallM genomes --fna-dir /path/to/genomes/ --rep-fna /path/to/rep_genome.fna --out-dir /path/CallM/output/ --min-prev 0.95 --snp-freq 0.001`  
