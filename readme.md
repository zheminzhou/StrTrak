# StrTrak:Strain Tracker for metagenomic data

StrTrak consists of three major modules: 
(a) StrBase: builds a genus-wide database
(b) StrProfile: runs species profiling and reconstructs MACGs
(c) StrPhylo: integrates genomic data to build phylogeny

# INSTALLATION 
## Dependencies

StrTrak was developed and tested in Python >=3.8. It depends on several Python libraries: 
~~~~~~~~~~
click
scipy
ete3
numba
numpy
pandas
json
igraph
biopython
pyarrow
fastparquet
~~~~~~~~~~

All libraries can be installed using pip: 

~~~~~~~~~~
pip install click numba numpy pandas biopython ete3 json igraph pyarrow fastparquet
~~~~~~~~~~

StreTrak also calls two 3rd party programs:

~~~~~~~~~~
samtools
iqtree
bindash
minimap2
HMMER
EMBOSS
~~~~~~~~~~

Both can be install via 'apt' in  UBUNTU:
~~~~~~~~~~
sudo apt install -y samtools iqtree bindash minimap2 hmmer emboss
~~~~~~~~~~

The whole environment can also be installed in conda:


~~~~~~~~~~
conda create --name strtrak python==3.11
conda activate strtrak
conda install -c conda-forge click numba numpy pandas biopython ete3 json igraph pyarrow fastparquet
conda install -c bio-conda samtools iqtree bindash minimap2 hmmer emboss
~~~~~~~~~~

The installation process normally finishes in <30 minutes.

## Install StrTrak

~~~~~~~~~~
git clone https://github.com/Zhou-lab-SUDA/StrTrak.git
~~~~~~~~~~

# Quick Start (with examples)
## Build a genus-wide database.

### Pre-built databases download
The table below offers information about the pre-built databases of 8 bacterial genus used in the paper. Users can download these databases and use them to identify strains directly (saved in USCGs_db folder).


Genus   |	Source  | Number of species |	Number of representative genomes
------------ | -------------| ------------- | ------------- 
Aeromonas |  NCBI | 33 | 454  
Elizabethkingia |  NCBI | 9 | 45  
Neisseria |  NCBI | 41 | 240  
Serratia |  NCBI | 27 | 194  
Staphylococcus |  NCBI | 71 | 260  
Stenotrophomonas |  NCBI | 49 | 304  
Streptococcus |  NCBI | 195 | 1557  
Veillonella |  NCBI | 33 | 78  

You can also use ncbi-genome-download tool to download genomes from NCBI

### Or use StrTrak to build your own custom database.

~~~~~~~~~~~
$ cd StrTrak/test
$ ../StrTrak base -i Streptococcus.USCG_test.txt -p Streptococcus.USCG_test
~~~~~~~~~~~

## Use StrTrak to screen species and reconstruct the metagenome-assembled core genotypes (MACGs) in metagenome short reads.

~~~~~~~~~~~
$ ../StrTrak profile -q simulate_metagenomeReads_1.fastq.gz -q simulate_metagenomeReads_2.fastq.gz -d  Streptococcus.USCG_test -o test_metagenomeReads
~~~~~~~~~~~

## Use StrTrak to reconstruct phylogeny based on MACGs

~~~~~~~~~~~
$ ../StrTrak phylo -d Streptococcus.USCG_test -r GCF_002055535 -o Streptococcus.USCG_test Streptococcus.USCG_test.json
~~~~~~~~~~~
## Full command-line options

~~~~~~~~~~
$ ../StrTrak --help

Usage: StrTrak [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  dbase
  phylo
  profile

~~~~~~~~~~

~~~~~~~~~~
$ ../StrTrak dbase --help

Usage: StrTrak dbase [OPTIONS]

Options:
  -i, --fna_list TEXT    comma-delimited list of the included fasta files in a
                         format of "file_name,accession,species".  [required]
  -d, --min_dist FLOAT   minimum distance between representatives. [default:
                         0.01]
  -c, --min_cov FLOAT    minimum presences of genes and genomes. [default:
                         0.75]
  -p, --prefix TEXT      prefix for the output [required].  [required]
  -t, --threads INTEGER  number of threads to use [Default: 20]
  --help                 Show this message and exit.
~~~~~~~~~~

~~~~~~~~~~
$ ../StrTrak profile --help

Usage: StrTrak profile [OPTIONS]

Options:
  -q, --query TEXT           fastq file [specify multiple times for multiple
                             fastq files]  [required]
  -d, --dbname TEXT          name of the USCG database [required]  [required]
  -o, --output TEXT          prefix of the output. [required]  [required]
  -t, --num_threads INTEGER  number of threads [Default: 16]
  -m, --max_dist FLOAT       maximum distance of alignment [Default: 0.08]
  --min_depth INTEGER        minimum number of coverage to call a base.
                             [Default: 3]
  --min_consensus FLOAT      minimum proportion of consensus base to call
                             [Default: 0.8]
  --help                     Show this message and exit.

~~~~~~~~~~

~~~~~~~~~~
$ ../StrTrak phylo --help

Usage: StrTrak phylo [OPTIONS] [USCGS]...

Options:
  -d, --dbname TEXT       absolute path of the database. [required]
                          [required]
  -r, --ref_acc TEXT      accession in the database to be used as reference.
                          [required]  [required]
  --min_identity FLOAT    minumum identity of a sequence comparing to ref_acc.
                          default: 0.94 [0. - 1.]
  --min_presence FLOAT    minumum coverages of genes for a strain to be
                          evaluated. [default: 0.2]
  --risky                 use low depth sites for tree. default: False
  -N, --no_db             do not include genomes in the database
  -g, --genome_list TEXT  additional genomes to be included in the analysis
  -j, --json_list TEXT    list of uscgs as a file [default: None]
  -o, --output TEXT       prefix of the outputs.  [required]
  -n, --n_proc INTEGER    number of processes [default: 8]
  --help                  Show this message and exit.

~~~~~~~~~~

