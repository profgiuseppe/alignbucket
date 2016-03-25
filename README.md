# Documentation

AlignBucket is a software for splitting a fasta file into smaller pieces suitable for
alignment with BLAST. The constraint used to optimize the result is the required
minimum alignment coverage.

Source code and example files are available also at the
[Bologna Biocomputing Group](http://www.biocomp.unibo.it/~giuseppe/partitioning.html) website.

## Citation
G. Profiti, P. Fariselli, and R. Casadio. ["AlignBucket: a tool to speed up 'all-against-all' protein sequence alignments optimizing length constraints"](http://bioinformatics.oxfordjournals.org/content/31/23/3841). *Bioinformatics*, 31 (23): 3841-3843, 2015.

## Requirements

The program was tested on a GNU/Linux Debian 7 system.

The required libraries are:
- [GNU Multiple Precision library](https://gmplib.org/) (gmp)
- [boost library](http://www.boost.org/) (at least version 1.49)

Also Python version >= 2.6.6 is required for some helper functions.

## Compilation

To compile the program, you need to install the required libraries.
On Debian and Ubuntu you can use the following commands:
```bash
     sudo apt-get install g++
     sudo apt-get install libgmp3-dev
     sudo apt-get install libboost-dev
     sudo apt-get install libboost-program-options-dev
```
Then you must run:
```bash
     make
```   
or, if you don't have Automake in your system:
```bash
     g++ -o alignbucket src/alignbucket.cpp -lgmpxx -lgmp -lboost_program_options
```     
Then, you have to set execution permission on the program:
```bash
     chmod u+x alignbucket
```

##  Usage

To split your fasta file into optimized buckets for 90% coverage, just run the
following command:
```bash
    ./alignbucket.sh <fasta file>
```

Example with Swissprot sequences:
```bash
    ./alignbucket.sh uniprot_sprot.fasta
```

The output will be a set of fasta file named *xxx-yyy.fasta*, with *xxx* and *yyy*
being the minimum and the maximum length of the sequences inside the file.


For different coverage ratios, just run
```bash
    ./alignbucket.sh <fasta file> <percentage>
```

Example, for 75% coverage use:
```
    ./alignbucket.sh uniprot_sprot.fasta 75
```

## Authors

 - [Giuseppe Profiti](https://sites.google.com/site/peppeprofiti)
 - [Piero Fariselli](http://www.biocomp.unibo.it/piero/)

## License

The program is released under a GNU Public License version 2. Please see the
LICENSE.txt file
