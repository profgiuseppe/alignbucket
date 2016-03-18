Documentation
------------------

AlignBucket is a software for splitting a fasta file into smaller pieces suitable for
alignment with BLAST. The constraint used to optimize the result is the required
minimum alignment coverage.


Requirements
------------------
The program was tested on a GNU/Linux Debian 7 system.

The required libraries are:
- GNU Multiple Precision library (gmp)
- boost (at least version 1.49)

Also Python version >= 2.6.6 is required for some helper functions.

Archive content
------------------
The archive alignbucket.zip contains:
  - src/alignbucket.cpp (source code)
  - bucketize.py (helper functions)
  - alignbucket.sh (main executable script)
  - Makefile
  - LICENSE.txt
  - README.txt (this file)


Compilation
------------------

To compile the program, you need to install the required libraries.
On Debian and Ubuntu you can use the following commands:

     sudo apt-get install g++
     sudo apt-get install libgmp3-dev
     sudo apt-get install libboost-dev
     sudo apt-get install libboost-program-options-dev

Then you must run:

     make
     
or, if you don't have Automake in your system:
     g++ -o alignbucket src/alignbucket.cpp -lgmpxx -lgmp -lboost_program_options
     
Then, you have to set execution permission on the program:
     chmod u+x alignbucket


Usage
------------------
To split your fasta file into optimized buckets for 90% coverage, just run the
following command:
    ./alignbucket.sh <fasta file>

Example with Swissprot sequences:
    ./alignbucket.sh uniprot_sprot.fasta

The output will be a set of fasta file named xxx-yyy.fasta, with xxx and yyy
being the minimum and the maximum length of the sequences inside the file.


For different coverage ratios, just run
    ./alignbucket.sh <fasta file> <percentage>

Example, for 75% coverage use:
    ./alignbucket.sh uniprot_sprot.fasta 75


Authors
----------------

Giuseppe Profiti <gamma2@users.sourceforge.net>
Piero Fariselli <piero@biocomp.unibo.it>

License
----------------
The program is released under a GNU Public License version 2. Please see the
LICENSE.txt file
