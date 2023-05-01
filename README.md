# ESS-color


ESS-Color is a bioinformatics tool for constructing compressed representation of sets of k-mer sets.


## Requirements

- Linux operating system (64 bit)
- GCC >= 4.8 or a C++11 capable compiler
- Snakemake
- Git
- CMake 3.12+
- [ESSCompress](https://github.com/medvedevgroup/ESSCompress)
- [KMC](https://github.com/refresh-bio/KMC)


## Quick start

First, install all the pre-requisites and make sure the executables are in your `PATH`. Then, install additional executables from source:

    git clone https://github.com/medvedevgroup/ESSColor.git
    cd ESSColor
    bash compile.sh
    
You can move/copy ALL the executables in `ESSColor/bin` to the bin directory that is already in your PATH. For instance, considering `/usr/bin` is already in PATH, you need to run the command `mv ESSColor/bin/* /usr/bin` to move all executables for ESS-Color software. An alternative to moving/copying executables is adding the location of `ESSColor/bin` to your PATH.

After compiling, set up the `config.yaml` to test compression. (Just follow the examples/sample3 folder structure) Change to the directory `examples/sample3/k5` and use [TODO: add details of config file]

    snakemake -s pipeline_compress.smk -j 8
    
[It generates a bunch of files, but if the pipeline finished successfully, compression worked.]


## Command lines

### Color matrix generation

Generates a color matrix from a KMC database list.\\
WARNING: In the current version of the software, kmer size must be 32 at maximum.

    genmatrix [OPTION...]
    -c, --count-list arg  [Mandatory] Path to KMC database files. One line per database
    -d, --debug-verif     Debug flag to verify if the output coresponds to the input (Time consuming).
    -o, --outmatrix arg   [Mandatory] Path to the output color matrix
    -l, --kmer-list arg   [Mandatory] Path to the output ordered kmer list file
    -s, --strout          String output


Command example for a 100 ecoli matrix:

    genmatrix -c db_list.txt -o matrix.bin -l kmers.bin

The file db_list.txt must contain the paths to the 100 KMC databases (one by ecoli).
The path can be absolute or relative to the exec directory.
The software is expecting one path per line.

The file matrix.bin contains the color matrix.
The matrix has one row per kmer and 100 column (1 per sample).
The columns have the same order than the databases in the db_list.txt file.
In string format rows are separated using '\n' chars.
Each row is composed of 100 chars that are 0 or 1 depending on the presence/absence of the row kmer in the column sample.\\
In binary format, a row is a large enough multiple of 64 bits.
For our 100 samples, a row is composed of 128 bits (16 Bytes).
The xth bit of the yth byte correspond to the sample $y*8+x$.
There is no separator between successive rows.

The file kmer.bin contains the kmer list corresponding to the matrix.
In string format, there is one kmer per line.\\
In binary format, all the values inside of the file are 64 bits.
Each 64 bit is decomposed in 8 bytes little endian ordered.
First value is k, second is the number n of kmers, then are n values that are kmers.
