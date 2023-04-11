# ESS-color


ESS-Color is a bioinformatics tool for constructing compressed representation of sets of k-mer sets (colored dBG).

[WARNING: current main branch is not updated to work according to the README.]


## Requirements

- Linux operating system (64 bit)
- GCC >= 4.8 or a C++11 capable compiler
- Snakemake
- Git
- CMake 3.12+


## Quick start

First, install all the pre-requisites and make sure the executables are in your `PATH`. Then, install additional executables from source:

    git clone --recursive https://github.com/medvedevgroup/ESSColor.git
    cd ESSColor
    bash compile.sh
    
You can move/copy ALL the executables in `ESSColor/bin` to the bin directory that is already in your PATH. For instance, considering `/usr/bin` is already in PATH, you need to run the command `mv ESSColor/bin/* /usr/bin` to move all executables for ESS-Color software. An alternative to moving/copying executables is adding the location of `ess-color/bin` to your PATH.

    



## Usage details

### ESSColorCompress: compression of set of k-mer set
```
Syntax: ./ESSColorCompress [parameters] 

mandatory arguments:
-k [int]          k-mer size (must be >=4)
-i [input-file]   Path to input file. Input file is a single text file containing the list of multiple fasta/fastq files (one file per line)
-o [output-dir]   Path to output directory.

optional arguments:
-r [int]          Default=16. runLength: parameter runLength as described in paper
-a [int]          Default=1. Sets a threshold X, such that k-mers that appear less than X times in the input dataset are filtered out. 
-o [output-dir]   Specify output directory
-t [int]          Default=1. Number of threads.   
-v                Enable verbose mode: print more useful information.
-h                Print this Help
-V                Print version number
```
Upon successful completion, the output directory will contain 5 different files for ESS-Color representation (union.essc file, main.rrr, local.rrr, freq.txt.gz, meta.txt.gz).

#### Output files description   

* union.essc  
    * ESS-Compress representaion of union SPSS.
* meta.txt.gz   
    * Contains a gzipped text file with a header indicating value of k, value of runLength, followed by C rows each indicating name of the samples. 
* main.rrr  
   * rrr-compressed m vector as described in the paper   
* local.rrr  
   * rrr-compressed local_class_table as described in the paper
* freq.txt.gz  
   * frequency of color classes in a gzipped text file   


### ESSColorDecompress: decompression of ESSColor representation

```
Syntax: ./ESSColorDecompress [parameters] 

mandatory arguments:
-i [input-folder]    Path to input folder. Input folder contains ESS-Color representation (files: union.essc file, main.rrr, local.rrr, freq.gz, meta.txt.gz)   

optional arguments:
-f [output-type]     Default="bin". Values=["bin", str","index"]. Details are described in README.
-v                   Enable verbose mode: print more useful information.
-h                   Print this Help
-V                   Print version number
```


#### Different types of output
* When --bin flag is used (default), output files are
    * union.ess.fa.gz
    * meta.txt 
    * colmatrix.bin
* When --str flag is used, output files are 
    * union.ess.fa.gz
    * meta.txt 
    * colmatrix.txt
* When --index flag is used, output files are 
    * meta.txt 
    * labeled_colmatrix.index


#### Output files description   

* union.ess.fa.gz  
    * A FASTA file (gzipped) containing union SPSS. It contains one simplitig per entry. Header of each entry is just ">".
* meta.txt   
    * Contains: A header with value of k, followed by C rows each indicating name of the samples. 
* colmatrix.bin  
   * binary file containing non-labeled color matrix
* colmatrix.txt  
   * text file containing non-labeled color matrix
* labeled_colmatrix.index   
   * MPHF based index of SPSS which can be queried to retrieve the color vector given a list of k-mers.


## Quick start with a step-by-step example

### Compression

#### Preparing the input    
`cd example/mini_k18c4m7/`  

Let's say your 4 gzipped FASTA files are stored in folder `example/mini_k18c4m7`, named
```
sample0.fa.gz
sample1.fa.gz
sample2.fa.gz
sample3.fa.gz
```

If you wish to run`ESSColorCompress` on all 4 ".fa.gz" files, first make a list named `list_mini_k18c4m7` containing the path to all 4 files in each line.

`ls *.fa > list_mini_k18c4m7`

#### Performing the compression given the input list

Now, to compress this list using 8 threads, runLengh=16, k-mer size 18 and output to directory compress_mini_k18c4m7, run the following command:

`ESSColorCompress -i list_mini_k18c4m7 -k 18 -o compress_mini_k18c4m7/ -r 16 -t 8`

Upon successful completion, the output directory will contain 5 files for ESS-Color representation (files: union.essc file, main.rrr, local.rrr, freq.gz, meta.txt.gz).



### Decompression

#### Generating text output 
Run `ESSColorDecompress -i compress_mini_k18c4m7/ -f str`

Output `colmatrix.txt` contains the non-labeled color matrix (ESS order). Output `union.ess.fa.gz` contains gzipped simplitigs in a FASTA file.

Let's look at the first simplitig 
`$ zcat union.ess.fa.gz | head -n 2`

The first simplitig looks like this:
```
>
AAAAACAAAAAAAAAAAATTT
```

Let's look at the first 4 rows of the non-labeled color matrix
`$ head -n 4 colmatrix.txt`

The first 4 rows of the non-labeled color matrix looks like this:
```
1100
0001
0001
0001
```

The color vectors in text should be read in MSB order. So, color vector 1100 indicates that the first k-mer AAAAACAAAAAAAAAAAA is present in sample0.fa.gz and sample1.fa.gz and absent in other two. The 2nd to 4th k-mers (AAAACAAAAAAAAAAAAT, AAACAAAAAAAAAAAATT, AACAAAAAAAAAAAATTT) are present only in sample3.fa.gz.


#### Generating index and query
Run `ESSColorDecompress -i compress_mini_k18c4m7/ -f index` to build the index.

Now, perform query `ESSColorQuery -i labeled_colmatrix.index -q example/query.fa`
* It prints to stdout the color vectors of k-mers (in order) in query.fa





## Bonus sub-programs 


### ColMatDiskCompress : compress non-labeled color matrix

* Input: 
    * non-labeled color matrix (in binary/txt) in ESS order
    * Simplitig boundary vector

* Output: 
    * compressed non-labeled color matrix preserving ESS order


### genmatrix : tool for color matrix generation (copied from Yoann's branch)

Generates a color matrix from a KMC database list.

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
Each row is composed of 100 chars that are 0 or 1 depending on the presence/absence of the row kmer in the column sample.

In binary format, a row is a large enough multiple of 64 bits.
For our 100 samples, a row is composed of 128 bits (16 Bytes).
The xth bit of the yth byte correspond to the sample $y*8+x$.
There is no separator between successive rows.

The file kmer.bin contains the kmer list corresponding to the matrix.
In string format, there is one kmer per line.   

In binary format, all the values inside of the file are 64 bits.
Each 64 bit is decomposed in 8 bytes little endian ordered.
First value is k, second is the number n of kmers, then are n values that are kmers.
