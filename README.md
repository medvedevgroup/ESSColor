# ESS-color


ESS-Color is a bioinformatics tool for constructing compressed representation of sets of k-mer sets.


## Requirements

- Linux operating system (64 bit)
- GCC >= 4.8 or a C++11 capable compiler
- Snakemake
- Git
- CMake 3.1+
- [ESSCompress](https://github.com/medvedevgroup/ESSCompress)
- [joinCounts](https://github.com/Transipedia/dekupl-joinCounts)
- [KMC](https://github.com/refresh-bio/KMC)
- sdsl-lite
- cmph


## Quick start

First, install all the pre-requisites and make sure the executables are in your `PATH`. Then, install additional executables from source:

    git clone https://github.com/medvedevgroup/ESSColor
    cd ESSColor
    bash compile.sh
    
You can move/copy ALL the executables in `ess-color/bin` to the bin directory that is already in your PATH. For instance, considering `/usr/bin` is already in PATH, you need to run the command `mv ess-color/bin/* /usr/bin` to move all executables for ESS-Color software. An alternative to moving/copying executables is adding the location of `ess-color/bin` to your PATH.


After compiling, set up the `config.yaml`. Change to the directory and use [TODO: add details of config file]

    snakemake -s Snakefile.smk -j 8
    
    

## Example to test compression and decompression 
    
`cd example/mini_k18c4m7/`  

It takes a list of ordered color vectors (`col_bitmatrix`), list of unique color vectors (`uniq_ms.txt`) [TODO: we can just pass list, and unique can be found in code], number of colors (-c), number of unique colors (-m), number of k-mers (-k) or size of col_bitmatrix, ess/simplitig boundary file (`ess_boundary_bit.txt`) where 1 indicates start of simplitig.

[TODO: instead of `ess_boundary_bit.txt`, just essCompressed file should be provided]

Parameter (-x): from paper: value `runDivisor`


To compress:   
     
`ess_color_compress -d col_bitmatrix -i uniq_ms.txt -c 4 -m 7 -k 18 -s ess_boundary_bit.txt -x 16`   


Now, to test that decompression works, run decompress:  

[TODO: Ideally, decompress should only take the compressed bb_main, bb_local_table etc., along with a file where values for -c, -m, -k, -x options are provided. But for development phase (i.e. helpful for debugging), we pass all the information from command line]

`ess_color_decompress -d col_bitmatrix -i uniq_ms.txt -c 4 -m 7 -k 18 -s ess_boundary_bit.txt -x 16` 

Now, the decompressed color matrix is in `dec_ess_color`.

Then, check that `col_bitmatrix` and `dec_ess_color` are the same with unix `diff` command. This confirms decompression was correct.





