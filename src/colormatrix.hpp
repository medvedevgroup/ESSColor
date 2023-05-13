#pragma once 

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "kmc_file.h"
#include "dictionary.hpp"


#define DBIO_ERROR 10
#define OUTPUT_ERROR 11
#define MATRIX_OVERFLOW 12




class KmerMatrix
{
public:
	size_t num_datasets;
	uint64_t num_kmers;
	// The number of uint64 needed to store a row of the matrix
	uint64_t numbers_per_row;
	uint64_t k;
	std::vector<uint64_t> colors;
	sshash::dictionary * dict;
	
	/** Creates a kmer matrix with the right amount of memory to store num_kmers * num_samples bits
	 * @param k kmer size
	 * @param num_kmers Number of kmers to store in the matrix
	 * @param num_samples Number of samples (number of columns)
	 * @param dict sshash dictrionary that store the index of each kmer.
	 **/
	KmerMatrix(uint64_t k, uint64_t num_kmers, uint64_t num_samples, sshash::dictionary * dict);

	/** Set a 1 to each position of the matrix that corresponds to a kmer from the sample.
	 * @param kmers Kmers from the sample
	 **/
	void add_dataset(std::vector<uint64_t> kmers);

	uint64_t mphf_get_kmer_id(uint64_t kmer);

	/** Fill the to_fill vector with the uints from the row of row_idx idx.
	 * @param row_idx Index of the matrix row of interest
	 * @param to_fill A uint vector that will be modified to store the row
	 **/
	void get_row(uint64_t row_idx, std::vector<uint64_t>& to_fill);

	/** Generate a color matrix in binary format. Each row is composed of a multiple of 64 bits (the minimum needed to store all the datasets from the matrix).
	 * @param outfile The filme where the matrix will be written
	 **/
	void to_color_binary_file(const std::string& outfile);
	/** Generate a color matrix in textual format. Should only be used for debug only. The process is very slow. In normal usage, please use the binary version.
	 * @param outfile The filme where the matrix will be written
	 **/
	void to_color_string_file(const std::string& outfile);
};

/** Loads a kmer list from a KMC database and sort the kmers.
 * @param db_path path to kmer database
 * @param k kmer size. This value is filled during the loading process.
 * @return Sorted list of kmers (lexicographic order)
 **/
std::vector<uint64_t> load_from_file(const std::string db_path, uint64_t& k);


/** Translates a uint64_t kmer into a string
 * @param kmer the kmer integer version
 * @param k kmer size (in nucleotides)
 * @return string version of the kmer
 **/
std::string kmer2str(uint64_t kmer, uint64_t k);

/** Translate a kmer from the alphabetic order encoding {A, C, G, T} to ascii order {A, C, T, G}
 **/
uint64_t translate(uint64_t kmer);

/** Reverse a kmer without complementing it
 * @param kmer Kmer to revert
 **/
void reverse_kmer(uint64_t & kmer);
