#include <iostream>
#include <string>
#include <algorithm>

#include "kmc_file.h"


#define DBIO_ERROR 10
#define OUTPUT_ERROR 11


class KmerMatrix
{
public:
	size_t num_datasets;
	uint64_t k;
	std::vector<uint64_t> kmers;
	std::vector<uint64_t> colors;

	KmerMatrix(std::vector<uint64_t> & dataset, uint64_t k);
	KmerMatrix(KmerMatrix && other);
	KmerMatrix& operator=(KmerMatrix&& other);

	/** Fill the to_fill vector with the uints from the row of row_idx idx.
	 * @param row_idx Index of the matrix row of interest
	 * @param to_fill A uint vector that will be modified to store the row
	 **/
	void get_row(uint64_t row_idx, std::vector<uint64_t>& to_fill);

	/** Add a sorted kmer list to the matrix. Will add 1 bit in each color row and insert absent kmers in the matrix
	 * @param kmers sorted list of kmers
	 **/
	void merge (KmerMatrix & other);

	/** Generate a color matrix in binary format. Each row is composed of a multiple of 64 bits (the minimum needed to store all the datasets from the matrix).
	 * @param outfile The filme where the matrix will be written
	 **/
	void to_color_binary_file(const std::string& outfile);
	/** Generate a color matrix in textual format. Should only be used for debug only. The process is very slow. In normal usage, please use the binary version.
	 * @param outfile The filme where the matrix will be written
	 **/
	void to_color_string_file(const std::string& outfile);

	/** Generates a binary file containing all the sorted kmers that correspond to the matrix rows.
	 * All the values are 64bits little endian uints.
	 * The 2 first values are the k value and the number of kmers in the file
	 * Then each 64bits uint is a kmer with encoding A:0, C:1, G:2, T:3
	 **/
	void to_kmer_binary_file(const std::string& outfile);
	/** Generates a text file containing all the sorted kmers that correspond to the matrix rows.
	 * The file contains one line per kmer
	 **/
	void to_kmer_string_file(const std::string& outfile);
};


/** Merges a subvector pointed by the iterator to_merge inside of the colors vector. The subvector is inserted starting at the last uint of colors at bit first_idx. The function assumes that all the bits after that position are 0.
 * @param colors Vector of colors that will be modified by the merge process
 * @param first_index The first bit of colors where to_merge will be inserted [1-64].
 * @param to_merge Iterator that contains the colors to insert
 * @param size Size of the insertion (in bits)
 **/
void merge_colors(std::vector<uint64_t> & colors, size_t first_idx, std::vector<uint64_t>::iterator to_merge, size_t size);

/** Loads a kmer list from a KMC database and sort the kmers.
 * @param db_path path to kmer database
 * @param k kmer size. This value is filled during the loading process.
 * @return Sorted list of kmers (lexicographic order)
 **/
std::vector<uint64_t> load_from_file(const std::string db_path, uint64_t& k);


/** Translate a uint64_t kmer into a string
 * @param kmer the kmer integer version
 * @param k kmer size (in nucleotides)
 * @return string version of the kmer
 **/
std::string kmer2str(uint64_t kmer, uint64_t k);


/** This class is made to amortize mergings. We do not want to merge each new dataset directly.
 * Merging process is efficient on random distributions when both list have similar sizes.
 * So, this object keep an ordered list of matricies to merge. When a matrice at a position x has a similar number of rows that x+1, a merge opperation is triggered. If needed, this operation is cascaded.
 * When all the datasets are present in this structure, a forced cascading merge has to be performed.
 **/
class CascadingMergingMatrix
{
private:
	float trigger_ratio;
	std::vector<KmerMatrix> matricies;

	void cascade_merging();
	void force_merging();

public:
	CascadingMergingMatrix(float trigger_ratio);
	void add_matrix(KmerMatrix & matrix);
	KmerMatrix& get_matrix();
};