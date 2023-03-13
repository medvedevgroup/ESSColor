#include <iostream>
#include <string>
#include <algorithm>

#include "kmc_file.h"

using namespace std;


#define ERROR_DBIO 1


class KmerMatrix
{
public:
	size_t num_datasets;
	vector<uint64_t> kmers;
	vector<vector<uint64_t> > colors;

	KmerMatrix(vector<uint64_t> & dataset);
	KmerMatrix(KmerMatrix && other);
	KmerMatrix& operator=(KmerMatrix&& other);

	/** Add a sorted kmer list to the matrix. Will add 1 bit in each color row and insert absent kmers in the matrix
	 * @param kmers sorted list of kmers
	 **/
	void merge (const KmerMatrix & other);
};


/** Merges the to_merge vector inside of the colors vector. The to_merge vector is inserted from the position first_idx in the colors vector assuming that all the bits after that position are 0.
 * @param colors Vector of colors that will be modified by the merge process
 * @param to_merge Vector of colors to insert
 * @param first_index The first bit of colors where to_merge will be inserted.
 * @param size Size of the insertion (in bits)
 **/
void merge_colors(vector<uint64_t> & colors, const vector<uint64_t> to_merge, size_t first_idx, size_t size);

/** Loads a kmer list from a KMC database and sort the kmers.
 * @param db_path path to kmer database
 * @return Sorted list of kmers (lexicographic order)
 **/
vector<uint64_t> load_from_file(string db_path);


/** This class is made to amortize mergings. We do not want to merge each new dataset directly.
 * Merging process is efficient on random distributions when both list have similar sizes.
 * So, this object keep an ordered list of matricies to merge. When a matrice at a position x has a similar number of rows that x+1, a merge opperation is triggered. If needed, this operation is cascaded.
 * When all the datasets are present in this structure, a forced cascading merge has to be performed.
 **/
class CascadingMergingMatrix
{
private:
	float trigger_ratio;
	vector<KmerMatrix> matricies;

	void cascade_merging();
	void force_merging();

public:
	CascadingMergingMatrix(float trigger_ratio);
	void add_matrix(KmerMatrix & matrix);
	KmerMatrix get_matrix();
};