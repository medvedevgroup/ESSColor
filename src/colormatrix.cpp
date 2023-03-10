#include <iostream>
#include <string>
#include <algorithm>

#include "kmc_file.h"

using namespace std;


#define ERROR_DBIO 1



/** Loads a kmer list from a KMC database and sort the kmers.
 * @param db_path path to kmer database
 * @return Sorted list of kmers (lexicographic order)
 **/
vector<uint64_t> load_from_file(string db_path)
{
	vector<uint64_t> kmers;

	// Opening the database
	CKMCFile db;
	if (!db.OpenForListing(db_path))
	{
		cerr << "Impossible to open DB " << db_path << endl;
		exit(ERROR_DBIO);
	}

	// Get the metadata from the database
	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
	uint64 _total_kmers, _max_count;

	db.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
	assert(_kmer_length <= 32);
	kmers.resize(_total_kmers);
	cout << _kmer_length << " " << _total_kmers << endl;

	// Load the kmers
	size_t idx = 0;
	CKmerAPI kmer(_kmer_length);
	vector<uint64> kmer_uint;
	uint32 counter;
	while (db.ReadNextKmer(kmer, counter)) {
		kmer.to_long(kmer_uint);
		kmers[idx++] = kmer_uint[0];
	}

	sort(kmers.begin(), kmers.end());

	return kmers;
}



class KmerMatrix
{
public:
	size_t num_datasets;
	vector<uint64_t> kmers;
	vector<vector<uint64_t> > colors;

	KmerMatrix() : num_datasets(0) {};

	void add_dataset (vector<uint64_t> kmers)
	{
		// Verify colors overflow
		if (num_datasets % 64 == 0)
		{
			// Add uints for new colors
			for (vector<uint64_t> & vector : colors)
			{
				vector.push_back(0);
			}
		}

		// Add the kmers from the dataset using a sorted list fusion procedure
		size_t row_idx=0, dataset_idx=0;
	};
};


int main(int argc, char const *argv[])
{
	string db_name("data/sarscov2_k31");
	vector<uint64_t> kmers = load_from_file(db_name);

	for (uint64_t kmer : kmers)
		cout << kmer << endl;

	return 0;
}
