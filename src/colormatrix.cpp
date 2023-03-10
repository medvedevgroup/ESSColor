#include <iostream>
#include <string>
#include <algorithm>

#include "kmc_file.h"

using namespace std;


#define ERROR_DBIO 1

vector<uint64_t> load_from_file(string db_path);
void merge_colors(vector<uint64_t> & colors, const vector<uint64_t> to_merge, size_t first_idx, size_t size);



class KmerMatrix
{
public:
	size_t num_datasets;
	vector<uint64_t> kmers;
	vector<vector<uint64_t> > colors;

	KmerMatrix(vector<uint64_t> & dataset) : num_datasets(1)
	{
		this->kmers = dataset;
		this->colors.resize(dataset.size(), vector<uint64_t>(1, 1));
	};

	/** Add a sorted kmer list to the matrix. Will add 1 bit in each color row and insert absent kmers in the matrix
	 * @param kmers sorted list of kmers
	 **/
	void merge (const KmerMatrix & other)
	{
		// Reserve memory for the merge
		size_t max_expected_kmers = this->kmers.size() + other.kmers.size();
		vector<uint64_t> new_kmers;
		vector<vector<uint64_t> > new_colors;
		new_kmers.reserve(max_expected_kmers);
		new_colors.reserve(max_expected_kmers);
		size_t new_idx = 0;

		// Add the kmers from the dataset using a sorted list fusion procedure
		size_t my_idx = 0, other_idx = 0;
		size_t my_color_size = this->colors[0].size();
		size_t other_color_size = other.colors[0].size();
		while(my_idx < this->kmers.size() and other_idx < other.kmers.size())
		{
			if (this->kmers[my_idx] < other.kmers[other_idx])
			{
				new_kmers.push_back(this->kmers[my_idx]);
				new_colors.push_back(this->colors[my_idx]);
				merge_colors(
					new_colors[new_colors.size()-1],
					vector<uint64_t>(other_color_size, 0),
					this->num_datasets, other.num_datasets
				);
				my_idx += 1;
			}
			else if (this->kmers[my_idx] > other.kmers[other_idx])
			{
				new_kmers.push_back(other.kmers[other_idx]);
				new_colors.push_back(vector<uint64_t>(my_color_size, 0));
				merge_colors(
					new_colors[new_colors.size()-1],
					other.colors[other_idx],
					this->num_datasets, other.num_datasets
				);
				other_idx += 1;
			}
			else
			{
				new_kmers.push_back(this->kmers[my_idx]);
				new_colors.push_back(this->colors[my_idx]);
				merge_colors(
					new_colors[new_colors.size()-1],
					other.colors[other_idx],
					this->num_datasets, other.num_datasets
				);
				my_idx += 1;
				other_idx += 1;
			}

			new_idx += 1;
		}

		// Add the last values
		for (size_t idx=my_idx ; idx<this->kmers.size() ; idx++)
		{
			new_kmers.push_back(this->kmers[idx]);
			new_colors.push_back(this->colors[idx]);
			merge_colors(
				new_colors[new_colors.size()-1],
				vector<uint64_t>(other_color_size, 0),
				this->num_datasets, other.num_datasets
			);
		}
		for (size_t idx=other_idx ; idx<other.kmers.size() ; idx++)
		{
			new_kmers.push_back(other.kmers[idx]);
			new_colors.push_back(vector<uint64_t>(my_color_size, 0));
			merge_colors(
				new_colors[new_colors.size()-1],
				other.colors[idx],
				this->num_datasets, other.num_datasets
			);
			other_idx += 1;
		}

		// Replace previous vectors
		this->kmers = new_kmers;
		this->colors = new_colors;
		this->num_datasets += other.num_datasets;
	};
};


/** Merges the to_merge vector inside of the colors vector. The to_merge vector is inserted from the position first_idx in the colors vector assuming that all the bits after that position are 0.
 * @param colors Vector of colors that will be modified by the merge process
 * @param to_merge Vector of colors to insert
 * @param first_index The first bit of colors where to_merge will be inserted.
 * @param size Size of the insertion (in bits)
 **/
void merge_colors(vector<uint64_t> & colors, const vector<uint64_t> to_merge, size_t first_idx, size_t size)
{
	while (size > 0)
	{
		// TODO
	}
}

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


int main(int argc, char const *argv[])
{
	string db_name("data/sarscov2_k31");
	vector<uint64_t> kmers = load_from_file(db_name);

	for (uint64_t kmer : kmers)
		cout << kmer << endl;

	return 0;
}
