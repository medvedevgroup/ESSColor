#include <iostream>
#include <string>
#include <algorithm>

#include "kmc_file.h"
#include "colormatrix.hpp"

using namespace std;


KmerMatrix::KmerMatrix(vector<uint64_t> & dataset) : num_datasets(1)
{
	this->kmers = dataset;
	this->colors.resize(dataset.size(), vector<uint64_t>(1, 1));
};

KmerMatrix::KmerMatrix(KmerMatrix && other)
{
	this->num_datasets = other.num_datasets;
	this->kmers = std::move(other.kmers);
	this->colors = std::move(other.colors);
};

KmerMatrix& KmerMatrix::operator=(KmerMatrix&& other)
{

	this->num_datasets = other.num_datasets;
	this->kmers = std::move(other.kmers);
	this->colors = std::move(other.colors);

	return *this;
};

/** Add a sorted kmer list to the matrix. Will add 1 bit in each color row and insert absent kmers in the matrix
 * @param kmers sorted list of kmers
 **/
void KmerMatrix::merge (const KmerMatrix & other)
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


/** Merges the to_merge vector inside of the colors vector. The to_merge vector is inserted from the position first_idx in the colors vector assuming that all the bits after that position are 0.
 * @param colors Vector of colors that will be modified by the merge process
 * @param to_merge Vector of colors to insert
 * @param first_index The first bit of colors where to_merge will be inserted.
 * @param size Size of the insertion (in bits)
 **/
void merge_colors(vector<uint64_t> & colors, const vector<uint64_t> to_merge, size_t first_idx, size_t size)
{
	assert(first_idx != 0);

	size_t colors_idx = (first_idx - 1) / 64;
	size_t to_merge_idx = 0;

	// occupied bits from the last uint64_t of colors
	const uint64_t occupied_bits = ((first_idx - 1) % 64) + 1;
	const uint64_t free_bits = 64 - occupied_bits;
	const uint64_t first_half_mask = occupied_bits == 64 ? 0 : (1 << free_bits) - 1;
	const uint64_t second_half_mask = ~first_half_mask;

	// Merge the vector one uint64_t at a time
	while (size > 0)
	{
		// Extract the first part of the current merged uint
		uint64_t first_half = to_merge[to_merge_idx] & first_half_mask;
		// Insert the first half in the current color uint
		colors[colors_idx] |= first_half << occupied_bits;

		size -= min(size, free_bits);
		if (size == 0)
			break;
		colors_idx += 1;

		// Extract the last part of the current to_merge uint
		uint64_t second_half = to_merge[to_merge_idx] & second_half_mask;
		second_half >>= free_bits;
		colors.push_back(second_half);

		size -= min(size, occupied_bits);
		to_merge_idx += 1;
	}
};

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
};


CascadingMergingMatrix::CascadingMergingMatrix(float trigger_ratio) : trigger_ratio(trigger_ratio) {};

void CascadingMergingMatrix::add_matrix(KmerMatrix & matrix)
{
	this->matricies.push_back(move(matrix));
	this->cascade_merging();
};

void CascadingMergingMatrix::cascade_merging()
{
	while (true)
	{
		// No merging for 0 or 1 matrix
		if (this->matricies.size() < 2)
			return;

		// Only cascade on proper trigger ratio
		size_t n = this->matricies.size();
		if (this->matricies[n-1].kmers.size() < this->trigger_ratio * this->matricies[n-2].kmers.size())
			return;

		this->matricies[n-2].merge(this->matricies[n-1]);
		this->matricies.pop_back();
	}
};

void CascadingMergingMatrix::force_merging()
{
	for (size_t n(this->matricies.size()) ; n > 1 ; n--)
	{
		this->matricies[n-2].merge(this->matricies[n-1]);
		this->matricies.pop_back();
	}
};


