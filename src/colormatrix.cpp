#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "kmc_file.h"
#include "colormatrix.hpp"

using namespace std;


KmerMatrix::KmerMatrix(vector<uint64_t> & dataset, uint64_t k) : num_datasets(1), k(k)
{
	this->kmers = dataset;
	this->colors.resize(dataset.size(), 1);
};

KmerMatrix::KmerMatrix(KmerMatrix && other)
{
	this->num_datasets = other.num_datasets;
	this->k = other.k;
	this->kmers = std::move(other.kmers);
	this->colors = std::move(other.colors);
};

KmerMatrix& KmerMatrix::operator=(KmerMatrix&& other)
{
	this->num_datasets = other.num_datasets;
	this->k = other.k;
	this->kmers = std::move(other.kmers);
	this->colors = std::move(other.colors);

	return *this;
};

void KmerMatrix::get_row(uint64_t row_idx, vector<uint64_t>& to_fill)
{
	size_t nb_uints = (this->num_datasets + 63) / 64;
	to_fill.resize(nb_uints);

	size_t row_offset = nb_uints * row_idx;

	for (size_t i(0) ; i<nb_uints ; i++)
	{
		to_fill[i] = this->colors[row_offset + i];
	}
}

void KmerMatrix::to_color_string_file(const std::string& outfile)
{
	ofstream out(outfile);

	const uint64_t num_uints_per_row = (this->num_datasets + 63) / 64;
	for (uint64_t kmer_idx(0) ; kmer_idx<this->kmers.size() ; kmer_idx++)
	{
		for (uint64_t dataset_idx(0) ; dataset_idx<this->num_datasets ; dataset_idx++)
		{
			uint64_t subvector = this->colors[kmer_idx * num_uints_per_row + dataset_idx / 64];
			uint64_t color = (subvector >> (dataset_idx % 64)) & 0b1;
			out << (color == 0 ? '0' : '1');
		}
		out << endl;
	}

	out.close();
}


void KmerMatrix::to_color_binary_file(const std::string& outfile)
{
	ofstream out(outfile, ios::out | ios::binary);

	uint8_t array[8];
	for (const uint64_t color_subvector : this->colors)
	{
		array[0] =  color_subvector        & 0xFF;
		array[1] = (color_subvector >>  8) & 0xFF;
		array[2] = (color_subvector >> 16) & 0xFF;
		array[3] = (color_subvector >> 24) & 0xFF;
		array[4] = (color_subvector >> 32) & 0xFF;
		array[5] = (color_subvector >> 40) & 0xFF;
		array[6] = (color_subvector >> 48) & 0xFF;
		array[7] = (color_subvector >> 56) & 0xFF;
		out.write((char *)array, 8);
	}

	out.close();
}

/** Generates a binary file containing all the sorted kmers that correspond to the matrix rows.
 * All the values are 64bits little endian uints.
 * The 2 first values are the k value and the number of kmers in the file
 * Then each 64bits uint is a kmer with encoding A:0, C:1, G:2, T:3
 **/
void KmerMatrix::to_kmer_binary_file(const std::string& outfile)
{
	ofstream out(outfile, ios::out | ios::binary);

	uint8_t array[8];

	// k value
	array[0] =  this->k        & 0xFF;
	array[1] = (this->k >>  8) & 0xFF;
	array[2] = (this->k >> 16) & 0xFF;
	array[3] = (this->k >> 24) & 0xFF;
	array[4] = (this->k >> 32) & 0xFF;
	array[5] = (this->k >> 40) & 0xFF;
	array[6] = (this->k >> 48) & 0xFF;
	array[7] = (this->k >> 56) & 0xFF;
	out.write((char *)array, 8);

	// number of kmers
	array[0] =  this->kmers.size()        & 0xFF;
	array[1] = (this->kmers.size() >>  8) & 0xFF;
	array[2] = (this->kmers.size() >> 16) & 0xFF;
	array[3] = (this->kmers.size() >> 24) & 0xFF;
	array[4] = (this->kmers.size() >> 32) & 0xFF;
	array[5] = (this->kmers.size() >> 40) & 0xFF;
	array[6] = (this->kmers.size() >> 48) & 0xFF;
	array[7] = (this->kmers.size() >> 56) & 0xFF;
	out.write((char *)array, 8);

	for (const uint64_t kmer : this->kmers)
	{
		array[0] =  kmer        & 0xFF;
		array[1] = (kmer >>  8) & 0xFF;
		array[2] = (kmer >> 16) & 0xFF;
		array[3] = (kmer >> 24) & 0xFF;
		array[4] = (kmer >> 32) & 0xFF;
		array[5] = (kmer >> 40) & 0xFF;
		array[6] = (kmer >> 48) & 0xFF;
		array[7] = (kmer >> 56) & 0xFF;
		out.write((char *)array, 8);
	}

	out.close();
}

/** Generates a text file containing all the sorted kmers that correspond to the matrix rows.
 * The file contains one line per kmer
 **/
void KmerMatrix::to_kmer_string_file(const std::string& outfile)
{
	ofstream out(outfile);

	for (uint64_t kmer: this->kmers)
		out << kmer2str(kmer, this->k);

	out.close();
};



/** Add a sorted kmer list to the matrix. Will add 1 bit in each color row and insert absent kmers in the matrix. Both matrices are destroyed in the process
 * @param kmers sorted list of kmers
 **/
void KmerMatrix::merge (KmerMatrix & other)
{
	// cout << "MERGING" << endl;
	// Reserve memory for the merge
	// size_t max_expected_kmers = this->kmers.size() + other.kmers.size();
	size_t cleaning_threshold = max(this->kmers.size(), other.kmers.size()) / 10;
	// const size_t num_colors = this->num_datasets + other.num_datasets;
	// const size_t num_uints_per_row = (num_colors + 63) / 64;
	size_t my_color_uint_size = (this->num_datasets + 63) / 64;
	size_t other_color_uint_size = (other.num_datasets + 63) / 64;
	size_t my_color_offset = ((this->num_datasets - 1) % 64) + 1;
	// Empty vectors for the merges
	vector<uint64_t> my_empty_color_vector(my_color_uint_size, 0);
	vector<uint64_t> other_empty_color_vector(other_color_uint_size, 0);

	// New datastructs for the merge
	vector<uint64_t> new_kmers;
	vector<uint64_t> new_colors;
	// new_kmers.reserve(max_expected_kmers);
	// new_colors.reserve(max_expected_kmers * num_uints_per_row);
	size_t new_idx = 0;

	// Setup variable for the merges
	size_t my_idx = 0, other_idx = 0;
	vector<uint64_t>::iterator my_color_iter = this->colors.begin();
	vector<uint64_t>::iterator other_color_iter = other.colors.begin();

	// Add the kmers from the dataset using a sorted list fusion procedure
	while(my_idx < this->kmers.size() and other_idx < other.kmers.size())
	{
		if (this->kmers[my_idx] < other.kmers[other_idx])
		{
			new_kmers.push_back(this->kmers[my_idx]);
			// cout << "< " << new_colors.size() << " -> ";
			new_colors.insert(new_colors.end(), my_color_iter, my_color_iter + my_color_uint_size);
			// cout << new_colors.size() << " -> ";
			merge_colors(
				new_colors, my_color_offset,
				other_empty_color_vector.begin(), other.num_datasets
			);
			// cout << new_colors.size() << endl;
			my_idx += 1;
			my_color_iter += my_color_uint_size;
		}
		else if (this->kmers[my_idx] > other.kmers[other_idx])
		{
			new_kmers.push_back(other.kmers[other_idx]);
			// cout << "> " << new_colors.size() << " -> ";
			new_colors.insert(new_colors.end(), my_empty_color_vector.begin(), my_empty_color_vector.end());
			// cout << new_colors.size() << " -> ";
			merge_colors(
				new_colors, my_color_offset,
				other_color_iter, other.num_datasets
			);
			// cout << new_colors.size() << endl;
			other_idx += 1;
			other_color_iter += other_color_uint_size;
		}
		else
		{
			new_kmers.push_back(this->kmers[my_idx]);
			// cout << "= " << new_colors.size() << " -> ";
			new_colors.insert(new_colors.end(), my_color_iter, my_color_iter + my_color_uint_size);
			// cout << new_colors.size() << " -> ";
			merge_colors(
				new_colors, my_color_offset,
				other_color_iter, other.num_datasets
			);
			// cout << new_colors.size() << endl;
			my_idx += 1; my_color_iter += my_color_uint_size;
			other_idx += 1; other_color_iter += other_color_uint_size;
		}
		new_idx += 1;
		// cout << new_kmers.size() << " " << new_colors.size() << endl;

		// Memory saving procedure
		if ((my_idx + other_idx) >= cleaning_threshold)
		{
			// cout << "cleaning " << my_idx << "+" << other_idx << " <=> " << cleaning_threshold << endl;
			// cout << this->kmers.size() << " " << other.kmers.size() << " " << new_kmers.size() << endl;
			// Modify current matrix
			this->kmers.erase(this->kmers.begin(), this->kmers.begin()+my_idx);
			my_idx = 0;
			this->colors.erase(this->colors.begin(), my_color_iter);
			my_color_iter = this->colors.begin();

			// Modify the other matrix
			other.kmers.erase(other.kmers.begin(), other.kmers.begin()+other_idx);
			other_idx = 0;
			other.colors.erase(other.colors.begin(), other_color_iter);
			other_color_iter = other.colors.begin();
		}
	}

	// Add the last values
	for (size_t idx=my_idx ; idx<this->kmers.size() ; idx++)
	{
		new_kmers.push_back(this->kmers[idx]);
		new_colors.insert(new_colors.end(), my_color_iter, my_color_iter + my_color_uint_size);
		merge_colors(
			new_colors, my_color_offset,
			other_empty_color_vector.begin(), other.num_datasets
		);
		my_color_iter += my_color_uint_size;
	}
	for (size_t idx=other_idx ; idx<other.kmers.size() ; idx++)
	{
		new_kmers.push_back(other.kmers[idx]);
		new_colors.insert(new_colors.end(), my_empty_color_vector.begin(), my_empty_color_vector.end());
		merge_colors(
			new_colors, my_color_offset,
			other_color_iter, other.num_datasets
		);
		other_color_iter += other_color_uint_size;
	}

	// Replace previous vectors
	this->kmers = new_kmers;
	this->colors = new_colors;
	this->num_datasets += other.num_datasets;
	// cout << "/merging" << endl;
};


/** Merges a subvector pointed by the iterator to_merge inside of the colors vector. The subvector is inserted starting at the last uint of colors at bit first_idx. The function assumes that all the bits after that position are 0.
 * @param colors Vector of colors that will be modified by the merge process
 * @param first_index The first bit of colors where to_merge will be inserted [1-64].
 * @param to_merge Iterator that contains the colors to insert
 * @param size Size of the insertion (in bits)
 **/
void merge_colors(vector<uint64_t> & colors, size_t first_idx, vector<uint64_t>::iterator to_merge, size_t size)
{
	assert(first_idx != 0 and first_idx <= 64);

	size_t colors_idx = colors.size() - 1;

	// occupied bits from the last uint64_t of colors
	const uint64_t occupied_bits = first_idx;
	const uint64_t free_bits = 64 - occupied_bits;
	const uint64_t first_half_mask = occupied_bits == 64 ? 0 : (1 << free_bits) - 1;
	const uint64_t second_half_mask = ~first_half_mask;

	// Merge the vector one uint64_t at a time
	while (size > 0)
	{
		uint64_t to_merge_uint = *to_merge;
		// Extract the first part of the current merged uint
		uint64_t first_half = to_merge_uint & first_half_mask;
		// Insert the first half in the current color uint
		colors[colors_idx] |= first_half << occupied_bits;

		size -= min(size, free_bits);
		if (size == 0)
			break;
		colors_idx += 1;

		// Extract the last part of the current to_merge uint
		uint64_t second_half = to_merge_uint & second_half_mask;
		second_half >>= free_bits;
		colors.push_back(second_half);

		size -= min(size, occupied_bits);
		to_merge++;
	}
};

/** Loads a kmer list from a KMC database and sort the kmers.
 * @param db_path path to kmer database
 * @param k kmer size. This value is filled during the loading process
 * @return Sorted list of kmers (lexicographic order)
 **/
vector<uint64_t> load_from_file(const string db_path, uint64_t& k)
{
	vector<uint64_t> kmers;

	// Opening the database
	CKMCFile db;
	if (!db.OpenForListing(db_path))
	{
		cerr << "Impossible to open DB " << db_path << endl;
		exit(DBIO_ERROR);
	}

	// Get the metadata from the database
	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
	uint64 _total_kmers, _max_count;

	db.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
	assert(_kmer_length <= 32);
	kmers.resize(_total_kmers);
	k = _kmer_length;
	
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

KmerMatrix& CascadingMergingMatrix::get_matrix()
{
	this->force_merging();
	return this->matricies[0];
};



string kmer2str(uint64_t kmer, uint64_t k)
{
	static const char nucleotides[] = {'A', 'C', 'G', 'T'};
	stringstream ss;

	for (uint64_t i(0) ; i<k ; i++)
	{
		ss << nucleotides[kmer & 0b11];
		kmer >>= 2;
	}

	string s = ss.str();
	std::reverse(s.begin(), s.end());

	return s;
};
