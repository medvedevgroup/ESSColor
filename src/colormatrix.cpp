#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "kmc_file.h"
#include "colormatrix.hpp"

using namespace std;
using namespace sshash;


KmerMatrix::KmerMatrix(uint64_t k, uint64_t num_kmers, uint64_t num_samples, dictionary * dict)
{
	this->k = k;
	this->num_kmers = num_kmers;
	this->num_datasets = 0;
	this->numbers_per_row = (num_samples + 63) / 64;
	this->colors.resize(num_kmers * this->numbers_per_row, 0);
	this->dict = dict;
};


uint64_t KmerMatrix::mphf_get_kmer_id(uint64_t kmer)
{
	kmer = translate(kmer);
	reverse_kmer(kmer);
	kmer >>= 64 - this->k * 2;
	uint64_t position = this->dict->lookup_uint(kmer);
	if (position > this->colors.size())
	{
		cerr << "problem position " << position << endl;
		exit(1);
	}

	return position;
};


void KmerMatrix::add_dataset(vector<uint64_t> kmers)
{
	if (this->num_datasets >= 64 * this->numbers_per_row)
	{
		cerr << "Too many dataset entered in the kmer matrix" << endl;
		exit(MATRIX_OVERFLOW);
	}

	uint64_t uint_idx = this->num_datasets / 64;
	uint64_t bit_idx = this->num_datasets % 64;
	
	for (uint64_t kmer : kmers)
	{
		uint64_t row_idx = this->mphf_get_kmer_id(kmer);
		this->colors[row_idx * this->numbers_per_row + uint_idx] |= (1ul << bit_idx);
	}

	this->num_datasets += 1;
}


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
	for (uint64_t row(0) ; row<this->num_kmers ; row++)
	{
		uint64_t subvector = 0;
		for (uint64_t col(0) ; col<this->num_datasets ; col++)
		{
			if (col % 64 == 0)
			{
				subvector = this->colors[row * this->numbers_per_row + col / 64];
				// cout << "subvector " << subvector << endl;
			}

			out << ((subvector & 0b1) == 0 ? '0' : '1');
			subvector >>= 1;
		}
		out << endl;
	}
	out.close();

	// cout << "--- MATRIX ---" << endl;
	// cout << this->colors[0] << endl;
	// cout << this->colors[1] << endl;
	// cout << "XXX MATRIX XXX" << endl;
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


std::string kmer2str(uint64_t kmer, uint64_t k)
{
	static const char nucleotides[] = {'A', 'C', 'G', 'T'};
	std::stringstream ss;

	for (uint64_t i(0) ; i<k ; i++)
	{
		ss << nucleotides[kmer & 0b11];
		kmer >>= 2;
	}

	std::string s = ss.str();
	std::reverse(s.begin(), s.end());

	return s;
};


const uint64_t encode[256] = {0, 1, 3, 2, 4, 5, 7, 6, 12, 13, 15, 14, 8, 9, 11, 10, 16, 17, 19, 18, 20, 21, 23, 22, 28, 29, 31, 30, 24, 25, 27, 26, 48, 49, 51, 50, 52, 53, 55, 54, 60, 61, 63, 62, 56, 57, 59, 58, 32, 33, 35, 34, 36, 37, 39, 38, 44, 45, 47, 46, 40, 41, 43, 42, 64, 65, 67, 66, 68, 69, 71, 70, 76, 77, 79, 78, 72, 73, 75, 74, 80, 81, 83, 82, 84, 85, 87, 86, 92, 93, 95, 94, 88, 89, 91, 90, 112, 113, 115, 114, 116, 117, 119, 118, 124, 125, 127, 126, 120, 121, 123, 122, 96, 97, 99, 98, 100, 101, 103, 102, 108, 109, 111, 110, 104, 105, 107, 106, 192, 193, 195, 194, 196, 197, 199, 198, 204, 205, 207, 206, 200, 201, 203, 202, 208, 209, 211, 210, 212, 213, 215, 214, 220, 221, 223, 222, 216, 217, 219, 218, 240, 241, 243, 242, 244, 245, 247, 246, 252, 253, 255, 254, 248, 249, 251, 250, 224, 225, 227, 226, 228, 229, 231, 230, 236, 237, 239, 238, 232, 233, 235, 234, 128, 129, 131, 130, 132, 133, 135, 134, 140, 141, 143, 142, 136, 137, 139, 138, 144, 145, 147, 146, 148, 149, 151, 150, 156, 157, 159, 158, 152, 153, 155, 154, 176, 177, 179, 178, 180, 181, 183, 182, 188, 189, 191, 190, 184, 185, 187, 186, 160, 161, 163, 162, 164, 165, 167, 166, 172, 173, 175, 174, 168, 169, 171, 170};

/** Translate a kmer from the alphabetic order encoding {A, C, G, T} to ascii order {A, C, T, G}
 **/
uint64_t translate(uint64_t kmer)
{
	uint64_t translation = 0;

	for (uint64_t bx(0) ; bx<8 ; bx++)
	{
		uint64_t byte = kmer & 0xFF;
		kmer >>= 8;
		translation |= encode[byte] << (8 * bx);
	}

	return translation;
};



/** Reverse a kmer without complementing it. This function do not pas the bits on the right.
 * @param kmer Kmer to revert
 **/
void reverse_kmer(uint64_t & kmer)
{
	uint64_t shifted;
	// 1/2
	shifted  = kmer << 32;
	shifted &= 0xFFFFFFFF00000000;
	kmer   >>= 32;
	kmer    &= 0x00000000FFFFFFFF;
	kmer    |= shifted;

	// 1/4
	shifted  = kmer << 16;
	shifted &= 0xFFFF0000FFFF0000;
	kmer   >>= 16;
	kmer    &= 0x0000FFFF0000FFFF;
	kmer    |= shifted;

	// 1/8
	shifted  = kmer << 8;
	shifted &= 0xFF00FF00FF00FF00;
	kmer   >>= 8;
	kmer    &= 0x00FF00FF00FF00FF;
	kmer    |= shifted;

	// 1/16
	shifted  = kmer << 4;
	shifted &= 0xF0F0F0F0F0F0F0F0;
	kmer   >>= 4;
	kmer    &= 0x0F0F0F0F0F0F0F0F;
	kmer    |= shifted;

	// 1/32
	shifted  = kmer << 2;
	shifted &= 0xCCCCCCCCCCCCCCCC;
	kmer   >>= 2;
	kmer    &= 0x3333333333333333;
	kmer    |= shifted;
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
		kmers[idx] =  kmer_uint[0];;
		idx++;
	}

	return kmers;
};
