#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "kmc_file.h"
#include "colormatrix.hpp"

using namespace std;
// #include "common.hpp"
// #include "bench_utils.hpp"
// #include "check_utils.hpp"
// // #include "lib/sshash/src/build.cpp"
// // #include "lib/sshash/src/query.cpp"
// #include "permute.cpp"

// using namespace sshash;
// class MPHFComparatoror
// {
// 	public:
// 	int k;
// 	dictionary dict;

// 	MPHFComparatoror() {
// 	}

dictionary global_dict;    
	

//     // bool operator()(uint64_t kmer1, uint64_t kmer2)
//     // {
//     //     return get_kmer_id(kmer1) < get_kmer_id(kmer2);
//     // }

int KmerMatrix::MPHFCompare(uint64_t kmer1, uint64_t kmer2)
{
	if(kmer1 > kmer2){
		return -1;
	}else if (kmer1 < kmer2){
		return +1;
	}else{
		return 0;
	}      

	std::string kmer1_str=kmer2str(kmer1, k);
	std::string kmer2_str=kmer2str(kmer2, k);
	auto answer1 = dict.lookup_advanced(kmer1_str.c_str());
	auto answer2 = dict.lookup_advanced(kmer2_str.c_str());
	assert(answer1.kmer_id != constants::invalid_uint64);
	assert(answer2.kmer_id != constants::invalid_uint64);
	if(answer1.kmer_id < answer2.kmer_id){
		return -1;
	}else if (answer1.kmer_id > answer2.kmer_id){
		return +1;
	}else{
		return 0;
	}        
}
// };


	

std::string kmer2str(uint64_t x, uint64_t k)
{
	std::stringstream ss;
	static const char nucleotides[] = {'A', 'C', 'G', 'T'};
	//char nucleotides[4] = {'A', 'C', 'T', 'G'};
	//static void uint_kmer_to_string_no_reverse(kmer_t x, char* str, uint64_t k) {
	//std::string str (k, 'A');
    for (uint64_t i = 0; i != k; ++i) {
        ss<<nucleotides[x & 3];
        x >>= 2;
    }
	std::string s = ss.str();
	return s;
}
	// static const char nucleotides[] = {'A', 'C', 'G', 'T'};
	// std::stringstream ss;

	// for (uint64_t i(0) ; i<k ; i++)
	// {
	// 	ss << nucleotides[kmer & 0b11];
	// 	kmer >>= 2;
	// }

	// std::string s = ss.str();
	// std::reverse(s.begin(), s.end());

	// return s;
//};

uint64_t KmerMatrix::mphf_get_kmer_id(uint64_t kmer1){
		std::string kmer1_str=kmer2str(kmer1, k);
		auto answer1 = dict.lookup_advanced(kmer1_str.c_str());
		assert(answer1.kmer_id != constants::invalid_uint64);
		return answer1.kmer_id;
	}	
uint64_t global_mphf_get_kmer_id(uint64_t kmer1){
		std::string kmer1_str=kmer2str(kmer1, global_dict.k());
		auto answer1 = global_dict.lookup_advanced(kmer1_str.c_str());
		assert(answer1.kmer_id != constants::invalid_uint64);
		return answer1.kmer_id;
	}	

KmerMatrix::KmerMatrix(vector<uint64_t> & dataset, uint64_t k, dictionary& dict)  
{
	num_datasets = 1; 
	this-> k = k; 
	//std::string ess_order_file
	this->kmers = dataset;
	this->colors.resize(dataset.size(), 1);
	this->dict = dict;
	//this->ess_order_file = ess_order_file;
};

KmerMatrix::KmerMatrix(KmerMatrix && other)
{
	this->num_datasets = other.num_datasets;
	this->k = other.k;
	this->kmers = std::move(other.kmers);
	this->colors = std::move(other.colors);
	this->dict = other.dict;
	//this->ess_order_file = other.ess_order_file;
};

KmerMatrix& KmerMatrix::operator=(KmerMatrix&& other)
{
	this->num_datasets = other.num_datasets;
	this->k = other.k;
	this->kmers = std::move(other.kmers);
	this->colors = std::move(other.colors);
	this->dict = other.dict;
	//this->ess_order_file = other.ess_order_file;

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
	//, std::string ess_order_file
	// auto k = 5;
    // auto m = 3;
    // dictionary dict;

    // build_configuration build_config;
    // build_config.k = k;
    // build_config.m = m;

    // build_config.canonical_parsing = true;
    // build_config.verbose = true;
    // build_config.print();

    // dict.build("/home/aur1111/s/proj4/minireal/k5/mega.essd", build_config);
    // assert(dict.k() == k);
    // std::cout<<"dict built complete";
	// //dict.streaming_query_from_file("/home/aur1111/s/proj4/minireal/k5/mega.essd");

	ofstream out(outfile);

	const uint64_t num_uints_per_row = (this->num_datasets + 63) / 64;
	vector<uint64_t> mphf_mapping(this->kmers.size()); 
	cout<<"Total number of k-mers: "<<this->kmers.size()<<endl;
	for (uint64_t kmer_idx(0) ; kmer_idx<this->kmers.size() ; kmer_idx++)
	{
		//
		//string thekmer=kmer2str(kmers[kmer_idx], this->k);
		//auto answer = dict.lookup_advanced(thekmer.c_str());
		auto answer_id = dict.lookup_uint((kmer_t) kmers[kmer_idx]);
		//assert(answer.kmer_id != constants::invalid_uint64);
		assert(answer_id != constants::invalid_uint64);
		if(answer_id <0 || answer_id>=this->kmers.size()){
			cout<<"Erroneus mphf."<<endl;
			exit(3);
		}
        mphf_mapping[answer_id] = kmer_idx;	
	}
	std::cout<<"Finished mphf mapping"<<endl;

	uint64_t kmer_idx;
	for (uint64_t iter_idx(0) ; iter_idx<this->kmers.size() ; iter_idx++)
	{
		kmer_idx = mphf_mapping[iter_idx];
		//kmer_idx = iter_idx;
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



// void KmerMatrix::to_kmer_string_file_sorted(const std::string& outfile)
// {
// 	ofstream out(outfile);

// 	uint64_t actual_idx = 0;
// 	uint64_t mapped_idx = 0;
// 	for (uint64_t kmer: this->kmers){
// 		out << kmer2str(kmer, this->k);
// 		actual_idx++;
// 	}
		

// 	out.close();
// };

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
		/// amatur remove: if (this->kmers[my_idx] < other.kmers[other_idx])
		if (MPHFCompare(this->kmers[my_idx], other.kmers[other_idx]) == -1)
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
		else if (MPHFCompare(this->kmers[my_idx], other.kmers[other_idx]) == 1)
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
void load_dictionary_sshash(dictionary& dict, std::string const& index_filename, bool verbose) {
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    if (verbose) {
        std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB)
                  << " [MB] (" << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])"
                  << std::endl;
        dict.print_info();
    }
}



// bool fancy_comparison(uint64_t a, uint64_t b) {
//   return mphf_get_kmer_id(a) <  mphf_get_kmer_id(b) ;
// }
// [] (uint64_t const& s1, uint64_t const& s2) -> bool 
// {
//        return global_mphf_get_kmer_id(a, dict) <  global_mphf_get_kmer_id(b, dict);
// };
auto sortRuleLambda = [] (uint64_t const& s1, uint64_t const& s2) -> bool
    {
       return global_mphf_get_kmer_id(s1) <  global_mphf_get_kmer_id(s2);
    };

/** Loads a kmer list from a KMC database and sort the kmers.
 * @param db_path path to kmer database
 * @param k kmer size. This value is filled during the loading process
 * @return Sorted list of kmers (lexicographic order)
 **/
vector<uint64_t> load_from_file(const string db_path, uint64_t& k, dictionary& dict)
{
	vector<uint64_t> kmers;
	global_dict = dict;

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
		//amatur comments out //kmers[idx++] = kmer_uint[0];
		//amatur adding
		// string kmerstr=kmer2str(kmer, k);
		// size_t mphf_idx = mphf_query(kmer); //TODO
		// kmers[mphf_idx] = kmer_uint[0];
		kmers[idx] =  kmer_uint[0];;
		idx++;

	}

	//amatur  comments out  
	sort(kmers.begin(), kmers.end());

	//sort(kmers.begin(), kmers.end(), sortRuleLambda);

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



