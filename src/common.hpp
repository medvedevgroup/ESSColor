#pragma once

#ifndef COMMON_MYTYPES_H
#define COMMON_MYTYPES_H
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/dictionary.hpp"

#include <vector>

// #include "bench_utils.hpp"
// #include "check_utils.hpp"
// // #include "lib/sshash/src/build.cpp"
// // #include "lib/sshash/src/query.cpp"
// #include "permute.cpp"

using namespace sshash;

namespace sshash {

// void random_kmer(char* kmer, uint64_t k) {
//     for (uint64_t i = 0; i != k; ++i) kmer[i] = "ACGT"[rand() % 4];
// }

void load_dictionary(dictionary& dict, std::string const& index_filename, bool verbose) {
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    if (verbose) {
        std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB)
                  << " [MB] (" << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])"
                  << std::endl;
        dict.print_info();
    }
}

}  // namespace sshash


using namespace sshash;

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
class MPHFComparatoror
{
	public:
	int k;
	dictionary dict;

	MPHFComparatoror() {
	}

    MPHFComparatoror(int k) {
		this->k = k;
		k = 5;
		auto m = 3;
		//dictionary dict;

		build_configuration build_config;
		build_config.k = k;
		build_config.m = m;

		build_config.canonical_parsing = true;
		build_config.verbose = true;
		build_config.print();

		dict.build("/home/aur1111/s/proj4/minireal/k5/mega.essd", build_config);
		assert(dict.k() == k);
		std::cout<<"dict built complete";
		//dict.streaming_query_from_file("/home/aur1111/s/proj4/minireal/k5/mega.essd");

	}

    bool operator()(uint64_t kmer1, uint64_t kmer2)
    {
		std::string kmer1_str=kmer2str(kmer1, k);
		std::string kmer2_str=kmer2str(kmer2, k);
		auto answer1 = dict.lookup_advanced(kmer1_str.c_str());
		auto answer2 = dict.lookup_advanced(kmer2_str.c_str());
		assert(answer1.kmer_id != constants::invalid_uint64);
		assert(answer2.kmer_id != constants::invalid_uint64);
        return answer1.kmer_id < answer2.kmer_id;
    }

	int MPHFCompare(uint64_t kmer1, uint64_t kmer2)
	{
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
};

#endif