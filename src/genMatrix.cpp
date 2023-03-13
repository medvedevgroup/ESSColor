#include <string>
#include <iostream>
#include <vector>

#include "colormatrix.hpp"






int main(int argc, char const *argv[])
{
	std::string db_name("data/sarscov2_k31");
	std::vector<uint64_t> kmers = load_from_file(db_name);

	for (uint64_t kmer : kmers)
		std::cout << kmer << std::endl;

	return 0;
}

