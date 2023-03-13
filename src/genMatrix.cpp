#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "cxxopts.hpp"

#include "colormatrix.hpp"


#define OPT_PARSING_ERROR 01
#define DB_LIST_ERROR 02




bool test_mandatory_option(const cxxopts::ParseResult& res, const std::string& opt)
{
	if (res.count(opt) == 0) {
		std::cerr << opt << " missing option" << std::endl;
		std::cerr << "cmd: genmatrix -c <bd_list.txt>" << std::endl;
		return false;
	}
	return true;
}


cxxopts::ParseResult parse_args(int argc, char const *argv[])
{
	cxxopts::Options options("genmatrix", "Generates a kmer color matrix from a list of KMC databases. The outputed matrix has one column per database (in the same order that the input list). The rows of the matrix are kmers. The order of the kmers is the same than the KMC order");

	options.add_options()
		("c,count-list", "[Mandatory] Path to KMC database files. One line per database", cxxopts::value<std::string>())
		("d,debug-verif", "Debug flag to verify if the output coresponds to the input (Time consuming).", cxxopts::value<bool>()->default_value("false"))
		;

	cxxopts::ParseResult res = options.parse(argc, argv);
	if (not test_mandatory_option(res, "count-list")) {
		std::cerr << std::endl << options.help() << std::endl;
		exit(OPT_PARSING_ERROR);
	}

	return res;
}


std::vector<std::string> get_databases(const std::string& list_file)
{
	std::vector<std::string> db_list;

	std::ifstream file(list_file);
	if (file.is_open())
	{
	    std::string line;
	    while (std::getline(file, line)) {
	    	if (line.length() > 0)
		        std::cout <<  line << std::endl;
	    }
	    file.close();
	}
	// No file list opened
	else
	{
		std::cerr << "File " << list_file << " cannot be opened. Please verify the path." << std::endl;
		exit(DB_LIST_ERROR);
	}

	return db_list;
}


int main(int argc, char const *argv[])
{
	auto args = parse_args(argc, argv);
	std::vector<std::string> db_list = get_databases(args["count-list"].as<std::string>());

	std::string db_name("data/sarscov2_k31");
	std::vector<uint64_t> kmers = load_from_file(db_name);

	// for (uint64_t kmer : kmers)
	// 	std::cout << kmer << std::endl;

	return 0;
}

