#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "cxxopts.hpp"

#include "colormatrix.hpp"

using namespace std;
using namespace sshash;


#define OPT_PARSING_ERROR 01
#define DB_LIST_ERROR 02




bool test_mandatory_option(const cxxopts::ParseResult& res, const string& opt)
{
	if (res.count(opt) == 0) {
		cerr << opt << " missing option" << endl;
		cerr << "cmd: genmatrix -c <db_list.txt> -o <out_matrix> -l <out_kmer_list>" << endl;
		return false;
	}
	return true;
}


cxxopts::ParseResult parse_args(int argc, char const *argv[])
{
	cxxopts::Options options("genmatrix", "Generates a kmer color matrix from a list of KMC databases. The outputed matrix has one column per database (in the same order that the input list). The rows of the matrix are kmers. The order of the kmers is the same than the KMC order");

	options.add_options()
		("c,count-list", "[Mandatory] Path to KMC database files. One line per database", cxxopts::value<string>())
		("d,debug-verif", "Debug flag to verify if the output coresponds to the input (Time consuming).", cxxopts::value<bool>()->default_value("false"))
		("o,outmatrix", "[Mandatory] Path to the output color matrix", cxxopts::value<string>())
		("l,spss", "[Mandatory] Path to the corresponding union SPSS", cxxopts::value<string>())
		("s,strout", "String output", cxxopts::value<bool>()->default_value("false"))
		("k,kmer-size", "[Mandatory] k-mer size", cxxopts::value<uint64_t>())
		;

	cxxopts::ParseResult res = options.parse(argc, argv);
	if ((not test_mandatory_option(res, "count-list"))
		or (not test_mandatory_option(res, "outmatrix"))
		or (not test_mandatory_option(res, "spss"))
		or (not test_mandatory_option(res, "kmer-size"))
	) {
		cerr << endl << options.help() << endl;
		exit(OPT_PARSING_ERROR);
	}

	return res;
}


vector<string> get_databases(const string& list_file)
{
	vector<string> db_list;

	ifstream file(list_file);
	if (file.is_open())
	{
	    string line;
	    while (getline(file, line)) {
	    	if (line.length() > 0 and line[0] != '#')
		        db_list.push_back(line);
	    }
	    file.close();
	}
	// No file list opened
	else
	{
		cerr << "File " << list_file << " cannot be opened. Please verify the path." << endl;
		exit(DB_LIST_ERROR);
	}

	return db_list;
}


void verif(KmerMatrix & matrix, vector<string> db_list)
{
	cerr << "----- Start kmer presence checking -----" << endl;
	uint64_t missing_kmers = 0;
	uint64_t overpresent_kmers = 0;

	uint64_t dataset_idx = 0;
	uint64_t uint_idx = 0;
	uint64_t bit_idx = 0;
	for (string& db_name : db_list)
	{
		cerr << "* Checking file " << db_name << endl;
		uint64_t db_missing_kmers = 0;
		uint64_t db_overpresent_kmers = 0;

		uint64_t k;
		vector<uint64_t> kmers = load_from_file(db_name, k);

		for (uint64_t kmer: kmers)
		{
			uint64_t idx = matrix.mphf_get_kmer_id(kmer);
			if (idx >= matrix.num_kmers)
			{
				missing_kmers += 1;
				continue;
			}

			uint64_t subvector = matrix.colors[matrix.numbers_per_row * idx + uint_idx];
			uint64_t presence = (subvector >> bit_idx) & 0b1;
			if (presence == 0)
			{
				missing_kmers += 1;
				cerr << kmer2str(kmer, matrix.k) << " in dataset " << db_name << " absent from the matrix" << endl;
			}
		}

		// Update values for the next dataset
		dataset_idx += 1;
		bit_idx += 1;
		if (bit_idx == 64)
		{
			uint_idx += 1;
			bit_idx = 0;
		}
	}

	if (missing_kmers == 0 and overpresent_kmers == 0)
	{
		cerr << "All the kmers are present in the matrix" << endl;
	}
	else
	{
		cerr << (missing_kmers + overpresent_kmers) << " errors detected!" << endl;
		cerr << "  * " << missing_kmers << " missing kmers in the matrix" << endl;
		cerr << "  * " << overpresent_kmers << " present kmers while it should not" << endl;
	}
}




void build_mphf(uint64_t k, dictionary * dict, string spss_file) {
	//write_ess_boundary(spss_file, k);
	auto m = (int) k/2;
	if(m<10) m = k;
	//dictionary dict;

	build_configuration build_config;
	build_config.k = k;
	build_config.m = m;

	build_config.canonical_parsing = true;
	build_config.verbose = true;
	build_config.print();

	dict->build(spss_file, build_config);
	assert(dict->k() == k);
	std::cout<<"dict built complete";
}

int main(int argc, char const *argv[])
{
	auto args = parse_args(argc, argv);
	vector<string> db_list = get_databases(args["count-list"].as<string>());

	// Create the sshash dictionnary
	dictionary * dict = new dictionary();
	build_mphf(args["kmer-size"].as<uint64_t>(),dict, args["spss"].as<string>());

	// Creates a global matrix
	KmerMatrix global_matrix(args["k"].as<uint64_t>(), dict->size()/**/, db_list.size(), dict);

	// For each file, add it to the matrix
	int db_idx = 1;
	uint64_t k = 0;
	for (const string& db_name : db_list)
	{
		cout << "processing database " << db_idx++ << "/" << db_list.size() << endl;
		vector<uint64_t> kmers = load_from_file(db_name, k);
		global_matrix.add_dataset(kmers);
	}

	cout << "Writing the matrix on disk..." << endl;
	
	if (args["strout"].as<bool>())
	{
		global_matrix.to_color_string_file(args["outmatrix"].as<string>());
	}
	else
	{
		global_matrix.to_color_binary_file(args["outmatrix"].as<string>());
	}

	if (args["debug-verif"].as<bool>())
	{
		verif(global_matrix, db_list);
	}

	delete dict;

	return 0;
}