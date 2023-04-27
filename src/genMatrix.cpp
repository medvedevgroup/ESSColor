#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "cxxopts.hpp"

#include "colormatrix.hpp"

using namespace std;

#include "common.hpp"
#include "bench_utils.hpp"
#include "check_utils.hpp"
// #include "lib/sshash/src/build.cpp"
// #include "lib/sshash/src/query.cpp"
#include "permute.cpp"

using namespace sshash;


#define OPT_PARSING_ERROR 01
#define DB_LIST_ERROR 02


using namespace std;



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
		("l,kmer-list", "[Mandatory] Path to the output ordered kmer list file", cxxopts::value<string>())
		("s,strout", "String output", cxxopts::value<bool>()->default_value("false"))
		;

	cxxopts::ParseResult res = options.parse(argc, argv);
	if ((not test_mandatory_option(res, "count-list"))
		or (not test_mandatory_option(res, "outmatrix"))
		or (not test_mandatory_option(res, "kmer-list"))
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
	for (string& db_name : db_list)
	{
		cerr << "* Checking file " << db_name << endl;
		uint64_t db_missing_kmers = 0;
		uint64_t db_overpresent_kmers = 0;

		uint64_t k;
		vector<uint64_t> kmers = load_from_file(db_name, k);
		vector<uint64_t> row;

		// Looks for all the kmers from the db
		uint64_t list_idx=0, matrix_idx=0;
		while(list_idx < kmers.size() and matrix_idx < matrix.kmers.size())
		{
			// Missing kmer from the db
			if (kmers[list_idx] < matrix.kmers[matrix_idx])
			{
				cerr << kmer2str(kmers[list_idx], k) << " not in matrix while in db" << endl;
				db_missing_kmers += 1;
				list_idx += 1;
				// exit(1);
			}
			// Kmer from the matrix that is not present inside of the current db
			else
			{
				// Get the current_row
				matrix.get_row(matrix_idx, row);
				uint64_t sub_row = row[dataset_idx / 64];
				bool is_present = 1 == ((sub_row >> (dataset_idx % 64)) & 0b1);

				if (kmers[list_idx] > matrix.kmers[matrix_idx])
				{
					if (is_present)
					{
						db_overpresent_kmers += 1;
						cerr << kmer2str(matrix.kmers[matrix_idx], k) << " in matrix while not in db" << endl;
						// exit(1);
					}
					matrix_idx += 1;
				}
				else
				{
					if (not is_present)
					{
						db_missing_kmers += 1;
						cerr << kmer2str(kmers[list_idx], k) << " not in matrix row while in db" << endl;
						// exit(1);
					}
					matrix_idx += 1;
					list_idx += 1;
				}
			}
		}

		// Potential remaining kmers after full consumption of the matrix
		if (list_idx < kmers.size())
		{
			for (uint64_t i(list_idx) ; i<kmers.size() ; i++)
			{
				cerr << kmer2str(kmers[i], k) << " not in matrix while in db" << endl;
				db_missing_kmers += 1;
			}
		}

		// Error stats
		missing_kmers += db_missing_kmers;
		overpresent_kmers += db_overpresent_kmers;
		dataset_idx += 1;
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


int main(int argc, char const *argv[])
{
	auto args = parse_args(argc, argv);
	vector<string> db_list = get_databases(args["count-list"].as<string>());

	MPHFComparator mphfcomparator;
	CascadingMergingMatrix cmm(0.9);
	

	// For each file, add it to the cascading matrix constructor
	int db_idx = 1;
	uint64_t k = 0;
	for (const string& db_name : db_list)
	{
		cout << "processing database " << db_idx++ << "/" << db_list.size() << endl;
		vector<uint64_t> kmers = load_from_file(db_name, k, mphfcomparator);
		KmerMatrix matrix(kmers, k, mphfcomparator);
		cmm.add_matrix(matrix);
	}

	// Get the final matrix after the merging
	cout << "Finalization of the matrix..." << endl;
	KmerMatrix& final_matrix = cmm.get_matrix();
	cout << "Writing the matrix on disk..." << endl;
	if (args["strout"].as<bool>())
	{
		final_matrix.to_color_string_file(args["outmatrix"].as<string>());
		final_matrix.to_kmer_string_file(args["kmer-list"].as<string>());
	}
	else
	{
		final_matrix.to_color_binary_file(args["outmatrix"].as<string>());
		final_matrix.to_kmer_binary_file(args["kmer-list"].as<string>());
	}

	if (args["debug-verif"].as<bool>())
	{
		verif(final_matrix, db_list);
	}

	return 0;
}
