#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "cxxopts.hpp"

#include "colormatrix.hpp"

using namespace std;


#define OPT_PARSING_ERROR 01
#define DB_LIST_ERROR 02




bool test_mandatory_option(const cxxopts::ParseResult& res, const string& opt)
{
	if (res.count(opt) == 0) {
		cerr << opt << " missing option" << endl;
		cerr << "cmd: genmatrix -c <bd_list.txt>" << endl;
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
		("s,strout", "String output", cxxopts::value<bool>()->default_value("false"))
		;

	cxxopts::ParseResult res = options.parse(argc, argv);
	if ((not test_mandatory_option(res, "count-list"))
		or (not test_mandatory_option(res, "outmatrix"))
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

}


int main(int argc, char const *argv[])
{
	auto args = parse_args(argc, argv);
	vector<string> db_list = get_databases(args["count-list"].as<string>());

	CascadingMergingMatrix cmm(0.9);

	// For each file, add it to the cascading matrix constructor
	int db_idx=1;
	for (const string& db_name : db_list)
	{
		cout << "processing database " << db_idx++ << "/" << db_list.size() << endl;
		vector<uint64_t> kmers = load_from_file(db_name);
		KmerMatrix matrix(kmers);
		cmm.add_matrix(matrix);
	}

	// Get the final matrix after the merging
	cout << "Finalization of the matrix..." << endl;
	KmerMatrix& final_matrix = cmm.get_matrix();
	if (args["strout"].as<bool>())
		final_matrix.to_color_string_file(args["outmatrix"].as<string>());
	else
		final_matrix.to_color_binary_file(args["outmatrix"].as<string>());


	if (args["debug-verif"].as<bool>())
	{
		verif(final_matrix, db_list);
	}

	return 0;
}

