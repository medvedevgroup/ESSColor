#include <iostream>
#include <string>

#include "kmc_file.h"

using namespace std;


int main(int argc, char const *argv[])
{
	// Open a kmc database
	string db_name("data/sarscov2_k31");
	CKMCFile db;

	if (!db.OpenForListing(db_name))
	{
		cerr << "Impossible to open DB " << db_name << endl;
		exit(1);
	}

	// Get metadata
	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
	uint64 _total_kmers, _max_count;

	db.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
	CKmerAPI kmer_object(_kmer_length);

	cout << _total_kmers << endl;

	CKmerAPI current_kmer(_kmer_length), prev_kmer(_kmer_length);
	uint32 counter; string kmer;
	while (db.ReadNextKmer(current_kmer, counter))
	{
		current_kmer.to_string(kmer);
		cout << kmer << " " << counter << " " << (current_kmer > prev_kmer ? ">" : "<") << endl;
		prev_kmer = current_kmer;
	}

	return 0;
}