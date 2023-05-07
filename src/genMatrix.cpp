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


#include "gz/zip_stream.hpp"

namespace sshash {
//@OVERRIDE


void parse_file(std::istream& is, parse_data& data, build_configuration const& build_config) {
	cout<<"OVERRIDE"<<endl;
    uint64_t k = build_config.k;
    uint64_t m = build_config.m;
    uint64_t seed = build_config.seed;
    uint64_t max_num_kmers_in_super_kmer = k - m + 1;
    uint64_t block_size = 2 * k - m;  // max_num_kmers_in_super_kmer + k - 1

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    /* fit into the wanted number of bits */
    assert(max_num_kmers_in_super_kmer < (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));

    compact_string_pool::builder builder(k);

    std::string sequence;
    uint64_t prev_minimizer = constants::invalid_uint64;

    uint64_t begin = 0;  // begin of parsed super_kmer in sequence
    uint64_t end = 0;    // end of parsed super_kmer in sequence
    uint64_t num_sequences = 0;
    uint64_t num_bases = 0;
    bool glue = false;

    auto append_super_kmer = [&]() {
        if (sequence.empty() or prev_minimizer == constants::invalid_uint64 or begin == end) {
            return;
        }

        assert(end > begin);
        char const* super_kmer = sequence.data() + begin;
        uint64_t size = (end - begin) + k - 1;
        assert(util::is_valid(super_kmer, size));

        /* if num_kmers_in_super_kmer > k - m + 1, then split the super_kmer into blocks */
        uint64_t num_kmers_in_super_kmer = end - begin;
        uint64_t num_blocks = num_kmers_in_super_kmer / max_num_kmers_in_super_kmer +
                              (num_kmers_in_super_kmer % max_num_kmers_in_super_kmer != 0);
        assert(num_blocks > 0);
        for (uint64_t i = 0; i != num_blocks; ++i) {
            uint64_t n = block_size;
            if (i == num_blocks - 1) n = size;
            uint64_t num_kmers_in_block = n - k + 1;
            assert(num_kmers_in_block <= max_num_kmers_in_super_kmer);
            data.minimizers.emplace_back(prev_minimizer, builder.offset, num_kmers_in_block);
            builder.append(super_kmer + i * max_num_kmers_in_super_kmer, n, glue);
            if (glue) {
                assert(data.minimizers.back().offset > k - 1);
                data.minimizers.back().offset -= k - 1;
            }
            size -= max_num_kmers_in_super_kmer;
            glue = true;
        }
    };

    uint64_t seq_len = 0;
    uint64_t sum_of_weights = 0;
    data.weights_builder.init();

    /* intervals of weights */
    uint64_t weight_value = constants::invalid_uint64;
    uint64_t weight_length = 0;

    auto parse_header = [&]() {
        if (sequence.empty()) return;

        /*
            Heder format:
            >[id] LN:i:[seq_len] ab:Z:[weight_seq]
            where [weight_seq] is a space-separated sequence of integer counters (the weights),
            whose length is equal to [seq_len]-k+1
        */

        // example header: '>12 LN:i:41 ab:Z:2 2 2 2 2 2 2 2 2 2 2'

        expect(sequence[0], '>');
        uint64_t i = 0;
        i = sequence.find_first_of(' ', i);
        if (i == std::string::npos) throw parse_runtime_error();

        i += 1;
        expect(sequence[i + 0], 'L');
        expect(sequence[i + 1], 'N');
        expect(sequence[i + 2], ':');
        expect(sequence[i + 3], 'i');
        expect(sequence[i + 4], ':');
        i += 5;
        uint64_t j = sequence.find_first_of(' ', i);
        if (j == std::string::npos) throw parse_runtime_error();

        seq_len = std::strtoull(sequence.data() + i, nullptr, 10);
        i = j + 1;
        expect(sequence[i + 0], 'a');
        expect(sequence[i + 1], 'b');
        expect(sequence[i + 2], ':');
        expect(sequence[i + 3], 'Z');
        expect(sequence[i + 4], ':');
        i += 5;

        for (uint64_t j = 0; j != seq_len - k + 1; ++j) {
            uint64_t weight = std::strtoull(sequence.data() + i, nullptr, 10);
            i = sequence.find_first_of(' ', i) + 1;

            data.weights_builder.eat(weight);
            sum_of_weights += weight;

            if (weight == weight_value) {
                weight_length += 1;
            } else {
                if (weight_value != constants::invalid_uint64) {
                    data.weights_builder.push_weight_interval(weight_value, weight_length);
                }
                weight_value = weight;
                weight_length = 1;
            }
        }
    };

    while (!is.eof()) {
        std::getline(is, sequence);  // header sequence
        if (build_config.weighted) parse_header();

        std::getline(is, sequence);  // DNA sequence
        if (sequence.size() < k) continue;

        if (++num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                      << data.num_kmers << " kmers" << std::endl;
        }

        begin = 0;
        end = 0;
        glue = false;  // start a new piece
        prev_minimizer = constants::invalid_uint64;
        num_bases += sequence.size();

        if (build_config.weighted and seq_len != sequence.size()) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << sequence.size() << std::endl;
            throw std::runtime_error("file is malformed");
        }

        while (end != sequence.size() - k + 1) {
            char const* kmer = sequence.data() + end;
            assert(util::is_valid(kmer, k));
            kmer_t uint_kmer = util::string_to_uint_kmer_no_reverse(kmer, k);
            uint64_t minimizer = util::compute_minimizer(uint_kmer, k, m, seed);

            if (build_config.canonical_parsing) {
                kmer_t uint_kmer_rc = util::compute_reverse_complement(uint_kmer, k);
                uint64_t minimizer_rc = util::compute_minimizer(uint_kmer_rc, k, m, seed);
                minimizer = std::min<uint64_t>(minimizer, minimizer_rc);
            }

            if (prev_minimizer == constants::invalid_uint64) prev_minimizer = minimizer;
            if (minimizer != prev_minimizer) {
                append_super_kmer();
                begin = end;
                prev_minimizer = minimizer;
                glue = true;
            }

            ++data.num_kmers;
            ++end;
        }

        append_super_kmer();
    }

    data.minimizers.finalize();
    builder.finalize();
    builder.build(data.strings);

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
              << data.num_kmers << " kmers" << std::endl;
    std::cout << "num_kmers " << data.num_kmers << std::endl;
    std::cout << "num_super_kmers " << data.strings.num_super_kmers() << std::endl;
    std::cout << "num_pieces " << data.strings.pieces.size() << " (+"
              << (2.0 * data.strings.pieces.size() * (k - 1)) / data.num_kmers << " [bits/kmer])"
              << std::endl;
    assert(data.strings.pieces.size() == num_sequences + 1);

    if (build_config.weighted) {
        std::cout << "sum_of_weights " << sum_of_weights << std::endl;
        data.weights_builder.push_weight_interval(weight_value, weight_length);
        data.weights_builder.finalize(data.num_kmers);
    }
}


}  // namespace sshash



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