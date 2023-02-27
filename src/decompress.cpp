//version: feb 6:pausing the combo
#include <cmph.h> //#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sdsl/bit_vectors.hpp>
#include <unordered_map>
#include <sstream>
#include <string>
#include <queue>
#include <map>
#include <climits> // for u_int32_t_BIT
#include <iterator>
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace std;
using namespace sdsl;

uint64_t written_kmer = 0;
#include <unordered_map>


void dbg(){

}
uint64_t CATEGORY_RUN=(uint64_t) 3;
uint64_t CATEGORY_COLCLASS=(uint64_t) 0;
uint64_t CATEGORY_COLVEC=(uint64_t) 2;
uint64_t CATEGORY_COLVEC_ONE = (uint64_t) 4; //100
uint64_t CATEGORY_COLVEC_TWO = (uint64_t) 5; //101

bool DEBUG_MODE = true;

namespace TimeMeasure
{
	double t_begin,t_end; struct timeval timet;
	void time_start(){
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	}
	void time_end(string msg){
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
		cout<<msg<<" time = ";
		printf("%.2fs\n",t_end - t_begin);
	}
} using namespace TimeMeasure;

// namespace CMPH{
// 	cmph_t *hash_cmph = NULL;
// 	void create_table(string filename, int nelem ){
// 		FILE * keys_fd = fopen(filename.c_str(), "r");
		
// 		if (keys_fd == NULL) 
// 		{
// 		fprintf(stderr, "File not found\n");
// 		exit(1);
// 		}	
// 		// Source of keys
// 		cmph_io_adapter_t *source = cmph_io_nlfile_adapter(keys_fd);
	
// 		cmph_config_t *config = cmph_config_new(source);
// 		cmph_config_set_algo(config, CMPH_CHM);
// 		hash_cmph = cmph_new(config);
// 		cmph_config_destroy(config);
		
// 		cmph_io_nlfile_adapter_destroy(source);   
// 		fclose(keys_fd);
// 	}

// 	unsigned int lookup(string str){	
// 		const char *key = str.c_str(); 
// 		//Find key
// 		unsigned int id = cmph_search(hash_cmph, key, (cmph_uint32)strlen(key));
// 		// fprintf(stderr, "Id:%u\n", id);
// 		//Destroy hash
// 		//cmph_destroy(hash);
// 		return id;
// 	}

// 	void mphf_destroy(){
// 		cmph_destroy(hash_cmph);
// 	}
// }
// using namespace CMPH;

class OutputFile{
	public:
		string filename;
		std::ofstream fs;
	OutputFile(){

	}
	OutputFile(string filename){
		this->filename = filename;
		fs.open (filename.c_str(),  std::fstream::out );
	}
	void init(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::out);
	}
	void write(string towrite){
		fs << towrite; // <<endl;
	}
    void close(){
        fs.close();
    }
	~OutputFile(){
		fs.close();
	}
};

class InputFile{
	public:
	string filename;
	std::fstream fs;
	InputFile(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::in);
	}
	InputFile(){
	}
	void init(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::in);
	}
	void rewind(){
		this->fs.close();
		this->fs.open(this->filename, fstream::in);
	}
    void close(){
        fs.close();
    }
	~InputFile(){
		fs.close();
	}
};

class DebugFile : public OutputFile	//derived class
{
	public:
		DebugFile(string filename){
					this->filename = filename;
					fs.open (filename.c_str(),  std::fstream::out );

		}
		DebugFile(){

		}
		void init(const std::string filename)
		{
				this->filename = filename;
				this->fs.open(this->filename, fstream::out);
		}
};

typedef std::vector<bool> HuffCode;
typedef std::map<u_int32_t, HuffCode> HuffCodeMap;

namespace Huffman{
	/// @brief source rosetta code
	class INode
	{
		public:
			const int f;
			virtual ~INode() {}
		protected:
			INode(int f) : f(f) {}
	};
	class InternalNode : public INode
	{
	public:
		INode *const left;
		INode *const right;

		InternalNode(INode* c0, INode* c1) : INode(c0->f + c1->f), left(c0), right(c1) {}
		~InternalNode()
		{
			delete left;
			delete right;
		}
	};
	class LeafNode : public INode
	{
	public:
		u_int32_t c;

		LeafNode(int f, u_int32_t c) : INode(f), c(c) {}
	};

	struct NodeCmp
	{
		bool operator()(const INode* lhs, const INode* rhs) const { return lhs->f > rhs->f; }
	};

	INode* BuildTree(u_int32_t* frequencies, u_int32_t UniqueSymbols)
	{
		std::priority_queue<INode*, std::vector<INode*>, NodeCmp> trees;
	
		for (u_int32_t i = 0; i < UniqueSymbols; ++i)
		{
			if(frequencies[i] != 0)
				trees.push(new LeafNode(frequencies[i], (u_int32_t)i));
		}
		while (trees.size() > 1)
		{
			INode* childR = trees.top();
			trees.pop();

			INode* childL = trees.top();
			trees.pop();

			INode* parent = new InternalNode(childR, childL);
			trees.push(parent);
		}
		return trees.top();
	}

	void GenerateCodes(const INode* node, const HuffCode& prefix, HuffCodeMap& outCodes)
	{
		if (const LeafNode* lf = dynamic_cast<const LeafNode*>(node))
		{
			outCodes[lf->c] = prefix;
		}
		else if (const InternalNode* in = dynamic_cast<const InternalNode*>(node))
		{
			HuffCode leftPrefix = prefix;
			leftPrefix.push_back(false);
			GenerateCodes(in->left, leftPrefix, outCodes);

			HuffCode rightPrefix = prefix;
			rightPrefix.push_back(true);
			GenerateCodes(in->right, rightPrefix, outCodes);
		}
	}

    u_int32_t HuffDecode(const INode *root, string s, int &start_loc)
    {
        string ans = "";
        u_int32_t ansint = 0;
        const INode *curr = root;
        for (int i = start_loc; i < s.size() + 1; i++)
        {
            if (const LeafNode *lf = dynamic_cast<const LeafNode *>(curr))
            {
                ansint = (u_int32_t)(lf->c);
                start_loc = i;
                return ansint;
                curr = root;
            }
            else if (const InternalNode *internal = dynamic_cast<const InternalNode *>(curr))
            {
                if (s[i] == '0')
                    curr = internal->left;
                else
                    curr = internal->right;

                if (const LeafNode *lf = dynamic_cast<const LeafNode *>(curr))
                {
                    ansint = (u_int32_t)(lf->c);
                    // cout<< ansint;
                }
            }
            else
            {
                cout << "FAIL" << endl;
            }
        }
        // cout<<ans<<endl;

        return ansint;
    }

    vector<int> read_l_huff_codes(int l, string s, uint64_t& b_it, INode* root){
        DebugFile logfile_huff_decode("logfile_huff_decode");
		 vector<int> v ;
        int loc  = b_it;
        if(DEBUG_MODE) logfile_huff_decode.fs<<l<<":";
		while(l){
			u_int32_t decoded_col_class = HuffDecode(root, s, loc);
            v.push_back(decoded_col_class);
            if(DEBUG_MODE) logfile_huff_decode.fs<<decoded_col_class<<",";
			l--;
		}
        b_it = loc;
        if(DEBUG_MODE) logfile_huff_decode.fs<<endl;
        return v;
	}
}
using namespace Huffman;
HuffCodeMap huff_code_map;

class COLESS_Decompress{
public:
    uint64_t num_kmers;
    int max_run = 16;
    int lmaxrun = 4;
    
    int C;
    int M;
    int lm, lc;

    OutputFile dec_ess_color;
    InputFile spss_boundary_file;
    
    string str_local;
    u_int64_t b_it_local = 0;

    string str_map;  //reuse
    uint64_t b_it = 0; //reuse

    //global
    string* global_table;
    vector<char> spss_boundary;
    INode* huff_root;

    //per simplitig
    vector<int> local_hash_table;
    int l_of_curr_simplitig;
    char per_simplitig_use_local_id;
    uint64_t per_simplitig_bigD;
    
    bool USE_LOCAL_TABLE = true;
    bool USE_HUFFMAN = true;
    bool ALWAYS_LOCAL_OR_GLOBAL = false;

 

    COLESS_Decompress(uint64_t num_kmers, int M, int C, string spss_boundary_fname, int max_run)
    {
        this->num_kmers = num_kmers;
        this->C = C;
        this->M = M;
        this->max_run = max_run;
        
        lmaxrun = ceil(log2(max_run));
        lm = ceil(log2(M));
        lc = ceil(log2(C));

        global_table = new string[M];

        dec_ess_color.init("dec_ess_color");
        spss_boundary_file.init(spss_boundary_fname);
    }

    void load_rrr_into_string(string rrr_filename, string &where_to_load){
        // rrr_vector<256> rrr_map; // rrr_vector rrr_map = rrr_vector<256>();
        // load_from_file(rrr_map, "rrr_map"); //sdsl namespace
        // rrr_vector<256> rrr_main; // rrr_main =  rrr_vector<256>();
        // load_from_file(rrr_main, "rrr_main"); //sdsl namespace  
        rrr_vector<256> rrr_bv;      
        load_from_file(rrr_bv, rrr_filename);  //sdsl namespace
        // std::ofstream out("str_bv_mapping.txt");
        // out << rrr_map;
        // out.close();
        stringstream buffer;
        buffer << rrr_bv;
        where_to_load = buffer.str();
    }
    void load_bb_into_string(string filename_bb, string &where_to_load){
        InputFile file_bb(filename_bb);//i.e. bb_main
        getline(file_bb.fs, where_to_load);
    }

    INode* build_huff_tree(string freq="frequency_sorted"){
        
		time_start();
		InputFile infile_freq(freq);
		string line;

		// Build frequency table
		u_int32_t *frequencies = new u_int32_t[M]; // M -> no. of unique symbols
		std::fill_n(frequencies, M, 0);
		u_int32_t i = 0;
		while(getline(infile_freq.fs, line)){
			stringstream ss(line);
			u_int32_t a ;
			ss >> a; 
			frequencies[i++]= a;
		}		
		time_end("Read freq for "+to_string(M)+" values.");


		time_start();
		INode* root = BuildTree(frequencies, M);
        DebugFile huff_codes("huff_codes");
        GenerateCodes(root, HuffCode(), huff_code_map); // huff_code_map is filled: uint32t colclassid-> vector bool
		delete frequencies;
        for (HuffCodeMap::const_iterator it = huff_code_map.begin(); it != huff_code_map.end(); ++it)
    {
        if(DEBUG_MODE) huff_codes.fs << it->first << " ";
        std::copy(it->second.begin(), it->second.end(),
                  std::ostream_iterator<bool>(huff_codes.fs));
        if(DEBUG_MODE) huff_codes.fs << std::endl;
    }
		//delete root;
		time_end("Build huffman tree on" +to_string(M)+" values.");
        return root;
    }


    uint64_t convert_binary_string_to_uint(string &str, int start, int end, int block_sz2)
    { // convert_binary_string_to_uint
        uint64_t res = 0;
        int block_sz = end - start + 1;
        // 		assert(block_sz==block_sz2);
        int i = 0;
        for (int64_t j = end; j >= start; j--)
        {
            if (str[j] == '1')
            {
                res |= 1 << i;
            }
            i += 1;
        }
        return res;
    }

    char read_one_bit(string& str, uint64_t& b_it){ //convert_binary_string_to_uint
        return str[b_it++];
    }

    int read_number_encoded_in_unary_zero(string& str, uint64_t& b_it){ //convert_binary_string_to_uint
        int length = 0;
        while(str[b_it++]=='0'){
            length+=1;
        }
        return length;
    }
    int read_number_encoded_in_unary_one(string& str, uint64_t& b_it){ //convert_binary_string_to_uint
        int length = 0;
        while(str[b_it++]=='1'){
            length+=1;
        }
        return length;
    }

    
    string read_color_vector(string& str, uint64_t& b_it){
        string col_vec = str.substr(b_it, C);
        b_it+=C;
        return col_vec;
    }

    void flip_bit(string& s, int pos){
        if(s[C-pos-1] == '1')  {
            s[C-pos-1]='0';
        } else{
            s[C-pos-1]='1';
        }
    }
    uint64_t read_uint(string& str, uint64_t& b_it, int block_sz){ //convert_binary_string_to_uint
        uint64_t res = 0;
        //int block_sz = end - start + 1;
        uint64_t end = block_sz + b_it - 1;
// 		assert(block_sz==block_sz2);
        uint64_t i = 0;
        uint64_t j = end;

        while(true){
            if (str[j]=='1') {
                res |= 1 << i;
            }
            i+=1;
            if(j!=b_it){
                j--;
            }else{
                break;
            }
        }
        b_it += block_sz;
        return res;
    }

    void read_local_hash_table_per_simplitig(string str_local, u_int64_t& b_it){
        per_simplitig_bigD = read_uint(str_local, b_it, 2); //0, 1, 2
        per_simplitig_use_local_id = read_one_bit(str_local, b_it);
        cout<<"curr: " << per_simplitig_bigD<<" "<<per_simplitig_use_local_id<<" ";
        
        if(per_simplitig_use_local_id == '1'){
            l_of_curr_simplitig = read_uint(str_local, b_it, lm);
            //int ll = ceil(log2(l));
            local_hash_table = read_l_huff_codes(l_of_curr_simplitig, str_local, b_it, huff_root); //0->(0,M-1), 1->(0,M-1) ... l*lm bits
            cout<<l_of_curr_simplitig;
        }else{
            
        }
        cout<<endl;
    }

    bool start_of_simplitig(int written_kmer_idx){
        return spss_boundary[written_kmer_idx] == '1';
    }

    bool end_of_simplitig(int written_kmer_idx){
        return spss_boundary[(written_kmer_idx+1)%num_kmers] == '1';
    }

    void flush_skip_and_del(vector<int> &differ_run, string& last_col_vector, OutputFile& dec_ess_color){
        if (differ_run.size())
        {
            for (int d : differ_run)
            {
                flip_bit(last_col_vector, d);
            }
            dec_ess_color.fs << last_col_vector << endl;
            // if(start_of_simplitig(written_kmer)){ 
            //     read_local_hash_table_per_simplitig(str_local, b_it_local);
            // }
            written_kmer+=1;
            
            differ_run.clear();
        }         
    }

    void run()
    {   
        huff_root = build_huff_tree();
        load_bb_into_string("bb_map", str_map);
        b_it = 0;
        DebugFile color_global("color_global"); //M color vectors //DEBUGFILE
        for (int i = 0; i < M; i++)
        {
            string col_vector = read_color_vector(str_map, b_it);
            if(DEBUG_MODE) color_global.fs << col_vector << endl;
            global_table[i] = col_vector;
        }
        color_global.fs.close();

        int num_simplitig = 0;
        // Load SPSS boundary file
        for (uint64_t i=0; i < num_kmers; i+=1){ //load spss_boundary vector in memory from disk
			string spss_line;
			getline (spss_boundary_file.fs,spss_line); 
			spss_boundary.push_back(spss_line[0]); //this kmer starts a simplitig
            if(spss_line[0]=='1'){
                num_simplitig+=1;
            }
		}

        //read local table, 
        load_bb_into_string("bb_local_table", str_local);
        
        // time_start();
        // create_table(color_global.filename, M);
        // time_end("CMPH table create for "+to_string(M)+" keys.");

        load_bb_into_string("bb_main", str_map);
        b_it = 0;
        vector<int> differ_run;
        string last_col_vector = "";
        while (b_it < str_map.length())
        {
            char c = read_one_bit(str_map, b_it);
            if (c == '0')
            {
                flush_skip_and_del(differ_run, last_col_vector,dec_ess_color);
                if(start_of_simplitig(written_kmer)){ 
                    read_local_hash_table_per_simplitig(str_local, b_it_local); //changes l_of_curr_simplitig
                }
                if(per_simplitig_use_local_id == '1'){//using local table
                    int local_id = 0;
                    if(ceil(log2(l_of_curr_simplitig)) != 0){
                        local_id = read_uint(str_map, b_it, ceil(log2(l_of_curr_simplitig)));
                    }
                     
                    int col_class = local_hash_table[local_id]; 
                    last_col_vector = global_table[col_class];
                    dec_ess_color.fs << last_col_vector << endl;
                    written_kmer+=1;
                }else{
                    if(USE_HUFFMAN==false){
                        uint64_t col_class = read_uint(str_map, b_it, lm);
                        last_col_vector = global_table[col_class];
                        dec_ess_color.fs << last_col_vector << endl;
                        written_kmer+=1;
                    }else{
                        //TODO -- untested
                        uint64_t col_class = read_l_huff_codes(1, str_map, b_it, huff_root)[0];
                        last_col_vector = global_table[col_class];
                        dec_ess_color.fs << last_col_vector << endl;
                        written_kmer+=1;
                    }

                }
            }
            else if (c == '1')
            {
                char c2 = '1';
                if(per_simplitig_bigD != 0){
                    c2 = read_one_bit(str_map, b_it);;
                }
                
                if (c2 == '1')
                { // run
                    flush_skip_and_del(differ_run, last_col_vector,dec_ess_color);

                    {//paul method
                        int q = read_number_encoded_in_unary_one(str_map, b_it);
                        assert(read_one_bit(str_map, b_it) == '0');
                        int rem = read_uint(str_map, b_it, lmaxrun);
                        int skip = q * max_run + rem;
                        while (skip)
                        {
                            dec_ess_color.fs << last_col_vector << endl;
                            written_kmer+=1;
                            skip--;
                        }
                    }
                    {//my method
                        // int len_of_logskip = read_number_encoded_in_unary_one(str_map, b_it);
                        // assert(read_one_bit(str_map, b_it) == '0');
                        // int skip = read_uint(str_map, b_it, len_of_logskip);
                    
                        // while (skip)
                        // {
                        //     dec_ess_color.fs << last_col_vector << endl;
                        //     written_kmer+=1;
                        //     skip--;
                        // }
                    }
                }
                else
                { // lc 10
                    flush_skip_and_del(differ_run, last_col_vector,dec_ess_color);
                    if(per_simplitig_bigD == 1){
                        int differing_bit = read_uint(str_map, b_it, lc);
                        differ_run.push_back(differing_bit);
                    }else if(per_simplitig_bigD == 2){
                        char c3 = read_one_bit(str_map, b_it);;
                        if(c3 == '1'){ // 101 // read two
                            int differing_bit = read_uint(str_map, b_it, lc);
                            differ_run.push_back(differing_bit);

                            differing_bit = read_uint(str_map, b_it, lc);
                            differ_run.push_back(differing_bit);
                        }else{//c3 = 0, 100 // read one
                            int differing_bit = read_uint(str_map, b_it, lc);
                            differ_run.push_back(differing_bit);
                        }

                    }
                    
                }
            }
        } 
        flush_skip_and_del(differ_run, last_col_vector,dec_ess_color);
    }
};


int main (int argc, char* argv[]){
	vector<string> args(argv + 1, argv + argc);
    string dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname; //string tmp_dir;
    int M, C;
    int max_run = 16;
	uint64_t num_kmers=0;
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: tool -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run>" << endl;
            return 0;
        } else if (*i == "-i") {
            dedup_bitmatrix_fname = *++i;
        } else if (*i == "-d") {
            dup_bitmatrix_fname = *++i;
        }else if (*i == "-c") {
            C = std::stoi(*++i);
        }else if (*i == "-m") {
            M = std::stoi(*++i);
        }else if (*i == "-k") {
            num_kmers = std::stol(*++i);
        }else if (*i == "-s") {
            spss_boundary_fname = *++i;
		}else if (*i == "-x") {
            max_run = std::stoi(*++i);
		}
		// else if (*i == "-t") {
        //     tmp_dir  = *++i;
		// }
    }

	COLESS_Decompress cdec(num_kmers, M, C,  spss_boundary_fname, max_run);
    cdec.run();
	return EXIT_SUCCESS;
}