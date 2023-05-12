//version: mar 1: trying to fix gut
#define VERSION_NAME "APR11,FIXING_GUT DECOM; stream read"
#include<cmph.h> //#include "BooPHF.h"
#include<csignal>
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
#include<unordered_map>
#include<sstream>
#include<string>
#include <queue>
#include <map>
#include <climits> // for u_int32_t_BIT
#include <iterator>
#include <algorithm>
#include <iostream>
#include<sstream>

#include<unistd.h>

using namespace std;
using namespace sdsl;

const int MAX_BUFFER_STRING=2048;
uint64_t written_kmer = 0;
#include <unordered_map>

uint64_t CATEGORY_RUN=(uint64_t) 3;
uint64_t CATEGORY_COLCLASS=(uint64_t) 0;
uint64_t CATEGORY_COLVEC=(uint64_t) 2;
uint64_t CATEGORY_COLVEC_ONE = (uint64_t) 4; //100
uint64_t CATEGORY_COLVEC_TWO = (uint64_t) 5; //101


bool TESTING_SPEED=false;
bool DEBUG_MODE = false;

namespace TimeMeasure
{
	double t_begin,t_end; struct timeval timet;
	void time_start(){
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	}
	void time_end(string msg){
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
		cout<<msg<<" time = ";
		printf("%.8fs\n",t_end - t_begin);
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
}
using namespace Huffman;
HuffCodeMap huff_code_map;


class BlockStream{
    public:

        std::fstream fs;
        uint64_t b_it;
        string str;
        string filename;
         char str_c[MAX_BUFFER_STRING] = "";
    
    BlockStream(string filename){
       
        str = "";
        b_it = 0;
        //fs as file input
        this->filename = filename;
        fs.open(filename, std::fstream::in);
        //fs >> setw(MAX_BUFFER_STRING) >> str;

        //         char b[3] = "";
        // ifstream f("prad.txt");

        str_c[0] = '\0';
        fs.read(str_c, sizeof(str_c)); // Read one less that sizeof(b) to ensure null
        // istream& get (char* s, streamsize n );
        string s2(str_c);
        str = s2;
        //cout<<str<<endl;
        //std::cout.write(reinterpret_cast<char*>(&str_c), sizeof(str_c));
        cout<<"Block size: "<<MAX_BUFFER_STRING<<endl;
    }

    ~BlockStream(){
        fs.close();
    }
    bool end_reached(){
        if(b_it==0 && str.length()==0){
            return true;
        }
        if(str.length()< MAX_BUFFER_STRING && b_it == str.length()){
            return true;
        }else{
            return false;
        }
    }
    void load_string_if_max_exceeds(){
        //char str_c[MAX_BUFFER_STRING] = "";

        if(b_it >= MAX_BUFFER_STRING){
            str = "";
            fs.read(str_c, sizeof(str_c)); // Read one less that sizeof(b) to ensure null
            string s2(str_c);
            str = s2;
            //fs >> setw(MAX_BUFFER_STRING) >> str;
            // if(str.length()< MAX_BUFFER_STRING){
            //     end_reached = true;
            // }
            b_it = 0;
            //cout<<str<<endl;
             
             //std::cout.write(reinterpret_cast<char*>(&str_c), sizeof(str_c));

        }
    }
    

    uint64_t read_uint(int block_sz){ //convert_binary_string_to_uint
        char* str_copy = new char[block_sz+1];
        str_copy[block_sz]='\0';
        uint64_t b_it_copy = b_it;
        for(int i = 0 ; i<block_sz; i++){
            str_copy[i] = str_c[b_it++];
            load_string_if_max_exceeds();
        }
        //cout<<"uint:";
       //      std::cout.write(reinterpret_cast<char*>(&str_copy), sizeof(str_copy));

        uint64_t res = 0;
        //int block_sz = end - start + 1;
        uint64_t end = block_sz - 1;
// 		assert(block_sz==block_sz2);
        uint64_t i = 0;
        uint64_t j = end;


        while(true){
            if (str_copy[j]=='1') {
                res |= 1 << i;
            }
            i+=1;
            if(j!=0){
                j--;
            }else{
                break;
            }
        }

        return res;
    }

    char peek(){
        if(b_it < str.length()){
            //cout<<"peeking str:"<< str_c[b_it]<<endl;
            return str_c[b_it];
        }else if(b_it == str.length() ){
            //cout<<"peeking fs:"<< fs.peek()<<endl;
            return fs.peek();
        }else{
            cout<<"ERROR"<<endl;
            exit(5);
        }
        
    }

    char read_one_bit(){ //convert_binary_string_to_uint
        char toreturn = str_c[b_it];
        b_it++;
        load_string_if_max_exceeds();
        return toreturn;
    }
    u_int32_t HuffDecode(const INode *root)
    {
        string ans = "";
        u_int32_t ansint = 0;
        const INode *curr = root;

        while(true)
        {
            if (const LeafNode *lf = dynamic_cast<const LeafNode *>(curr))
            {
                ansint = (u_int32_t)(lf->c);
                return ansint;
                curr = root;
            }
            else if (const InternalNode *internal = dynamic_cast<const InternalNode *>(curr))
            {
                if (str_c[b_it++] == '0')
                    curr = internal->left;
                else
                    curr = internal->right;
                load_string_if_max_exceeds();

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

    vector<int> read_l_huff_codes(INode* root, int l){
        vector<int> huffcodes;
        while(l){
            u_int32_t decoded_col_class = HuffDecode(root);
            huffcodes.push_back(decoded_col_class);
            l--;
        }
        return huffcodes;
	}


    int read_number_encoded_in_unary_zero(){ 
        int length = 0;
        while(true){
            if(str_c[b_it]=='0'){
                length+=1;
                b_it+=1;
                load_string_if_max_exceeds();
            }else{
                break;
            }
        }
        return length;
    }

    int read_number_encoded_in_unary_one(){ 
        int length = 0;
        while(true){
            if(str_c[b_it]=='1'){
                length+=1;
                b_it+=1;
                load_string_if_max_exceeds();
            }else{
                break;
            }
        }
        return length;
    }

    
    void read_string_of_length(char* col_vec, int C){
        int i = 0;
        while(C){
            col_vec[i++] = str_c[b_it];
            b_it++;
            load_string_if_max_exceeds();
            C--;
        }
        return;
    }
};




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

        //global_table = new string[M];

        dec_ess_color.init("dec_ess_color");
        spss_boundary_file.init(spss_boundary_fname);
    }
    COLESS_Decompress()
    
    {
        max_run = 16;
        dec_ess_color.init("dec_ess_color");
    }
    void dump_rrr_into_bb(string rrr_filename, string bb_file){
        rrr_vector<256> rrr_bv;      
        load_from_file(rrr_bv, rrr_filename);  //sdsl namespace
        std::ofstream out(bb_file);
        out << rrr_bv;
        out.close();
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
        // std::string result(50, '\0');
        // if (!inss.read(&result[0], result.size()))
        // throw std::runtime_error("Could not read enough characters.\n");



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


    void get_ess_boundary_vector_from_pos(std::vector<uint64_t>& positions, uint64_t bv_size, vector<char> &vec){
        //sort(positions.begin(), positions.end());
        uint64_t bvi = 0;
        for (uint64_t k = 0; k<bv_size; k++){
            if(bvi < positions.size()){
                if(positions[bvi]==k){
                    vec[k]='1'
                    bvi++;
                }else{
                    vec[k]='0';
                }
            }else{
                vec[k]='0';
            }	
        }
    }
    void get_ess_boundary_from_essd(string megaessd, int num_simplitig, int kmer_size, uint64_t num_kmers, vector<char> ess_boundary_vector ){
        ifstream is(megaessd);
        ess_boundary_vector.push_back(0);
        num_kmers = 0;
        num_simplitig=0;
        
        uint64_t last_written_pos = 0;
        string sequence;
        while (!is.eof()) {
            std::getline(is, sequence);  // header sequence
            std::getline(is, sequence);  // DNA sequence
            uint64_t numkmer_in_simplitig = sequence.size() - kmer_size + 1;
            last_written_pos += numkmer_in_simplitig;
            ess_boundary_vector.push_back(last_written_pos);
            num_simplitig+=1;
        }
        uint64_t boundary_num_kmer = ess_boundary_vector.back();
        num_kmers = boundary_num_kmer;
        ess_boundary_vector.pop_back();
        vector<char> vec(boundary_num_kmer);
        get_ess_boundary_vector_from_pos(ess_boundary_vector, boundary_num_kmer, vec);
    }
    void read_meta(string file_meta, int& kmer_size, int& num_colors){
        ifstream is(file_meta);
        string line;
        getline(is, line);
        kmer_size = stoi(line);
        num_colors = 0;
        while(getline(is, line)){
            num_colors++;
        }
        is.close();
    }


    void read_local_hash_table_per_simplitig(BlockStream& bs_local){
        per_simplitig_bigD = bs_local.read_uint(2); //0, 1, 2 
        per_simplitig_use_local_id = bs_local.read_one_bit();
        
        if(DEBUG_MODE) cout<<"curr: " << per_simplitig_bigD<<" "<<per_simplitig_use_local_id<<" ";
        
        if(per_simplitig_use_local_id == '1'){
            l_of_curr_simplitig = bs_local.read_uint(lm);
            //int ll = ceil(log2(l));
            local_hash_table = bs_local.read_l_huff_codes(huff_root, l_of_curr_simplitig); //0->(0,M-1), 1->(0,M-1) ... l*lm bits
            if(DEBUG_MODE)  cout<<l_of_curr_simplitig;
        }else{
            
        }
        if(DEBUG_MODE)  cout<<endl;
    }

    bool start_of_simplitig(uint64_t written_kmer_idx){
        return spss_boundary[written_kmer_idx] == '1';
    }

    bool end_of_simplitig(uint64_t written_kmer_idx){
        return spss_boundary[(written_kmer_idx+1)%num_kmers] == '1';
    }

    void flush_skip_and_del(vector<int> &differ_run, string& last_col_vector, OutputFile& dec_ess_color){
        if (differ_run.size())
        {
            for (int d : differ_run)
            {
                flip_bit(last_col_vector, d);
            }
            if(!TESTING_SPEED) dec_ess_color.fs << last_col_vector << endl;
            

            // if(start_of_simplitig(written_kmer)){ 
            //     read_local_hash_table_per_simplitig(str_local, b_it_local);
            // }
            written_kmer+=1;
            
            differ_run.clear();
        }         
    }

    // void test_run(){
    //     rrr_vector<256> rrr_bv; 
    //     time_start();     
    //     load_from_file(rrr_bv, "rrr_main");  //sdsl namespace
    //     b_it = 0;
    //      for (int i = 0; i < M; i++)
    //     {
    //         read_color_vector_binary(rrr_bv, b_it);
    //     }
    //     time_end("binary read");
    // }

    void flip_bit(string& s, int pos){
        if(s[C-pos-1] == '1')  {
            s[C-pos-1]='0';
        } else{
            s[C-pos-1]='1';
        }
    }
    
    // void read_from_stream(std::fstream& fs, int block_sz, string& str, uint64_t& b_it){
    //     //str is always a string of length bsz
    //     // at the end, b_it marks the loc of str from where it should be read
    //     if(b_it==0){
    //         fs >> setw(block_sz) >> str;
    //     }else{
    //         string str1, str2;
    //         if(remainder_str.length()-b_it >= block_sz){//remainder str is large enough
    //            // read fully from remainder, update b_it, remainder_str
    //            str = remainder_str.substr(b_it, block_sz);
    //            remainder_str = remainder_str.substr(b_it+block_sz, remainder_str.length()-b_it-block_sz);
    //            b_it += block_sz;
               
    //         }else{
    //             //read first part from remainder
    //             str1 = remainder_str.substr(b_it, remainder_str.length()-b_it);
               
    //             //read second part from fs 
    //             fs >> setw(block_sz - remainder_str.length() + b_it ) >> str2;
    //             str = str1+str2;
    //             b_it=0;
    //         }
    //     }
    // }

    void run()
    {   
        dump_rrr_into_bb("rrr_main", "bb_main");
        dump_rrr_into_bb("rrr_map", "bb_map");
        dump_rrr_into_bb("rrr_local_table", "bb_local_table");
       

        BlockStream bs_main("bb_main");
        time_start();
        huff_root = build_huff_tree();
        DebugFile color_global("color_global"); //M color vectors //DEBUGFILE
        
        bool USING_NONMST=true;
        if(!USING_NONMST){

            BlockStream bs_map("bb_map");
            global_table = new string[M];
            for (int i = 0; i < M; i++)
            {
                char* col_vector = new char[C+1];
                col_vector[C] = '\0';
                bs_map.read_string_of_length(col_vector, C);
                if(DEBUG_MODE) {
                    color_global.fs << col_vector << endl;
                    //color_global.fs.write(reinterpret_cast<char*>(&col_vector), sizeof(col_vector));
                }
                global_table[i] = string(col_vector);
                delete col_vector;
            }
            color_global.fs.close();
            cout<<"Global table done."<<endl;
        }
        if(USING_NONMST){
            int lc =ceil(log2(C));
            rrr_vector<256> rrr_map_hd_boundary;      
            load_from_file(rrr_map_hd_boundary, "rrr_map_hd_boundary");  //sdsl namespace
    
            string last_colvector(C, '0');
            //dump_rrr_into_bb("rrr_map_hd_boundary", "bb_map_hd_boundary");
            dump_rrr_into_bb("rrr_map_hd", "bb_map_hd");
            BlockStream bs_map_hd("bb_map_hd");
    
            size_t ones = rrr_vector<256>::rank_1_type(&rrr_map_hd_boundary)(rrr_map_hd.size()); 
            rrr_vector<256>::select_1_type rrr_hd_sel(&rrr_map_hd_boundary);
            M = ones;
            global_table = new string[M];
            lm = ceil(log2(M));
            lc = ceil(log2(C));
            
            size_t prev_begin = 0;
            for (size_t i=1; i <= M; ++i){
                size_t prev_end = rrr_hd_sel(i);
                size_t block_len = prev_end - prev_begin + 1;
                size_t numblocks = (block_len / lc);
                //howmanydelta = (pos+1)/lc;
                vector<int> flip_loc;
                while(numblocks){
                    uint64_t fliploc = bs_map_hd.read_uint(lc);
                    if(last_colvector[fliploc]=='1'){
                        last_colvector[fliploc] = '0';
                    }else{
                        last_colvector[fliploc] = '1';
                    }
                    numblocks--;
                }
                global_table[i] = last_colvector;
                if(DEBUG_MODE) {
                    color_global.fs << last_colvector << endl;
                }
                prev_begin = prev_end+1;
            }
            cout<<"Global table done (NONMST)."<<endl;

        }
        //exit(0);
        int num_simplitig = 0;
        // Load SPSS boundary file
        time_start();
        // for (uint64_t i=0; i < num_kmers; i+=1){ //load spss_boundary vector in memory from disk
		// 	string spss_line;
		// 	getline (spss_boundary_file.fs,spss_line); 
		// 	spss_boundary.push_back(spss_line[0]); //this kmer starts a simplitig
        //     if(spss_line[0]=='1'){
        //         num_simplitig+=1;
        //     }
		// }
        int kmer_size;
        read_meta("meta.txt", kmer_size, C);
        get_ess_boundary_from_essd("mega.essd", num_simplitigs, kmer_size, num_kmers,spss_boundary);
        time_end("SPSS boundary read "+to_string(num_kmers)+" bits.");

        //read local table, 
        BlockStream bs_local("bb_local_table");

        // time_start();
        // create_table(color_global.filename, M);
        // time_end("CMPH table create for "+to_string(M)+" keys.");

        vector<int> differ_run;
        string last_col_vector = "";
        while(written_kmer < num_kmers )
        {
            char c = bs_main.read_one_bit();
            if (c == '0')
            {
                flush_skip_and_del(differ_run, last_col_vector,dec_ess_color);
                if(start_of_simplitig(written_kmer)){ 
                    read_local_hash_table_per_simplitig(bs_local); //changes l_of_curr_simplitig
                }
                if(per_simplitig_use_local_id == '1'){//using local table
                    int local_id = 0;
                    if(ceil(log2(l_of_curr_simplitig)) != 0){
                       local_id = bs_main.read_uint(ceil(log2(l_of_curr_simplitig)));
                    }
                     
                    int col_class = local_hash_table[local_id]; 
                    last_col_vector = global_table[col_class];
                    if(!TESTING_SPEED) dec_ess_color.fs << last_col_vector << endl;
                    written_kmer+=1;
                }else{
                    if(USE_HUFFMAN==false){
                        uint64_t col_class = bs_main.read_uint(lm);
                        last_col_vector = global_table[col_class];
                        if(!TESTING_SPEED) dec_ess_color.fs << last_col_vector << endl;
                        written_kmer+=1;
                    }else{
                        //TODO -- untested
                        uint32_t col_class = bs_main.HuffDecode(huff_root);
                        last_col_vector = global_table[col_class];
                        if(!TESTING_SPEED) dec_ess_color.fs << last_col_vector << endl;
                        written_kmer+=1;
                    }
                }
            }
            else if (c == '1')
            {
                char c2 = '1';
                if(per_simplitig_bigD != 0){
                    c2 = bs_main.read_one_bit();
                }
                
                if (c2 == '1')
                { // run
                    flush_skip_and_del(differ_run, last_col_vector,dec_ess_color);

                    {//paul method
                        int q = bs_main.read_number_encoded_in_unary_one();
                        char asst = bs_main.read_one_bit(); //should be  == '0');
                        int rem = bs_main.read_uint(lmaxrun);

                        int skip = q * max_run + rem;
                        while (skip)
                        {
                            if(!TESTING_SPEED) dec_ess_color.fs << last_col_vector << endl;
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
                        int differing_bit = bs_main.read_uint(lc);
                        differ_run.push_back(differing_bit);
                    }else if(per_simplitig_bigD == 2){
                        char c3 = bs_main.read_one_bit();;
                        if(c3 == '1'){ // 101 // read two
                            int differing_bit = bs_main.read_uint(lc);
                            differ_run.push_back(differing_bit);
                            differing_bit = bs_main.read_uint(lc);
                            differ_run.push_back(differing_bit);
                        }else{//c3 = 0, 100 // read one
                            int differing_bit = bs_main.read_uint(lc);
                            differ_run.push_back(differing_bit);
                        }
                    }
                }
            }
        } 
        flush_skip_and_del(differ_run, last_col_vector,dec_ess_color);
        time_end("Total run of decompression");
    }
};


int main (int argc, char* argv[]){
	cout<<"Version: "<<VERSION_NAME<<endl;

	vector<string> args(argv + 1, argv + argc);
    //string dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname; //string tmp_dir;
    string folder_name;
    int M, C;
    int max_run = 16;
	uint64_t num_kmers=0;
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: decompress -i <folder with mega.essd, meta.txt>" << endl;
            //meta gives kmer_size, colors, max_run
            //cout << "Syntax: tool -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run> -p <debug-mode>" << endl;
            return 0;
        } else if (*i == "-i") {
            folder_name = *++i;
        } 
        // if (*i == "-h" || *i == "--help") {
        //     cout << "Syntax: decompress -i <folder>" << endl;
        //     //meta gives kmer_size, colors, max_run
        //     cout << "Syntax: tool -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run> -p <debug-mode>" << endl;
        //     return 0;
        // } else if (*i == "-i") {
        //     dedup_bitmatrix_fname = *++i;
        // } else if (*i == "-d") {
        //     dup_bitmatrix_fname = *++i;
        // }else if (*i == "-c") {
        //     C = std::stoi(*++i);
        // }else if (*i == "-m") {
        //     M = std::stoi(*++i);
        // }else if (*i == "-k") {
        //     num_kmers = std::stol(*++i);
        // }else if (*i == "-s") {
        //     spss_boundary_fname = *++i;
		// }else if (*i == "-x") {
        //     max_run = std::stoi(*++i);
		// }else if (*i == "-p") {
        //     DEBUG_MODE = true;
		// }
		// // else if (*i == "-t") {
        // //     tmp_dir  = *++i;
		// // }
    }
    chdir(folder_name);
    //COLESS_Decompress cdec(num_kmers, M, C,  spss_boundary_fname, max_run);
    COLESS_Decompress cdec;
    
    //cdec.test_run();
    cdec.run();
    
	return EXIT_SUCCESS;
}