#include <sstream>
#include <iostream>
#include <string>
#include<fstream>
#include<vector>

using namespace std;
// k: 5
// ab: 1
// mem: 100
// type: "m"
// option: ""
// matrix_generator: "genmatrix"

// SAMPLES: ["c1","c2","c3"]
// EXTENSION: ".fa"

// const char config[] = "url=http://example.com\n"
//                       "file=main.exe\n"
//                       "true=0";

// std::istringstream is_file(config);

// std::string line;
// while( std::getline(is_file, line) )
// {
//   std::istringstream is_line(line);
//   std::string key;
//   if( std::getline(is_line, key, '=') )
//   {
//     std::string value;
//     if( std::getline(is_line, value) ) 
//       store_line(key, value);
//   }
// }
string extension_remove(string pathname){
    //string filename = "C:\\MyDirectory\\MyFile.bat";

    // Remove directory if present.
    // Do this before extension removal incase directory has a period character.
    string filename = pathname;
    const size_t last_slash_idx = filename.find_last_of("/");
    if (std::string::npos != last_slash_idx)
    {
        filename.erase(0, last_slash_idx + 1);
    }

    // Remove extension if present.
    const size_t period_idx = filename.rfind('.');
    if (std::string::npos != period_idx)
    {
        filename.erase(period_idx);
    }
    return filename;
}

ofstream config;
string enclose_listelem_with_brackets(string filewithlist){
    ifstream ifs(filewithlist);
    string line="";
    string lastline="";
    string str = "[";
    while(getline(ifs, line)){
        if(lastline!=""){
            str+=lastline+",";
        }
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
        lastline=line;
    }
    str+=lastline;
    str+="]";
    ifs.close();
    return str;
}
void write_config(string field, string value){
    config<<field<<": "<<value<<endl;
}

#include <iostream>

// mandatory arguments:
// -k [int]          k-mer size (must be >=4)
// -i [input-file]   Path to input file. Input file is a single text file containing the list of multiple fasta/fastq files (one file per line)
// -o [output-dir]   Path to output directory.

// optional arguments:
// -r [int]          Default=16. runLength: parameter runLength as described in paper
// -a [int]          Default=1. Sets a threshold X, such that k-mers that appear less than X times in the input dataset are filtered out. 
// -o [output-dir]   Specify output directory
// -t [int]          Default=1. Number of threads.   
// -v                Enable verbose mode: print more useful information.
// -h                Print this Help
// -V                Print version number
int main(int argc, char **argv) {
    vector<string> args(argv + 1, argv + argc);
    string list_file;
    uint64_t NUM_STRING=0;
    uint64_t NUM_KMER=0;
    int k;

    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: ./config_writer -k <kmer-size> -i <Path to input file: list of samples>" << endl;
            return 0;
        } else if (*i == "-k") {
            k = std::stoi(*++i);
        } else if (*i == "-i") {
            list_file = *++i;
        } 
    }
    config.open("config.yaml");
    write_config("k", "31");
    write_config("SAMPLES", enclose_listelem_with_brackets(list_file));
    // EXTENSION: ".fa"

    return 0;
}
