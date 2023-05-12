#include <sstream>
#include <iostream>
#include <string>
#include<fstream>
#include<vector>
#include <cstring>
#include <unistd.h>
#include<stdio.h>
#include <sys/stat.h>
//path2, name2 is the name of the file
//created, name1 is the string used in creating the symbolic  link

// int symlink(const char  *oldname, const char *symlinknewname);
#define PATH_MAX 4096

namespace CustomSysFns{
    const std::string WHITESPACE = " \n\r\t\f\v";
    
    std::string ltrim(const std::string &s)
    {
        size_t start = s.find_first_not_of(WHITESPACE);
        return (start == std::string::npos) ? "" : s.substr(start);
    }
    
    std::string rtrim(const std::string &s)
    {
        size_t end = s.find_last_not_of(WHITESPACE);
        return (end == std::string::npos) ? "" : s.substr(0, end + 1);
    }
    
    std::string trim(const std::string &s) {
        return rtrim(ltrim(s));
    }
    int basename_r(const char* path, char*  buffer, size_t  bufflen)
    {
        const char *endp, *startp;
        int         len, result;

        /* Empty or NULL string gets treated as "." */
        if (path == NULL || *path == '\0') {
            startp  = ".";
            len     = 1;
            goto Exit;
        }

        /* Strip trailing slashes */
        endp = path + strlen(path) - 1;
        while (endp > path && *endp == '/')
            endp--;

        /* All slashes becomes "/" */
        if (endp == path && *endp == '/') {
            startp = "/";
            len    = 1;
            goto Exit;
        }

        /* Find the start of the base */
        startp = endp;
        while (startp > path && *(startp - 1) != '/')
            startp--;

        len = endp - startp +1;

    Exit:
        result = len;
        if (buffer == NULL) {
            return result;
        }
        if (len > (int)bufflen-1) {
            len    = (int)bufflen-1;
            result = -1;
            errno  = ERANGE;
        }

        if (len >= 0) {
            memcpy( buffer, startp, len );
            buffer[len] = 0;
        }
        return result;
    }


    // string removeExtFromBase(string basefilename, string ext){
    //     std::string last_element(str.substr(str.rfind(":") + 2))

    // }

    const char *get_filename_ext(const char *filename) {
        const char *dot = strrchr(filename, '.');
        if(!dot || dot == filename) return "";
        return dot + 1;
    }

    

    char* BaseName(const char*  path) {
        static char* bname = NULL;
        int ret;

        if (bname == NULL) {
            bname = (char *)malloc(PATH_MAX);
            if (bname == NULL)
                return(NULL);
        }
        ret = basename_r(path, bname, PATH_MAX);
        return (ret < 0) ? NULL : bname;
    }

    int
    dirname_r(const char*  path, char*  buffer, size_t  bufflen)
    {
        const char *endp, *startp;
        int         result, len;

        /* Empty or NULL string gets treated as "." */
        if (path == NULL || *path == '\0') {
            startp = ".";
            len  = 1;
            goto Exit;
        }

        /* Strip trailing slashes */
        endp = path + strlen(path) - 1;
        while (endp > path && *endp == '/')
            endp--;

        /* Find the start of the dir */
        while (endp > path && *endp != '/')
            endp--;

        /* Either the dir is "/" or there are no slashes */
        if (endp == path) {
            startp = (*endp == '/') ? "/" : ".";
            len  = 1;
            goto Exit;
        }

        do {
            endp--;
        } while (endp > path && *endp == '/');

        startp = path;
        len = endp - startp +1;

    Exit:
        result = len;
        if (len+1 > PATH_MAX) {
            errno = ENAMETOOLONG;
            return -1;
        }
        if (buffer == NULL)
            return result;

        if (len > (int)bufflen-1) {
            len    = (int)bufflen-1;
            result = -1;
            errno  = ERANGE;
        }

        if (len >= 0) {
            memcpy( buffer, startp, len );
            buffer[len] = 0;
        }
        return result;
    }

    char* DirName(const char*  path) {
        static char*  bname = NULL;
        int           ret;

        if (bname == NULL) {
            bname = (char *)malloc(PATH_MAX);
            if (bname == NULL)
                return(NULL);
        }

        ret = dirname_r(path, bname, PATH_MAX);
        return (ret < 0) ? NULL : bname;
    }
    std::string getFullExt(const char * filename, std::string& stripped_tokenizer){
        stripped_tokenizer = "";
        std::string fullext="";
        filename = BaseName(filename);
        std::vector <std::string> tokens;
        // stringstream class check1
        std::stringstream check1(filename);
        
        std::string intermediate;
        
        // Tokenizing w.r.t. space ' '
        while(getline(check1, intermediate, '.'))
        {
            tokens.push_back(intermediate);
        }
    
        for(int i = 0; i < tokens.size()-1; i++){
            stripped_tokenizer = ""+tokens[tokens.size()-i-1] + stripped_tokenizer;
            fullext = "."+tokens[tokens.size()-i-1] + fullext;
        }
        
        const char* allowed_ext[] = {"fastq.gz", "fastq", "fq","fq.gz", "fa","fa.gz","fna","fna.gz", "fasta","fasta.gz","unitigs.fa","unitigs.fa.gz","essd","essd.gz"};
        for(auto s : allowed_ext){
            std::string st(s);
            if(fullext == "."+st){
                std::cout<<"Checked that input file extension is correct... success: "<< fullext<<std::endl;
                return fullext;
            }
        }
        std::cout<< "file extension not supported: "<<fullext<<std::endl;
        exit(5);
    }

// int main(int argc, char **argv) {
//     if (argc != 2)
//         printf("usage: bin [path/to/name]\n");
//     else {
//         printf("basename=%s\n", BaseName(argv[1]));
//         char cmd[PATH_MAX];
//         sprintf(cmd, "basename %s", argv[1]);
//         system(cmd);

//         printf("dirname=%s\n", DirName(argv[1]));
//         sprintf(cmd, "dirname %s", argv[1]);
//         system(cmd);
//         printf("original:%s\n", argv[1]);
//     }
//     return 0;
// }
};
using namespace CustomSysFns;

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
    string str = "[\"";
    while(getline(ifs, line)){
        if(lastline!=""){
            str+=lastline+"\",\"";
        }
        lastline=trim(line);
    }
    str+=lastline;
    str+="\"]";
    ifs.close();
    return str;
}
string enclose_listelem_with_brackets(vector<string> list){
    string line="";
    string lastline="";
    string str = "[\"";
    int i = 0;
    while(i!=list.size()){
        line = list[i];
        if(lastline!=""){
            str+=lastline+"\",\"";
        }
        lastline=trim(line);
        i++;
    }
    str+=lastline;
    str+="\"]";
    return str;
}
void write_config(string field, string value){
    config<<field<<": "<<value<<endl;
}
void write_config_str(string field, string value){
    config<<field<<": "<<"\""<<value<<"\""<<endl;
}

vector<string> read_list(string filename){
    vector<string> filenames;
    ifstream is(filename);
    string line;
    while(getline(is, line)){
        filenames.push_back(line);
    }
    return filenames;
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
// -j [int]          Default=1. Number of threads.   
// -v                Enable verbose mode: print more useful information.
// -h                Print this Help
// -V                Print version number

std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer, sizeof buffer, pipe) != NULL) {
        result += buffer;
    }
    pclose(pipe);
    return result;
}

void check_if_program_exist(const char* file){
    string ss = "which "+string(file)+" > /dev/null 2>&1";     //test using const char* file = "nano";
    string ss2 = "which "+string(file);
    if(system(ss.c_str())){  //program doesnt exist
        cout<<"Error: could not find progam '"<<file<<"' in path."<<endl;
        exit(2);
    }else{
        cout<<"Found '"<<file<<"' in path: "<<exec(ss2.c_str());
    }
}
int main(int argc, char **argv) {
    #if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
    #else
        cout<<"This program does not support this OS."<<endl;
        exit(2);
    #endif

    cout<<"Checking if all dependencies are installed... "<<endl;
    check_if_program_exist("sort");
    // check_if_program_exist("ggcat");
    // check_if_program_exist("essCompress");
    // check_if_program_exist("kmc");
    // check_if_program_exist("snakemake");

    vector<string> args(argv + 1, argv + argc);
   
    uint64_t NUM_STRING=0;
    uint64_t NUM_KMER=0;
    
    //mandatory input
    int k=-1;
    string list_file = "";

    //default value provided
    int ab=1;
    string OUT_DIR = ".";

    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: ./essColor -k <kmer-size> -a <min-abundance> -j <num-threads> -i <Path to input file: list of samples> -o <output directory>" << endl;
            return 0;
        } else if (*i == "-k") {
            k = std::stoi(*++i);
        } else if (*i == "-i") {
            list_file = *++i;
        } else if (*i == "-o") {
            OUT_DIR = *++i;
        } 
    }

    if(k==-1 || ab==-1 || list_file==""){
        cout<<"Error: mandatory input args not provided"<<endl;
        exit(4);
    }

    if(ab<1){
        cout<<"Error: min abundance should be at least 1."<<endl;
        exit(3);
    }
    // if (chdir(OUT_DIR.c_str()) != 0){
    //     cout<<("chdir() failed");
    //     exit(3);
    // }

    //cout<<getFullExt("/Users/aur1111/projects/git-medvedevgroup/color-compression/a.b.fastq.gz")<<endl;

  
    config.open("config.yaml");
    write_config("k", to_string(k));
    write_config("ab", to_string(ab));
    write_config("mem", to_string(100));
    write_config_str("unitig", "ggcat");
    write_config_str("ess", "notip");
    write_config_str("option", "");
    write_config_str("matrix_generator", "genmatrix");




    // create symlink
    vector<string> inputnames = read_list(list_file);

    string prev_ext="";
    string ext="";
    string stripped_ext="";
    for (int i = 0; i< inputnames.size(); i++){
        string inputname = inputnames[i];
        string basename(BaseName(inputname.c_str()));
        ext=getFullExt(inputname.c_str(),stripped_ext);
        if(i!=0){
            if(prev_ext!=ext){
                cout<<"All input files must have the same extension."<<endl;
                exit(3);
            }
        }else{
             system(("rm -rf "+stripped_ext).c_str());
            system(("mkdir -p "+stripped_ext).c_str());
            system("mkdir -p esscolor");

        }
        //cout<<inputname.c_str()<<" "<<(stripped_ext+"/"+basename).c_str()<<endl;
        system (("ln -s "+inputname+" "+stripped_ext+"/"+basename).c_str());
       
        //cout<<ext<<endl;
        prev_ext = ext;
        basename.replace(basename.find(ext), sizeof(ext) - 1, "");
        inputnames[i] = basename;
        
        
    }
    write_config("SAMPLES", enclose_listelem_with_brackets(inputnames)); 
    write_config_str("EXTENSION", ext);  // EXTENSION: ".fa"
    if(ext==".fastq" || ext==".fastq.gz" || ext==".fq"  || ext==".fq.gz" ){
        write_config_str("src", "q");
    }else if(ext==".fna" || ext==".fna.gz" || ext==".fa"  || ext==".fa.gz" || ext==".fasta" || ext==".fasta.gz" ){
        write_config_str("src", "m");
    }

    //system("snakemake -s ");
    return 0;
}


// fa, fq, fastq, fastq.gz 


