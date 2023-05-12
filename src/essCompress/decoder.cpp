//
//  decoder.h
//
//  Created by Amatur Rahman on 28/11/19.
//  Copyright © 2019 Penn State. All rights reserved.
//

#include <ctype.h>
#include<string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//

#include "decoder.hpp"
using namespace std;

string remove_extension(const string& filename) {
    string out;
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) {
        out = filename;
    }else{
        out = filename.substr(0, lastdot);
    }
    if(out.length()==0){
        out = "out";
    }
    return out;
}

string basename(const string& s) {

   char sep = '/';

#ifdef _WIN32
   sep = '\\';
#endif

   size_t i = s.rfind(sep, s.length());
   if (i != string::npos) {
      return(s.substr(i+1, s.length() - i));
   }

   return("");
}

bool isFileExist(string& fileName) {
    FILE *fp = fopen(fileName.c_str(), "r");
    if (fp) {
        fclose(fp);
        return true;
    }
    return errno != ENOENT;
}


vector<string> split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

void parseHeader(string pathname, int &K, int& MODE){
    ifstream infile(pathname);
    string sLine;
    getline(infile, sLine);
    
    vector<string> v = split(sLine, '_');
    K = std::stoi(v[1]);
    MODE = std::stoi(v[2]);
    infile.close();
}

int main(int argc, char** argv) {
    DDEBUG = 0;
    bool HEADER_FIX=false;
    vector<string> args(argv + 1, argv + argc);
    std::string pathname;
    // = argv[1];

        for (auto i = args.begin(); i != args.end(); ++i) {
            if (*i == "-h" || *i == "--help") {
                cout << "Syntax: ./essAuxDecompress -i <ess-compressed-file> [-f]. -f flag indicates header should be a full line " << endl;
                return 0;
            } else if (*i == "-f") {
                HEADER_FIX = true;
            } else if (*i == "-i") {
                pathname = *++i;
            } 
        }
    cout<<pathname<<endl;
    cout<<HEADER_FIX<<endl;

    //if(DDEBUG) cout<<"----------------RUNNING IN DEBUG MODE--------------------"<<endl;
    /*
    const char* nvalue = "" ;
    int K=0;
    string filename = "ust_ess_abs.txt";
    bool tip_mode = 0;

        int c ;

            while( ( c = getopt (argc, argv, "i:k:t:") ) != -1 )
            {
                switch(c)
                {
                    case 'i':
                        if(optarg) nvalue = optarg;
                        break;
                    case 't':
                        if(optarg) {
                            tip_mode = static_cast<bool>(std::atoi(optarg));
                        }
                        break;
                    case 'k':
                        if(optarg) {
                            K = std::atoi(optarg) ;
                            if(K<=0){
                                fprintf(stderr, "Error: Specify a positive k value.\n");
                                exit(EXIT_FAILURE);
                            }
                        }else{
                            fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                                    argv[0]);
                            exit(EXIT_FAILURE);
                        }
                        break;
                    default: //
                        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-path>\n",
                                argv[0]);
                        exit(EXIT_FAILURE);

                }
            }

    
    

    if(K==0 || strcmp(nvalue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    // if(tip_mode){
    //     filename = "ust_ess_tip.txt";
    // }
     
     string pathname = string(nvalue);
     */
    
    
    
    if(!isFileExist(pathname)){
        cout<<"File \""<<pathname<<"\""<<"does not exist."<<endl;
        exit(EXIT_FAILURE);
    }

    //pathname+="/"+filename;
    //int K=31, string ENCODED_FILE= "tipOutput.txt"

    int K;
    int MODE;
    //string outputFilename = remove_extension(basename(pathname));
    string outputFilename = "kmers.";
    parseHeader(pathname, K, MODE);
    bool tip_mode = MODE;
    
    if(tip_mode){
        decodeTip(K, pathname, outputFilename+"esstip.spss", HEADER_FIX);
        //cout<<"ESS-Tip-Compress (core) decoding done!"<<endl;
        //cout<<"Output SPSS is in file \""<<outputFilename+"esstip.spss"<<"\""<<endl;
    }else{
        decodeOneAbsorb(K, pathname, outputFilename+"ess.spss", HEADER_FIX);
        //cout<<"ESS-Compress  (core)  decoding done!"<<endl;
        //cout<<"Output SPSS is in file \""<<outputFilename+"ess.spss"<<"\""<<endl;
    }


    // cout<<GetStdoutFromCommand("top -p "+ std::to_string(pid)+" -b -n1 | tail -1 |  awk '{ print $5}'")<<endl;
    // cout<<"printed memory requirement"<<endl;
    //

    return 0;
}

