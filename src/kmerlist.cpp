//g++ kmerlist.cpp -o kmerlist
#include <iostream>
#include <string>
#include <stdint.h>
#include <vector>
#include <fstream>
#include<algorithm>
#include<cmath>

#include <string>
#include <iostream>                                                                                                                                                                                          
#include <algorithm>
#include <cassert>

using namespace std;

inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

char complement(char n)
{   
    switch(n)
    {   
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }   
    assert(false);
    return ' ';
}   

string retsmall(string nucs){
  string revs=nucs;
  bool rev=false;
  size_t n = nucs.length();
  for (size_t i = 0; i < n; i++){
    char comp = complement(nucs[n-i-1]);
    if(!rev){
      if(comp  < nucs[i]){
          rev=true;
      }
      if(nucs[i] < comp){
        break;
      }
    }
    revs[i] = comp;
  }
  if(rev){
    return revs;
  }
  return nucs;
}

void revcomp(string& nucs)
{ 
  transform(
        begin(nucs),
        end(nucs),
        begin(nucs),
        complement);
  std::reverse(nucs.begin(), nucs.end());
}
 
string canon(string str){
  string fwd = str;
  revcomp(str);
  if(str<fwd ){
    return str;
  }
  return fwd ;
}


class OutFile{
  public:
  string filename;
  std::ofstream fs;

  OutFile(string filename){
    this->filename = filename;
    this->fs.open (filename,  std::fstream::out );
  }
  void writeLine(string towrite){
      this->fs << towrite <<endl;
  }
  ~OutFile(){
    this->fs.close();
  }
};


class InFile{
  public:
  string filename;
  std::ifstream fs;

  InFile(string filename){
    this->filename = filename;
    this->fs.open (filename,  std::fstream::in );
  }
  ~InFile(){
    this->fs.close();
  }
};

int main(int argc, char **argv) {
    vector<string> args(argv + 1, argv + argc);
    string in_fname;
    uint64_t NUM_STRING=0;
    uint64_t NUM_KMER=0;
    int k;

    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: ./kmerlist -k <kmer-size> -u <path to union spss>" << endl;
            return 0;
        } else if (*i == "-k") {
            k = std::stoi(*++i);
        } else if (*i == "-u") {
            in_fname = *++i;
        } 
    }

    OutFile out_bv("ess_boundary.txt");
    OutFile out_bv1("ess_boundary_bit.txt");
    OutFile out_kmer_order("ess_kmer_id.txt");
    InFile in_f(in_fname);

    string seq;

    uint64_t i=0;
    while (getline(in_f.fs, seq)){
      seq=trim(seq);
      if(seq[0]!='>'){
        out_bv.fs<<i<<endl;
        for(int j=0; j<seq.length()-k+1; j++){
              string kmer = seq.substr(j, k);
              if(j!=0){
                out_bv1.fs<<"0"<<endl;
              }else{
                out_bv1.fs<<"1"<<endl;
              }
              out_kmer_order.fs<<retsmall(kmer)<<" "<<i<<endl;
              i=i+1;
        }
        NUM_STRING += 1;
        NUM_KMER += seq.length() - k + 1;    
      }      
    }
    cout<<"NUM_STRING="<<NUM_STRING<<endl;
    cout<<"NUM_KMER="<<NUM_KMER<<endl;
}



