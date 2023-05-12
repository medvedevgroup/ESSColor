// VERSION 2.0

//
//  decoder.hpp
//
//  Created by Amatur Rahman on 28/11/19.
//  Copyright © 2019 psu. All rights reserved.
//
//

#ifndef decoder_h
#define decoder_h

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdint.h>
#include <unistd.h>

#include "misc.hpp"


using namespace std;

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

#define NONDNA_START_TIP_SUF '('
#define NONDNA_END_TIP_SUF ')'
#define NONDNA_START_TIP_PRE '{'
#define NONDNA_END_TIP_PRE '}'
#define NONDNA_PLUS "+"
#define NONDNA_MINUS "-"
#define NONDNA_START "["
#define NONDNA_END "]"

//
//string reverseComplement(string base)
//{
//    size_t len = base.length();
//    char *out = new char[len + 1];
//    out[len] = '\0';
//    for (int i = 0; i < len; i++)
//    {
//        if (base[i] == 'A')
//            out[len - i - 1] = 'T';
//        else if (base[i] == 'C')
//            out[len - i - 1] = 'G';
//        else if (base[i] == 'G')
//            out[len - i - 1] = 'C';
//        else if (base[i] == 'T')
//            out[len - i - 1] = 'A';
//    }
//    string outString(out);
//    free(out);
//    return outString;
//}

int decodeTip(int K, string UNITIG_FILE="ust_ess_tip.txt", string OUT_FILE = "absorbDecompressed.fa" ){
    ifstream unitigFile;
    unitigFile.open(UNITIG_FILE);

    ofstream outFile;
    outFile.open(OUT_FILE);

    string line;
    bool startPrefCut = false;
    bool startSufCut = false;
    string pref;
    string suf;
    int walkid = 0;

    while (getline(unitigFile, line)) {
        //cout<<line<<endl;
        if (line.empty() || line.substr(0, 1).compare(">") == 0) {

        } else {
            walkid++;

            startPrefCut = false;
            startSufCut = false;

            string tip = "";
            string sbuf = "";
            string lastk = "";
            string pref = "";
            string suf = "";

            //cout<<walkid-1<<" "<<line[0]<<endl;
            for(int i = 0; i<line.length(); i++){
                if(line[i]=='A'|| line[i]=='C' || line[i]=='G' || line[i]=='T'){
                    if(startPrefCut){
                        tip += line[i];
                    }else if(startSufCut){
                        tip += line[i];
                    }else{
                        sbuf +=   line[i];
                        if(lastk.length() < K - 1){
                            lastk += line[i];
                        }else{
                            lastk = lastk.substr(1, K-2) + line[i];
                        }
                    }
                }else if(line[i]=='(' || line[i]==')'){
                    if(startPrefCut){ //already has one
                        startPrefCut = false;
                        outFile<< ">" ;
                        if(DDEBUG){
                            outFile<< walkid << " pref"; //for Debug
                        }
                        outFile<< endl;
                        outFile<< (pref + tip) << endl;
                        tip = "";
                    }else if(!startPrefCut){ //prefix is cut starts
                        startPrefCut = true;
                        pref = lastk;
                    }
                }else if(line[i]=='{' || line[i]=='}'){ //suffix is cut
                    if(startSufCut){ //already has one
                        startSufCut = false;
                        
                        outFile<< ">" ;
                        if(DDEBUG){
                            outFile<< walkid << " suf"; //for Debug
                        }
                        outFile<< endl;
                        outFile<< (tip + suf) << endl;
                        tip = "";
                    }else if(!startSufCut){
                        startSufCut = true;
                        suf = lastk;
                    }
                }
            }
            //print the tips before printing the contig
            //outFile<<">"<<walkid-1<<" : "<<"\n"; for Debug
            
            outFile<< ">" ;
            if(DDEBUG){
                outFile<< walkid; //for Debug
            }
            outFile<<"\n";
            outFile<<sbuf << endl;
        }
    }

    unitigFile.close();
    outFile.close();

    //cout << "Complete conversion." << endl;
    return EXIT_SUCCESS;
}


void updateCurrOverlap(char c, string& lastk1, int K){
    if(lastk1.length() < K - 1){
        lastk1 += c;
    }else{
        lastk1 = lastk1.substr(1, K-2) + c;
    }
}

void decompressEnclosed(string& s, int& i, string overlapFromParent, int K, ofstream & FOUT, int stringID){
    //assert(overlapFromParent.length()==K-1);
    ASSERT(overlapFromParent.length()==K-1, "The +/- should be replace by a string of length " << K-1 << ", but we receive string of length " << overlapFromParent.length() << "\n");

    string sOut = "";
    string currOverlap;

    stack<tuple<string, string, string> > decstack;

    while(i < s.length()){
        char c = s[i++];
        if(c=='['){
            decstack.push(make_tuple(currOverlap, overlapFromParent, sOut));
            overlapFromParent=currOverlap;
            sOut="";
        }else if(c==']'){
            if(DDEBUG){
            FOUT<<">"<<stringID<<"\n"<<sOut<<endl;
                
            }else{
            FOUT<<">\n"<<sOut<<endl;

            }
            if(!decstack.empty()){
                currOverlap =  get<0>(decstack.top());
                overlapFromParent = get<1>(decstack.top());
                sOut = get<2>(decstack.top());
                decstack.pop();
            }else{
                return;
            }
        }else if(c=='+'){
            sOut += overlapFromParent;
            //update currOverlap
            currOverlap = overlapFromParent;
            //at this point overlapFromParent is useless
            
        }else if(c=='-'){
            sOut += reverseComplement(overlapFromParent);
            //update currOverlap
            currOverlap = reverseComplement(overlapFromParent);
            //at this point overlapFromParent is useless
        }else{
            sOut += c;
            updateCurrOverlap(c,currOverlap, K);
            //update currOverlap
        }
    }
}

//vector<string> decompressEnclosed(string& s, int& i, string overlapFromParent, int K, ofstream & FOUT){
//    vector<string> S;
//    string sOut;
//    string currOverlap;
//
//    while(i < s.length()){
//        char c = s[i++];
//        if(c=='['){
//            vector<string> S1 = decompressEnclosed(s, i,currOverlap, K, FOUT);
//            for(int j = 0; j < S1.size(); j++){
//                //S.push_back(S1.at(j));
//            }
//        }else if(c==']'){
//            //S.push_back(sOut);
//            FOUT<<">\n"<<sOut<<endl;
//            return S;
//        }else if(c=='+'){
//            sOut += overlapFromParent;
//            //update currOverlap
//            currOverlap = overlapFromParent;
//        }else if(c=='-'){
//            sOut += reverseComplement(overlapFromParent);
//            //update currOverlap
//            currOverlap = reverseComplement(overlapFromParent);
//        }else{
//            sOut += c;
//            updateCurrOverlap(c,currOverlap, K);
//            //update currOverlap
//        }
//
//    }
//    return S;
//}

void decompress(string& s, int K, ofstream & FOUT, int stringID){
    vector<string> S;
    string sOut;
    string currOverlap;
    int i = 0;
    while(i < s.length()){
        char c = s[i++];
           if(c=='['){
               decompressEnclosed(s, i,currOverlap, K, FOUT, stringID);
           }else{
               sOut += c;
               updateCurrOverlap(c, currOverlap, K);  //update currOverlap
           }
    }
   // cout<< "Complete conversion of: " <<s << endl;
    
    if(DDEBUG){
        FOUT<<">"<<stringID<<"\n"<<sOut<<endl;  
    }else{
        FOUT<<">\n"<<sOut<<endl;
    }
    
}
void decodeOneAbsorb(int K, string ENCODED_FILE= "ust_ess_abs.txt", string OUT_FA_FILE="absorbDecompressed.fa"){
    ifstream encodedFile;
    encodedFile.open(ENCODED_FILE);

    ofstream outFile;
    outFile.open(OUT_FA_FILE);
    string line;
    int nWalks = 0;
    while (getline(encodedFile, line)) {
        //cout<<line<<endl;
        if (line.empty() || line.substr(0, 1).compare(">") == 0) {

        } else {
            nWalks++;
            decompress(line, K, outFile, nWalks);
        }
    }
    encodedFile.close();
    outFile.close();
}



#endif /* decoder_h */
