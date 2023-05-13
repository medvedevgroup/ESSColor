//
//  stat.hpp
//  UST
//
//  Created by Amatur Rahman on 23/5/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//

#ifndef stat_h
#define stat_h

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

class Stat{
public:
    typedef int32_t NODE_T;
    typedef long long unsigned int TOTALCHAR_T;
    
    TOTALCHAR_T nKmers = 0;
    TOTALCHAR_T C_bcalm = 0;
    TOTALCHAR_T E_bcalm = 0;
    NODE_T V_bcalm = 0;
    
    template <class T>
    void statPrinter(ofstream& FOUT, string label, T a, bool consoleOutput = false){
        FOUT<<label<<"="<<a<<endl;
        if(consoleOutput){
            cout<<label<<"="<<a<<endl;
        }
    }
};

class ProfileStat: public Stat{
public:
    TOTALCHAR_T charLowerbound = 0;
    NODE_T isolated_node_count = 0;
    NODE_T sink_count = 0;
    NODE_T source_count = 0;
    NODE_T sharedparent_count = 0;
    NODE_T sharparentCntRefined = 0;
    NODE_T onecount = 0;
    
    NODE_T walkstarting_node_count = 0;
    float upperbound = 0;
    
};

class USTStat: public Stat{
    public:
    NODE_T V_ust = 0;
    TOTALCHAR_T C_ust = 0;
};

class ESSTipStat: public Stat{
    public:
    NODE_T V_esstip = 0;
    TOTALCHAR_T C_esstip = 0;

    NODE_T V_nontip = 0;
    TOTALCHAR_T C_nondna_esstip = 0;
};


class ESSStat: public Stat{
    public:
    NODE_T V_ess = 0;
    TOTALCHAR_T C_ess = 0;
    
    NODE_T V_ust = 0;
    
    NODE_T V_non_absorbed = 0;
    TOTALCHAR_T C_nondna_ess = 0;
    
    NODE_T absorbGraphNumCC_endAbosrb = 0;
};

vector<string> split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

//void substat(string globalStatFileName){
//    ifstream gl(globalStatFileName);
//    string fileline;
//    while (getline (gl, fileline)) {
//        vector<string> v = split (fileline, '=');
//        v[0]
//        for (auto i : v) cout << i << endl;
//        result.push_back (item);
//    }
//
//    getline(unitigFile, line);
//
//
//}



#endif /* stat_h */
