//
//  param.hpp
//  UST
//
//  Created by Amatur Rahman on 15/6/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//

#ifndef param_h
#define param_h

#include "misc.hpp"
#include<string>

using namespace std;
class Param{
public:
    bool VERBOSE_MODE = false;
    bool VALIDATE = false;
    bool PROFILE_AND_STAT = false;
//    string APP_PATH_DSK="bin/dsk";
//    string APP_PATH_DSK2ASCII="bin/dsk2ascii";
//    string APP_PATH_BCALM="bin/bcalm";
    
   // string DSK_PATH="/Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/";
    //string DSK_PATH="/home/aur1111/w/dsk/build/bin";
    ///home/aur1111/w/dsk/build/bin

    string DSK_PATH = "";
    string UNITIG_FILE;
    string OUTPUT_FILENAME;
    
    string VERSION = "2.0";
    
    
    inline string getBcalmFileBasename(){
        return "kmers";
        string a = getFileName(UNITIG_FILE).substr(0, getFileName(UNITIG_FILE).length()-11);
        if(a.empty()){
            a = "file";
        }
        return a;
    }
};

class ESSTipParam:public Param{
public:
    
};

class ESSParam : public Param{
public:
    bool GETLOWERBOUND_CC = false;
    bool PROFILE_ONLY_ABS = false;
    bool EARLYABSORB = true;
    //bool VALIDATE = true;
    
    /*FILENAMES*/
    string ofileTipOutput = "kmers.ess";
};
#endif /* param_h */
