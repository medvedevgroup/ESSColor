//  Last modified by Amatur Rahman on 19/7/20.

#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <list>
#include <stack>
#include <unordered_map>
#include <utility>
#include <queue>
#include <deque>
#include <tuple>
#include <cstring>
#include <string>

#ifndef misc_h
#define misc_h

using namespace std;


#ifdef DEBUGMODE
	int DDEBUG=1;
#else
	int DDEBUG=0;
#endif
//int DDEBUG = 1;


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



double readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}

string getFileName(const string& s) {
    char sep = '/';
    
    size_t i = s.rfind(sep, s.length());
    if (i != string::npos) {
        return(s.substr(i+1, s.length() - i));
    }
    
    return("");
}

// string getParentPath(string str){   
//     char * lastSlash = strrchr( str, '/');
//     if ( *lastSlash != '\n') *(lastSlash +1) = '\n';
//     string s(lastSlash);
//     return s;
// }

string delSpaces(string &str) {
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

/**
int equal_files(){
    int result = system("diff /Users/Sherlock/Documents/bcl/bcl/a.txt /Users/Sherlock/Documents/bcl/bcl/output-my.txt >nul 2>nul");
    //They are different if result != 0
    return result;
}
**/


inline string currentDateTime() {
    // Get current date/time, format is YYYY-MM-DD HH:mm:ss
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof (buf), "%Y-%m-%d %X\n", &tstruct);
    return buf;
}

string reverseComplement(string base) {
    size_t len = base.length();
    char* out = new char[len + 1];
    out[len] = '\0';
    for (int i = 0; i < len; i++) {
        if (base[i] == 'A') out[len - i - 1] = 'T';
        else if (base[i] == 'C') out[len - i - 1] = 'G';
        else if (base[i] == 'G') out[len - i - 1] = 'C';
        else if (base[i] == 'T') out[len - i - 1] = 'A';
    }
    string outString(out);
    free(out);
    return outString;
}

inline string plus_strings(const string& a, const string& b, size_t kmersize) {
    if (a == "") return b;
    if (b == "") return a;
    string ret = a + b.substr(kmersize - 1, b.length() - (kmersize - 1));
    return ret;
}

#endif /* misc_h */
