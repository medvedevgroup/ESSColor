// Last Edited: July 19, 2020

#include "ust.hpp"
#include "tip.hpp"
#include "ess.hpp"
#include "profileGraph.hpp"

using namespace std;

int main(int argc, char** argv) {
#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#else
    cout<<"ESS-Compress does not support this OS."<<endl;
    return EXIT_FAILURE; //exit(2);
#endif

    if(DDEBUG) cout<<"----------------RUNNING IN DEBUG MODE--------------------"<<endl;
    bool PROFILE_AND_STAT = false;
    bool VERBOSE_MODE = false;

    int TYPE_0_1_2 = 0;
    int K = 0;
    
    char * nvalue = NULL;
    string graph_file_name = "";
    
    char c;
    while( ( c = getopt (argc, argv, "i:k:t:p:v:") ) != -1 )
    {
        switch(c)
        {
            case 't':
                if(optarg) {
                    if(strcmp(optarg, "0")==0 || strcmp(optarg, "1")==0 ||  strcmp(optarg, "2")==0){
                        TYPE_0_1_2 = static_cast<int>(std::atoi(optarg));
                    }else{
                        fprintf(stderr, "Error: use either -t 0 or -t 1 or -t 2 \n");
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'p':
                if(optarg) {
                    if(strcmp(optarg, "0")==0 || strcmp(optarg, "1")==0){
                        PROFILE_AND_STAT = static_cast<bool>(std::atoi(optarg));
                    }else{
                        fprintf(stderr, "Error: use either -p 0 or -p 1 \n");
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'v':
                if(optarg) {
                    if(strcmp(optarg, "0")==0 || strcmp(optarg, "1")==0){
                        VERBOSE_MODE = static_cast<bool>(std::atoi(optarg));
                    }else{
                        fprintf(stderr, "Error: use either -v 0 or -v 1 \n");
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'i':
                if(optarg) nvalue = optarg;
                break;
            case 'k':
                if(optarg) {
                    K = std::atoi(optarg) ;
                    if(K<=0){
                        fprintf(stderr, "Error: Specify a positive k value.\n");
                        exit(EXIT_FAILURE);
                    }
                }else{
                    fprintf(stderr, "Usage: %s -k <kmer-size> -i <input-file-name> [-v <verbosity level; 0 (DEFAULT): non-verbose, 1: verbose>] [-t <0 (DEFAULT): max compression, 1: fast mode (low compression), 2: output spss only with no compression>]\n",
                            argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                fprintf(stderr, "Usage: %s -k <kmer-size> -i <input-file-name> [-v <verbosity level; 0 (DEFAULT): non-verbose, 1: verbose>] [-t <0 (DEFAULT): max compression, 1: fast mode (low compression), 2: output spss only with no compression>]\n",
                            argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    if(K==0 || strcmp(nvalue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-f <0: for default mode (maximum compression), 1: for fast mode (low compression), default = 0>]\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }
    graph_file_name = string(nvalue);
    
    
    if(PROFILE_AND_STAT){
        ProfileGraph().run(graph_file_name, K);
    }else if(TYPE_0_1_2 == 2){
        UST ust;
        ust.run(graph_file_name, K, VERBOSE_MODE);
    }else if(TYPE_0_1_2 == 1){
        ESSTip esstip;
        esstip.run(graph_file_name, K, VERBOSE_MODE);
    }else{
        AbsorbGraph ess;
        ess.run(graph_file_name, K, VERBOSE_MODE);
    }
}
