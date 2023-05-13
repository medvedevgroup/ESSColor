//
//  ust.h
//  UST
//
//  Created by Amatur Rahman on 21/5/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//
#ifndef ust_h
#define ust_h

#include<tuple>
#include<vector>
#include<string>
#include<map>
#include<unordered_map>

#include "param.hpp"
#include "stat.hpp"
#include "spss.hpp"

using namespace std;

class UST : public SPSS{
private: //different var/functions in UST and ESS
    typedef struct {
        bool isWalkEnd = false;
        unitig_t pos_in_walk = -100;
        unitig_t finalWalkId = -1; // renders some walkId as invalid
    } new_node_info_t;
    UST::new_node_info_t* oldToNew;
    
    unitig_t &V_bcalm = stat.V_bcalm;
    
    typedef tuple<unitig_t,unitig_t,unitig_t> threetuple; // uid, pathid, pos
    vector<threetuple> sorter;
    bool static sort_by_walkId (const threetuple &lhs, const threetuple &rhs){
        return get<1>(lhs) < get<1>(rhs);
    }
    bool static sort_by_pos (const threetuple &lhs, const threetuple &rhs){
        return get<2>(lhs) < get<2>(rhs);
    }
    
    USTStat stat;
    Param param;    /*!< holds all necessary parameters*/
    unitig_t countNewNode = 0; // number of paths in initial phase

    void collectInput(int argc, char** argv, string & graph_file_name, int & K, bool & FLG_ABUNDANCE);
    void readUnitigFile(const string& unitigFileName, vector<unitig_struct_t>& unitigs, vector<vector<edge_t> >& adjList);
    vector<threetuple> sorterMaker();
    void ustOutputToDisk(vector<threetuple>& sorter);

public:
    void run(string graph_file_name, int K, bool runFlag);
    ~UST();
};


void UST::run(string graph_file_name, int K, bool runFlag){
    param.VERBOSE_MODE = runFlag;
    cout<<"Running UST"<<endl;
    
    
    this->K= K;
    param.UNITIG_FILE = graph_file_name;
    //collectInput(argc, argv, graph_file_name, K, FLG_ABUNDANCE); //input: argc, argv
    
    if(param.VERBOSE_MODE) cout << "## [1] Please wait while input file (" << graph_file_name << ") is being read : k = "<<K<<endl;
    double startTime = readTimer();
    
    this->FLG_ABUNDANCE = false;
    this->readUnitigFile(param.UNITIG_FILE, unitigs, adjList); //input: graph_filename, output: last 2
    
    double TIME_READ_SEC = readTimer() - startTime;
    if(param.VERBOSE_MODE) cout<<"Done. TIME to read file "<<TIME_READ_SEC<<" sec."<<endl;
    
    //initialization phase
    param.OUTPUT_FILENAME = "kmers.ust.spss";
    globalStatFile.open("stat_ust.txt", std::fstream::out);
    
    stat.V_bcalm = adjList.size();
    nodeSign = new bool[stat.V_bcalm];
    oldToNew = new new_node_info_t[stat.V_bcalm];
    sortStruct = new struct node_sorter[stat.V_bcalm];
    
    for (unitig_t i = 0; i < stat.V_bcalm; i++) {
        nodeSign[i] = false;
        disSet.make_set(i);
        oldToNew[i].finalWalkId = -1;
        sortStruct[i].sortkey = 0;
        sortStruct[i].node = i;
    }
    
    if(param.PROFILE_AND_STAT){
        stat.E_bcalm = 0;
        for (unitig_t i = 0; i < V_bcalm; i++) { //count total number of edges
            stat.E_bcalm += adjList[i].size();
        }
    }
    
    stat.nKmers = 0;
    stat.C_bcalm = 0;
    for (unitig_struct_t unitig : unitigs) {    //count total number of kmers and characters
        stat.C_bcalm += unitig.ln;
        stat.nKmers +=  unitig.ln - K + 1;
    }
    stat.V_bcalm = V_bcalm;
    
    
    double time_a = readTimer();
    
    //DFS
    if(param.VERBOSE_MODE) cout<<"## [2] Please wait while DFS step is going on... "<<endl;
    sorter = this->sorterMaker();
    
    time_a = readTimer();
    
    
    // if(runFlag==1){
    //     return;
    // }
    
    if(param.VERBOSE_MODE) cout<<"## [3] Writing strings to file.... "<<endl;
    this->ustOutputToDisk(sorter);
    
    
    if(param.VERBOSE_MODE) cout<<"Done. TIME to output: "<<readTimer() - time_a<<" sec."<<endl;
    double TIME_TOTAL_SEC = readTimer() - startTime;
    
    
    if(param.VALIDATE){
        cout<<"## [4] Validating UST...\n";
        if(system((param.DSK_PATH +"dsk -file "+param.UNITIG_FILE+" -kmer-size "+to_string(K)+" -out list_reads.unitigs.h5 -abundance-min 1  -verbose 0").c_str())!=0) exit(3);
        if(system((param.DSK_PATH + "dsk -file absorbDecompressed.fa -kmer-size "+to_string(K)+" -abundance-min 1  -verbose 0").c_str())!=0) exit(3);
        if(system((param.DSK_PATH+"dsk2ascii -file list_reads.unitigs.h5 -out output-bcalm.txt  -verbose 0").c_str())!=0) exit(3);
        if(system((param.DSK_PATH + "dsk2ascii -file absorbDecompressed.h5 -out output-my.txt   -verbose 0").c_str())!=0) exit(3);
        //cout<<"doing highly  accurate validation................"<<endl;
        if(system("sort -k 1 -n output-bcalm.txt -o a.txt; sort -k 1 -n output-my.txt -o b.txt")!=0) exit(3);
        if(system("cmp a.txt b.txt && echo '### SUCCESS: Files Are Identical! ###' || echo '### WARNING: Files Are Different! ###'")!=0) exit(3);
        if(system("rm -rf a.txt b.txt output-bcalm.txt output-my.txt list_reads.unitigs.h5 absorbDecompressed.h5")!=0) exit(3);
    }else{
        //cout<<"## [4] Skipping output validation (please validate using validate.sh in scripts folder)\n";
    }
    

    float percent_saved_c = (1-(stat.C_ust*1.0/stat.C_bcalm))*1.0;
    float ustitchBitsPerKmer = stat.C_ust*2.0/stat.nKmers;
    
    stat.statPrinter(globalStatFile, "V_UST", stat.V_ust);
    stat.statPrinter(globalStatFile, "C_UST", stat.C_ust);
    stat.statPrinter(globalStatFile, "NT_PER_KMER_UST", ustitchBitsPerKmer/2.0);
    
    stat.statPrinter(globalStatFile, "TIME_TOTAL_SEC", TIME_TOTAL_SEC);

    globalStatFile.close();
    
    //cout << "------------ UST completed successfully.----------------- " << endl;
    //Output is in file "<<param.OUTPUT_FILENAME << " ------------"<<endl;
    cout << "Total number of unique "<<K<<"-mers " <<  " = " << stat.nKmers << endl;
    //cout<<"~~~"<<endl;
    cout << "Size of unitig based representation" <<  " = " <<stat.C_bcalm*1.0/stat.nKmers << " nucleotide/k-mer"<< endl;
    cout << "No. of strings in unitig based representation" <<  " = " <<stat.V_bcalm << endl;
    //cout<<"~~~"<<endl;
    cout << "Size of UST representation" <<  " = " <<ustitchBitsPerKmer/2.0 << " nucleotide/k-mer"<< endl;
    cout << "No. of strings in UST representation" <<  " = " <<stat.V_ust << endl;

    //cout << "------------------------------------------------------"<<endl;
}


void UST::readUnitigFile(const string& unitigFileName, vector<unitig_struct_t>& unitigs, vector<vector<edge_t> >& adjList)
{
    //bcalm_file_type = 0
    
    ifstream unitigFile;
    if(!unitigFile.good()){
        fprintf(stderr, "Error: File named \"%s\" cannot be opened.\n", unitigFileName.c_str());
        exit(EXIT_FAILURE);
    }
    
    unitigFile.open(unitigFileName);
    string line;
    
    unitig_t nodeNum;
    char lnline[20];
    char kcline[20];
    char kmline[20];
    char edgesline[100000];
    bool doCont = false;
    
    int smallestK = 9999999;
    
    {
        getline(unitigFile, line);
        do {
            unitig_struct_t unitig_struct;
            
            if(FLG_ABUNDANCE){
                //>3 LN:i:24 ab:Z:10 10 10 10 10 7 7 7 7 7 3 3 3 3   L:-:0:+ L:-:1:+  L:+:0:-
                edgesline[0] = '\0';
                sscanf(line.c_str(), "%*c %d %s", &unitig_struct.serial, lnline);
                
                if(    line.find("ab:Z") == string::npos){
                    cout<<"Incorrect input format. Check that input file matches the format of example cdbg."<<endl;
                    exit(3);
                }
                
                sscanf(lnline, "%*5c %llu", &unitig_struct.ln);
                
                int abpos = line.find("ab") + 5;
                int Lpos = line.find("L:");
                
                if(Lpos < 0){
                    Lpos = line.length() ;
                }
                // initialize string stream
                //cout<<line.substr(abpos, Lpos - abpos);
                stringstream ss(line.substr(abpos, Lpos - abpos));
                string abun;
                
                sscanf(line.substr(Lpos, line.length() - Lpos).c_str(), "%[^\n]s", edgesline);
                
                if(unitig_struct.ln < smallestK){
                    smallestK = unitig_struct.ln ;
                }
                if(unitig_struct.ln < K){
                    printf("Wrong k! Try again with correct k value. \n");
                    exit(2);
                }
                
            }else{
                edgesline[0] = '\0';
                sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &unitig_struct.serial, lnline, kcline, kmline, edgesline);
                
                if(    line.find("KC") == string::npos){
                    cout<<"Incorrect input format. Try using flag -a 1."<<endl;
                    exit(3);
                }
                
                //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
                sscanf(lnline, "%*5c %llu", &unitig_struct.ln);
                
                
                if(unitig_struct.ln < smallestK){
                    smallestK = unitig_struct.ln ;
                }
                if(unitig_struct.ln < K){
                    printf("Wrong k! Try again with correct k value. \n");
                    exit(2);
                }
            }
            
            char c1, c2;
            stringstream ss(edgesline);
            
            vector<edge_t> edges;
            while (getline(ss, line, ' ')) {
                if (delSpaces(line).length() != 0) {
                    if(DBGFLAG==VERIFYINPUT){
                        cout<<line<<endl;
                    }
                    
                    sscanf(line.c_str(), "%*2c %c %*c %d  %*c  %c", &c1, &nodeNum, &c2); //L:-:0:-
                    edge_t newEdge;
                    
                    bool DELSELFLOOP=true;
                    if(DELSELFLOOP){
                        if((unitig_struct.serial)!= nodeNum){
                            newEdge.left = charToBool(c1);
                            newEdge.right = charToBool(c2);
                            newEdge.toNode = nodeNum;
                            edges.push_back(newEdge);
                        }
                    }else{
                        newEdge.left = charToBool(c1);
                        newEdge.right = charToBool(c2);
                        newEdge.toNode = nodeNum;
                        edges.push_back(newEdge);
                    }
                }
            }
            adjList.push_back(edges);
            
            doCont = false;
            while (getline(unitigFile, line)) {
                if (line.substr(0, 1).compare(">")) {
                    //unitig_struct.sequence = unitig_struct.sequence + line;
                    unitigs.push_back(unitig_struct);
                } else {
                    doCont = true;
                    break;
                }
            }
        } while (doCont);
        unitigFile.close();
    }
    
    
    if(smallestK > K ){
        cout<<"\n :::: :::: :::: :::: !!!!!!!!! WARNING !!!!!!!!!!! :::: :::: :::: ::::"<<endl;
        cout<<"The length of the smallest string we found was " << smallestK << ". Please make sure you are using the correct value of 'k' to ensure correctness of output."<<endl;
        //cout << "------------------------------------------------------"<<endl;
    }
    //cout << "Complete reading input unitig file (bcalm2 file)." << endl;
}



void UST::collectInput(int argc, char** argv, string & graph_file_name, int & K, bool & FLG_ABUNDANCE){
    char * nvalue;
    char c;
    while( ( c = getopt (argc, argv, "i:k:a:") ) != -1 )
    {
        switch(c)
        {
            case 'a':
                if(optarg) {
                    if(strcmp(optarg, "0")==0 || strcmp(optarg, "1")==0){
                        FLG_ABUNDANCE = static_cast<bool>(std::atoi(optarg));
                    }else{
                        fprintf(stderr, "Error: use either -a 0 or -a 1 \n");
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
                    fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                            argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                        argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    if(K==0 || strcmp(nvalue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }
    graph_file_name = string(nvalue);
}


vector<UST::threetuple> UST::sorterMaker() {
    bool* saturated = new bool[V_bcalm];
    char* color = new char[V_bcalm];
    unitig_t * p_dfs = new unitig_t[V_bcalm];
    vector<list<unitig_t> > newToOld;
    double time_a = readTimer();
    
    for (unitig_t i = 0; i < V_bcalm; i++) {
        color[i] = 'w';
        p_dfs[i] = -1;
        saturated[i] = false;
    }
    
    if(param.VERBOSE_MODE) cout<<"Time to loop over all unitigs in DFS: "<<readTimer() - time_a<<" sec"<<endl;
    time_a = readTimer();
    
    for (unitig_t j = 0; j < V_bcalm; j++) {
        unitig_t u;
        
        if(0==1){
            u = sortStruct[j].node;
        }else{
            u = j;
        }
        
        if (color[u] == 'w') {  //DFS_visit(u)
            unordered_map<unitig_t, vector<edge_t> > sinkSrcEdges; //int is the unitig id (old id)
            stack<edge_t> s;
            edge_t uEdge;
            uEdge.toNode = u;
            s.push(uEdge);
            
            while (!s.empty()) {
                edge_t xEdge = s.top();
                
                unitig_t x = xEdge.toNode;
                s.pop();
                
                if (color[x] == 'w') {
                    color[x] = 'g';
                    s.push(xEdge);
                    vector<edge_t> adjx = adjList.at(x);
                    
                    // Now our branching code ::
                    // For a white x
                    // Consider 2 case:
                    // Case 1. p[x] = -1, it can happen in two way, x is the first one ever in this connected component, or no one wanted to take x
                    // either way, if p[x] = -1, i can be representative of a new node in new graph
                    // Case 2. p[x] != -1, so x won't be the representative/head of a newHome. x just gets added to its parent's newHome.
                    uint64_t u = unitigs.at(x).ln; //unitig length
                    assert(u >= K);
                    
                    if (p_dfs[x] == -1) {
                        list<unitig_t> xxx;
                        xxx.push_back(x);
                        newToOld.push_back(xxx);
                        oldToNew[x].finalWalkId = countNewNode++; // countNewNode starts at 0, then keeps increasing
                        oldToNew[x].pos_in_walk = 1;
                    } else {
                        newToOld[oldToNew[p_dfs[x]].finalWalkId].push_back(x);
                        oldToNew[x].finalWalkId = oldToNew[p_dfs[x]].finalWalkId;
                        
                        { //ALGOMODE==TWOWAYEXT || ALGOMODE==BRACKETCOMP
                            disSet.Union(x, p_dfs[x]);
                        }
                        
                        oldToNew[x].pos_in_walk = oldToNew[p_dfs[x]].pos_in_walk + 1;
                    }
                    
                    // x->y is the edge, x is the parent we are extending
                    for (edge_t yEdge : adjx) { //edge_t yEdge = adjx.at(i);
                        unitig_t y = yEdge.toNode;
                        
                        if (color[y] == 'w') { //Normal DFS
                            s.push(yEdge);
                        }
                        
                        if (y == x) { // self-loop
                            cout<<"FAIL: should not have self-loop."<<endl;
                            assert(false);
                            edge_both_t e;
                            e.edge = yEdge;
                            e.fromNode = x;
                        } else if (saturated[x]) {
                            // Since x is saturated, we only add resolveLater edges
                            // no need to check for consistency
                            if (y != p_dfs[x]) {
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                            }
                        } else {
                            // If x has space to take a child, meaning x is not saturated
                            // hunting for potential child
                            
                            if (color[y] == 'w' && p_dfs[y] == -1) {
                                // y has white color & got no parent => means it's homeless, so let's see if we can take it as a child of x
                                //But just see if it is eligible to be a child, i.e. is it consistent (sign check)?
                                
                                //2 case, Does x's child have grandparent?
                                // If No:
                                if (p_dfs[x] == -1) {
                                    // case 1: child has no grandparent
                                    // so extend path without checking any sign
                                    nodeSign[x] = yEdge.left;
                                    nodeSign[y] = yEdge.right;
                                    p_dfs[y] = x;
                                    saturated[x] = true; //found a child
                                    
                                } else if (nodeSign[x] == yEdge.left) {
                                    // case 2: child (=y) has grandparent, i.e. x's parent exists
                                    nodeSign[y] = yEdge.right;
                                    p_dfs[y] = x;
                                    saturated[x] = true; //found a child
                                    
                                } else {
                                    // do we reach this case?
                                    edge_both_t e;
                                    e.edge = yEdge;
                                    e.fromNode = x;
                                }
                                
                            } else {
                                //merger
                                {
                                    // y is not white
                                    bool consistentEdge = (nodeSign[y] == yEdge.right && (p_dfs[x]==-1 || (p_dfs[x]!=-1&& nodeSign[x] == yEdge.left)) );
                                    if(p_dfs[y]==-1 && consistentEdge && oldToNew[x].finalWalkId != oldToNew[y].finalWalkId){
                                        
                                        //cout<<"x: "<<x<<":" <<disSet.find_set(x)<<" ";
                                        //cout<<"y: "<<y<<":" <<disSet.find_set(y) <<endl;
                                        
                                        //not in same group already, prevent cycle
                                        if(disSet.find_set(x)!=disSet.find_set(y)){
                                            nodeSign[x] = yEdge.left;
                                            nodeSign[y] = yEdge.right;
                                            p_dfs[y] = x;
                                            saturated[x] = true; //found a child
                                            // oldToNew[y].serial
                                            
                                            disSet.Union(x, y);
                                            gmerge.connectGroups(oldToNew[x].finalWalkId,oldToNew[y].finalWalkId );
                                            
                                        }
                                    }
                                }
                                
                                if (y != p_dfs[x]) {
                                    edge_both_t e;
                                    e.edge = yEdge;
                                    e.fromNode = x;
                                }
                            }
                        }
                    }
                } else if (color[x] == 'g') {
                    color[x] = 'b';
                }
            }
        }
    }
    delete [] color;
    delete [] p_dfs;
    delete [] saturated;
    if(param.VERBOSE_MODE) cout<<"DFS time: "<<readTimer() - time_a<<" sec"<<endl;
    
    
    /***MERGE START***/
    bool* merged = new bool[countNewNode];  // now we will do union-find with path compresison for both way merge
    for (unitig_t i = 0; i<countNewNode; i++) {
        merged[i] = false;
    }
    
    { //both way merging
        for ( const auto& p: gmerge.fwdWalkId)
        {
            if(gmerge.fwdVisited[p.first] == false){
                unitig_t fromnode =p.first;
                unitig_t tonode = p.second;
                deque<unitig_t> lst;
                
                lst.push_back(fromnode);
                lst.push_back(tonode);
                
                gmerge.fwdVisited[fromnode] = true;
                gmerge.bwdVisited[tonode] = true;
                
                if(gmerge.fwdVisited.count(tonode)>0){
                    while(gmerge.fwdVisited[tonode] == false){
                        gmerge.fwdVisited[tonode] = true;
                        tonode = gmerge.fwdWalkId[tonode];
                        gmerge.bwdVisited[tonode] = true;
                        
                        lst.push_back(tonode);
                        if(gmerge.fwdVisited.count(tonode)==0)
                            break;
                    }
                }
                if(gmerge.bwdVisited.count(fromnode)>0){
                    while(gmerge.bwdVisited[fromnode] == false){
                        gmerge.bwdVisited[fromnode] = true;
                        fromnode = gmerge.bwdWalkId[fromnode];
                        gmerge.fwdVisited[fromnode] = true;
                        
                        lst.push_front(fromnode);
                        if(gmerge.bwdVisited.count(fromnode)==0)
                            break;
                    }
                }
                
                assert(!lst.empty());
                unitig_t commonWalkId = lst.at(0);
                unitig_t posOffset = 1;
                
                for(auto i: lst){
                    // i is new walk id before merging
                    merged[i] = true;
                    // travesing the walk list of walk ID i
                    for(unitig_t uid: newToOld[i]){
                        oldToNew[uid].finalWalkId = commonWalkId;
                        oldToNew[uid].pos_in_walk = posOffset++;
                    }
                }
                oldToNew[newToOld[lst.back()].back()].isWalkEnd = true;
                stat.V_ust++;
            }
        }
        
        for (unitig_t newNodeNum = 0; newNodeNum<countNewNode; newNodeNum++){
            if(merged[newNodeNum] == false){
                oldToNew[newToOld[newNodeNum].back()].isWalkEnd = true;
                stat.V_ust++;
            }
        }
        delete [] merged;
        vector<list<unitig_t> >().swap(newToOld);
    }
    
    //sorter of all walks and printing them
    vector<threetuple> sorter;
    for(unitig_t uid = 0 ; uid< stat.V_bcalm; uid++){
        new_node_info_t nd = oldToNew[uid];
        //sorter.push_back(make_tuple(uid, nd.finalWalkId, nd.pos_in_walk, nd.isTip));
        sorter.push_back(make_tuple(uid, nd.finalWalkId, nd.pos_in_walk));
    }
    
    stable_sort(sorter.begin(),sorter.end(),sort_by_pos);
    stable_sort(sorter.begin(),sorter.end(),sort_by_walkId);
    return sorter;
}

void UST::ustOutputToDisk(vector<threetuple>& sorter){
    string uidSeqFilename = "uidSeq.usttemp"; //"uidSeq"+ mapmode[ALGOMODE] +".txt"
    ofstream uidSequence;
    uidSequence.open(uidSeqFilename);
    
    unitig_t finalUnitigSerial = 0;
    for(threetuple n : sorter){
        unitig_t uid = get<0>(n);
        unitig_t bcalmid = unitigs.at(uid).serial;
        uidSequence << finalUnitigSerial <<" "<< bcalmid << endl;
        finalUnitigSerial++;
    }
    uidSequence.close();
    
    //keep the sequences only
    if(system(("awk '!(NR%2)' "+param.UNITIG_FILE+" > seq.usttemp").c_str())!=0) exit(3);
    if(system("sort -n -k 2 -o uidSeq.usttemp uidSeq.usttemp")!=0) exit(3);
    if(FLG_ABUNDANCE){
        if(system(("awk '(NR%2)' "+param.UNITIG_FILE+" | cut -f 5 -d ':' | cut -f 1 -d 'L' > count.usttemp").c_str())!=0) exit(3); // get a separate count file
        if(system("paste -d' ' uidSeq.usttemp seq.usttemp count.usttemp > merged.usttemp ")!=0) exit(3);
        if(system("sort -n -k 1 -o merged.usttemp merged.usttemp")!=0) exit(3);
        if(system(("cat  merged.usttemp  | awk '{for (i=4;i<=NF;i+=1) print $i}' > "+getFileName(param.UNITIG_FILE)+".ust.counts").c_str())!=0) exit(3);
    }else{
        if(system("paste -d' ' uidSeq.usttemp seq.usttemp > merged.usttemp ")!=0) exit(3);
        if(system("sort -n -k 1 -o merged.usttemp merged.usttemp")!=0) exit(3);
    }
    if(system("cat  merged.usttemp  | cut -d' ' -f3 >  seq.usttemp")!=0) exit(3);
    
    
    ifstream sequenceStringFile ("seq.usttemp");
    ofstream ustOutputFile (param.OUTPUT_FILENAME);
    
    
    //both string and abundance sort
    //keep string only and output
    //open the string file
    
    unitig_t lastWalk = -1;
    string walkString = "";
    string unitigString = "";
    for(threetuple n : sorter){
        unitig_t uid = get<0>(n);
        unitig_t finalWalkId = get<1>(n);
        
        //for each line in file
        string sequenceFromFile = "";//getline
        getline (sequenceStringFile,sequenceFromFile);
        if(nodeSign[uid] == false){
            unitigString =  reverseComplement(sequenceFromFile);
        }else{
            unitigString =  sequenceFromFile;
        }
        
        if(finalWalkId!=lastWalk){
            if(lastWalk != -1){
                //print previous walk
                if(walkString.length()>=K){
                    ustOutputFile<<">"<<endl;
                    stat.C_ust+=walkString.length();
                    ustOutputFile<< walkString<<endl;
                }
            }
            
            //start a new walk
            walkString = "";
            lastWalk = finalWalkId;
        }
        walkString = plus_strings(walkString, unitigString, K);
        //ustOutputFile<<">"<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<endl;
    }
    if(walkString.length()>=K){
        ustOutputFile<<">"<<endl;
        stat.C_ust+=walkString.length();
        ustOutputFile<< walkString<<endl;
    }else{
    }
    
    sequenceStringFile.close();
    //smallKFile.close();
    if(system("rm -rf *.usttemp")!=0) exit(3);
    if(system("rm -rf uidSeq.txt")!=0) exit(3);
    
    // UST (TWOWAYEXT) DONE!
    //delete []  global_issinksource;
}



UST::~UST(){
    delete [] nodeSign;
    delete [] oldToNew;
    delete [] sortStruct;
}


#endif /* ust_h */
