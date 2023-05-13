//
//  profileGraph.h
//  UST
//
//  Created by Amatur Rahman on 2/7/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//

#ifndef profileGraph_h
#define profileGraph_h


#include<tuple>
#include<vector>
#include<string>
#include<map>
#include<unordered_map>

#include "param.hpp"
#include "stat.hpp"
#include "spss.hpp"


using namespace std;

class ProfileGraph : public SPSS{
private: //different var/functions in UST and ESS
    unitig_t &V_bcalm = stat.V_bcalm;
    ProfileStat stat;
    Param param;
    void readUnitigFile(const string& unitigFileName, vector<unitig_struct_t>& unitigs, vector<vector<edge_t> >& adjList);
    void indegreePopulateAndLowerBound();
    unitig_t pullOutConnectedComponent(bool);
public:
    void run(string graph_file_name, int K);
};


/********************************************//**
*  ... Count the n. connected component in cdBG
***********************************************/
unitig_t ProfileGraph::pullOutConnectedComponent(bool pullit = false)
{
    unitig_t numCC_bcalm_graph = 0;
    bool *ccVisited = new bool[adjList.size()];
    vector<unitig_t> verticesInOneCC;
    map<unitig_t, unitig_t> uidRemapper;

    for (unitig_t i = 0; i < adjList.size(); i++)
    {
        ccVisited[i] = (false);
    }

    stack<unitig_t> stack;
    for (unitig_t ii = 0; ii < adjList.size(); ii++)
    {
        if (ccVisited[ii] == false)
        {
            stack.push(ii);
            numCC_bcalm_graph++;
            while (!stack.empty())
            {
                unitig_t s = stack.top();
                stack.pop();
                if (!ccVisited[s])
                {
                    //cout << s << " ";
                    ccVisited[s] = true;
                }

                for (auto i = adjList[s].begin(); i != adjList[s].end(); ++i){
                    int tonode = (*i).toNode;
                    if (!ccVisited[tonode])
                        stack.push(tonode);
                }
            }
        }
    }

   

    if (pullit)
    {
        //>1 LN:i:41 KC:i:178 km:f:16.2  L:+:2:-  L:-:1112198:+ L:-:1331912:+
        //ATCCTGAATATGGTTTTGAAAAAAACGCGGGTTACGGTACC
        ofstream subgraphFile;
        subgraphFile.open("subgraph.fa");
        for (unitig_t uid : verticesInOneCC)
        {
            subgraphFile << ">" << uidRemapper[uid] << " "
                         << "LN:i:" << unitigs.at(uid).ln << " KC:i:178 km:f:16.2 ";
            for (edge_t e : adjList[uid])
            {
                if (uidRemapper.count(e.toNode) > 0)
                {
                    subgraphFile << "L:" << boolToCharSign(e.left) << ":" << uidRemapper[e.toNode] << ":" << boolToCharSign(e.right) << " ";
                }
            }
            subgraphFile << endl;
//            subgraphFile << unitigs.at(uid).sequence << endl;
        }
        subgraphFile.close();
    }

    //cout << "number of connected components in bcalm graph: " << numCC_bcalm_graph << endl;
    delete[] ccVisited;
    return numCC_bcalm_graph;
}

void ProfileGraph::run(string graph_file_name, int K){
    cout<<"ESS-Compress v2.0 (cdBG profiler)"<<endl;
    cout<<"---------------------------------"<<endl;
    
    this->FLG_ABUNDANCE = false;
    this->K= K;
    param.UNITIG_FILE = graph_file_name;
    //collectInput(argc, argv, graph_file_name, K, FLG_ABUNDANCE); //input: argc, argv
    
    cout<<"--------------------"<<endl;
    cout << "## Please wait while input file (" << graph_file_name << ") is being read : k = "<<K<<endl;
    cout<<".........................."<<endl;
    double startTime = readTimer();
    
    this->readUnitigFile(param.UNITIG_FILE, unitigs, adjList); //input: graph_filename, output: last 2
    
    double TIME_READ_SEC = readTimer() - startTime;
    cout<<".........................."<<endl;
    cout<<"Done. TIME to read file "<<TIME_READ_SEC<<" sec."<<endl;
    cout<<"--------------------"<<endl;
    
    
    //TODO: ASSERT filename ends with .unitigs.fa
    string basename = getFileName(param.UNITIG_FILE).substr(0, getFileName(param.UNITIG_FILE).length()-11);
    
    //initialization phase
//    globalStatFile.open(("stat_cdbg_"+getFileName(param.UNITIG_FILE).substr(0, getFileName(param.UNITIG_FILE).length()-11)+".txt").c_str(), std::fstream::out);
    globalStatFile.open("stat_cdbg.txt", std::fstream::out);

    stat.statPrinter(globalStatFile, "CDBG_FILE", param.UNITIG_FILE);
    
    V_bcalm = adjList.size();
    
    //stat collector
    stat.E_bcalm = 0;
    for (unitig_t i = 0; i < V_bcalm; i++) { //count total number of edges
        stat.E_bcalm += adjList[i].size();
    }
    
    stat.nKmers = 0;
    stat.C_bcalm = 0;
    for (unitig_struct_t unitig : unitigs) {    //count total number of kmers and characters
        stat.C_bcalm += unitig.ln;
        stat.nKmers +=  unitig.ln - K + 1;
    }
    
    stat.V_bcalm = V_bcalm;
    
    cout<<"--------------------"<<endl;
    cout<<"## Please wait while info about lower bound is being gathered.... "<<endl;
    cout<<".........................."<<endl;
    double time_a = readTimer();
    
    indegreePopulateAndLowerBound();
    
    cout<<".........................."<<endl;
    cout<<"Done. TIME to gather information about unitig graph: "<<readTimer() - time_a<<" sec."<<endl;
    cout<<"--------------------"<<endl;
    
    stat.walkstarting_node_count = ceil((stat.sharedparent_count + stat.sink_count + stat.source_count)/2.0) + stat.isolated_node_count;
    stat.charLowerbound = stat.C_bcalm-(K-1)*(V_bcalm - stat.walkstarting_node_count*1);
    stat.upperbound = (1-((stat.C_bcalm-(K-1)*(V_bcalm - stat.walkstarting_node_count*1.0))/stat.C_bcalm))*100.0;
    stat.statPrinter(globalStatFile, "K", K);
    stat.statPrinter(globalStatFile, "N_KMER", stat.nKmers, true);
    stat.statPrinter(globalStatFile, "V_BCALM", stat.V_bcalm);
    stat.statPrinter(globalStatFile, "E_BCALM", stat.E_bcalm);
    stat.statPrinter(globalStatFile, "C_BCALM", stat.C_bcalm);
    stat.statPrinter(globalStatFile, "C_LB", stat.charLowerbound);
    stat.statPrinter(globalStatFile, "V_LB", stat.walkstarting_node_count, true);
    stat.statPrinter(globalStatFile, "N_ISOLATED",  stat.isolated_node_count);
    stat.statPrinter(globalStatFile, "N_SINK",  stat.sink_count);
    stat.statPrinter(globalStatFile, "N_SOURCE", stat.source_count);
    stat.statPrinter(globalStatFile, "N_SPECIAL", stat.sharedparent_count); //new bound
    
    unitig_t NUM_CC_CDBG = pullOutConnectedComponent();
    stat.statPrinter(globalStatFile, "NUM_CC_CDBG", NUM_CC_CDBG, true);
    
    
    //OPTIONAL;DERIVABLE
    stat.statPrinter(globalStatFile, "SPSS_LB", (stat.charLowerbound*1.0)/stat.nKmers, true);
    stat.statPrinter(globalStatFile, "PERCENT_UB", stat.upperbound/100.0);
    stat.statPrinter(globalStatFile, "PERCENT_N_SPECIAL", stat.sharedparent_count*1.0/stat.V_bcalm);
    stat.statPrinter(globalStatFile, "PERCENT_N_ISOLATED", stat.isolated_node_count*1.0/stat.V_bcalm);
    stat.statPrinter(globalStatFile, "PERCENT_N_DEADEND", (stat.sink_count+stat.source_count)*1.0/stat.V_bcalm);
    
//
    stat.statPrinter(globalStatFile, "ESS_LB", (stat.nKmers+3*stat.walkstarting_node_count+NUM_CC_CDBG*(K-4))/1.0/stat.nKmers, true);
    
    
    // Iterating the map and printing ordered values
    //    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
    //        cout << "(" << i->first.first<< ", "<< i->first.second << ")" << " := " << i->second << '\n';
    //    }
    
    //if(ALGOMODE == PROFILE_ONLY){
    //printf("\n");
    //exit(1);
    //}
    
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
    //##################################################//
    
    
    
    stat.statPrinter(globalStatFile, "TIME_READINPUT_SEC", TIME_READ_SEC);
    uint64_t maxulen = maximumUnitigLength();
    stat.statPrinter(globalStatFile, "MAX_UNITIG_LEN", maxulen);
    
    stat.statPrinter(globalStatFile, "MAXLEN_MB", maxulen*1.0/1024.0/1024.0);
    stat.statPrinter(globalStatFile, "E_MB", stat.E_bcalm*8.0/1024.0/1024.0);
    stat.statPrinter(globalStatFile, "V_MB", stat.V_bcalm*8.0/1024.0/1024.0);
    globalStatFile.close();
    
}


void ProfileGraph::readUnitigFile(const string& unitigFileName, vector<unitig_struct_t>& unitigs, vector<vector<edge_t> >& adjList)
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
                    cout<<"Incorrect input format. Try using flag -a 0."<<endl;
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
        cout << "------------------------------------------------------"<<endl;
    }
    //cout << "Complete reading input unitig file (bcalm2 file)." << endl;
}


void ProfileGraph::indegreePopulateAndLowerBound(){
    map<pair <int, int>, int> inOutCombo;
    
    int* global_indegree;
    int* global_outdegree;
    int* global_plusindegree;
    int* global_plusoutdegree;
    bool* global_issinksource;
    
    global_indegree = new int[V_bcalm];
    global_outdegree = new int[V_bcalm];
    global_plusindegree = new int[V_bcalm];
    global_plusoutdegree = new int[V_bcalm];
    global_issinksource = new bool[V_bcalm];
    bool* countedForLowerBound = new bool[V_bcalm];
    
    
    for (unitig_t i = 0; i < V_bcalm; i++) {
        global_indegree[i] = 0;
        global_outdegree[i] = 0;
        global_plusindegree[i] = 0;
        global_plusoutdegree[i] = 0;
        global_issinksource[i] = false;
        countedForLowerBound[i] = false;
    }
    
    int xc = 0;
    for(vector<edge_t> elist: adjList){
        for(edge_t e: elist){
            global_indegree[e.toNode] += 1;
            if(e.right == true){
                global_plusindegree[e.toNode] += 1;
            }
            if(e.left == true){
                global_plusoutdegree[xc] += 1;
            }
        }
        global_outdegree[xc] = elist.size();
        xc++;
    }
    
    for(int i = 0; i<5; i++){
        for(int j = 0; j<5; j++){
            inOutCombo[make_pair(i,j)] = 0;
        }
    }
    
    stat.sink_count  = 0;
    
    for(int i = 0; i<V_bcalm; i++){
        pair<int, int> a;
        a = make_pair(global_plusindegree[i], global_plusoutdegree[i]);
        inOutCombo[a] = (inOutCombo.count(a)  ? inOutCombo[a] + 1 : 1  );
        if(global_plusoutdegree[i] == 0 && global_plusindegree[i] != 0){
            stat.sink_count++;
            global_issinksource[i] = 1;
            countedForLowerBound[i] = true;
        }
        if(global_plusindegree[i] == 0 && global_plusoutdegree[i] != 0){
            stat.source_count++;
            global_issinksource[i] = 1;
            countedForLowerBound[i] = true;
        }
        if(global_indegree[i] == 0){
            global_issinksource[i] = 1;
            stat.isolated_node_count++;
        }
        if(global_indegree[i] == 1){
            stat.onecount++;
        }
    }
    
    
    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
        globalStatFile << "PERCENT_DEGREE_"<<(i->first).first << "_" << (i->first).second <<  "=" << (i->second)*100.0/V_bcalm <<"%" << endl;
    }
    
    xc = 0; // current vertex while traversing the adjacency list
    for(vector<edge_t> elist: adjList){
        int neighborCount = 0;
        int spNeighborCount[2];
        spNeighborCount[0]=0;
        spNeighborCount[1]=0;
        stack<int> countedNodes;
        set<pair<int, bool> > countedSides;
        //if(true){
        if(FLG_NEWUB == true){
            //ENDPOINT SIDE UPPER BOUND - improved
            for(edge_t e_xy: elist){    //central node: all neighbors of x
                int y = e_xy.toNode;
                vector<edge_t> adjY = adjList[y];
                bool eligibleSp = true;
                
                //pair<int, bool> pairr;
                for(edge_t e_from_y : adjY){    // check if this neighbor is speacial
                    //pairr =make_pair(e_from_y.toNode, sRight(e_xy) );
                    if(e_from_y.toNode!=xc){
                        
                        if(sRight(e_xy) == sLeft(e_from_y)){
                            eligibleSp = false;
                            break;
                        }
                    }
                }
                
                if(eligibleSp){
                    spNeighborCount[sLeft(e_xy)]++;
                }
            }
            if(spNeighborCount[0]>1){
                stat.sharedparent_count += spNeighborCount[0] - 1 ;
            }
            if(spNeighborCount[1]>1){
                stat.sharedparent_count += spNeighborCount[1] - 1 ;
            }
        }
        if(FLG_NEWUB == false){
            //ENDPOINT SIDE UPPER BOUND
            for(edge_t e_xy: elist){
                int y = e_xy.toNode;
                vector<edge_t> adjY = adjList[y];
                bool eligible = true;
                pair<int, bool> pairr;
                for(edge_t e_from_y : adjY){
                    pairr =make_pair(e_from_y.toNode, sRight(e_xy) );
                    if(e_from_y.toNode!=xc){
                        
                        if(sRight(e_xy) == sLeft(e_from_y)){
                            eligible = false;
                            break;
                        }
                        
                    }
                    
                }
                
                if(eligible){
                    neighborCount++;
                }
            }
            
            if(global_issinksource[xc] == 1){
                if(neighborCount>1){
                    stat.sharedparent_count += neighborCount - 1 ;
                }
            }else{
                if(neighborCount>2){
                    stat.sharedparent_count += neighborCount - 2 ;
                }
            }
        }
        //sharedparent_count_wrong =sharedparent_count;
        
        //if(true){
        if(1==0){
            // OLDER UPPER BOUND CALC
            
            int neighborCount = 0;
            for(edge_t e_xy: elist){
                int y = e_xy.toNode;
                
                if(!countedForLowerBound[y]){
                    vector<edge_t> adjY = adjList[y];
                    bool eligible = true;
                    for(edge_t e_from_y : adjY){
                        if(e_from_y.toNode!=xc){
                            if(sRight(e_xy) == sLeft(e_from_y) ){
                                eligible = false;
                                break;
                            }
                        }
                    }
                    if(eligible){
                        countedForLowerBound[y] = true;
                        //global_priority[y] = 4;
                        neighborCount++;
                        countedNodes.push(y);
                    }
                }
            }
            
            if(global_issinksource[xc] == 1){
                if(neighborCount>1){
                    stat.sharedparent_count += neighborCount - 1 ;
                }else{
                    while(!countedNodes.empty()){
                        countedForLowerBound[countedNodes.top()] = false;
                        countedNodes.pop();
                    }
                }
            }else{
                if(neighborCount>2){
                    stat.sharedparent_count += neighborCount - 2 ;
                }else{
                    while(!countedNodes.empty()){
                        countedForLowerBound[countedNodes.top()] = false;
                        countedNodes.pop();
                    }
                }
            }
        }
        
        xc++;
    }
    
    delete [] global_indegree;
    delete [] global_outdegree;
    delete [] global_plusindegree;
    delete [] global_plusoutdegree;
    delete [] countedForLowerBound;
}


#endif /* profileGraph_h */
