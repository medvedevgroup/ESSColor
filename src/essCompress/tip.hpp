//
//  tip.cpp
//  UST
//
//  Created by Amatur Rahman on 21/5/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//
#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include <cmath>
#include <cstring>
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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include<list>

#include "spss.hpp"
#include "decoder.hpp"

#include "param.hpp"


using namespace std;

class ESSTip : public SPSS{
public:
    typedef struct {
        bool isWalkEnd = false;
        unitig_t pos_in_walk = -100;
        unitig_t finalWalkId = -1; // renders some walkId as invalid
        int isTip = 0;
    } new_node_info_t;
    new_node_info_t* oldToNew;
    unitig_t &V_bcalm = stat.V_bcalm;
    
    unitig_t countNewNode = 0;
    bool* nodeSign;
    ESSTipStat stat;
    ESSTipParam param;
    
    bool* global_issinksource;
    
    typedef tuple<unitig_t,unitig_t,unitig_t, int> fourtuple; // uid, walkid, pos, isTip
    bool static sort_by_walkId (const fourtuple &lhs, const fourtuple &rhs){
        return get<1>(lhs) < get<1>(rhs);
    }
    bool static sort_by_pos (const fourtuple &lhs, const fourtuple &rhs){
        return get<2>(lhs) < get<2>(rhs);
    }
    bool static sort_by_tipstatus (const fourtuple &lhs, const fourtuple &rhs){
        return get<3>(lhs) < get<3>(rhs);
    }
    vector<fourtuple> sorter;
    
    void reassign(bool* obsoleteWalkId){
        unitig_t newPid = 0;
        map<unitig_t, unitig_t> reassignMap;
        for(unitig_t p = 0; p<countNewNode; p++){
            if(!obsoleteWalkId[p]){
                reassignMap[p] = newPid++;
            }
        }
        for(unitig_t u = 0; u<V_bcalm; u++){
            unitig_t pid = oldToNew[u].finalWalkId;
            oldToNew[u].finalWalkId = reassignMap[pid];
        }
        countNewNode = newPid;
    }
    
    void ustOutputToDisk(vector<fourtuple>& sorter){
        //OVERRIDE
        
        ofstream uidSequence;
        string uidSeqFilename = "uidSeq.usttemp"; //"uidSeq"+ mapmode[ALGOMODE] +".txt"
        uidSequence.open(uidSeqFilename);
        
        unitig_t finalUnitigSerial = 0;
        for(fourtuple& n : sorter){
            unitig_t uid = get<0>(n);
            unitig_t bcalmid = unitigs.at(uid).serial;
            //                        int finalWalkId = get<1>(n);
            //                        int pos_in_walk = get<2>(n);
            //                        int isTip = get<3>(n);
            //                        cout<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<" "<<oldToNew[uid].isWalkEnd<< " was merged: "<< merged[oldToNew[uid].finalWalkId]<< endl;
            uidSequence << finalUnitigSerial <<" "<< bcalmid << endl;
            finalUnitigSerial++;
        }
        uidSequence.close();
        
        
        //system();
        
        //keep the sequences only
        if(system(("awk '!(NR%2)' "+param.UNITIG_FILE+" > seq.usttemp").c_str())!=0) exit(3);
        if(system("sort -n -k 2 -o uidSeq.usttemp uidSeq.usttemp")!=0) exit(3);
        if(system("paste -d' ' uidSeq.usttemp seq.usttemp > merged.usttemp ")!=0) exit(3);
        if(system("mkdir -p tmpUSTdir")!=0) exit(3);
        if(system("sort -T tmpUSTdir/ -n -k 1 -o merged.usttemp merged.usttemp")!=0) exit(3);
        
        if(system("cat  merged.usttemp  | cut -d' ' -f3 >  seq.usttemp")!=0) exit(3);
        
        
        //string walkString = "";
        
        //param.OUTPUT_FILENAME = param.getBcalmFileBasename() + ".esstip";
        //param.OUTPUT_FILENAME = "kmers.esstip";
        ofstream tipFile(param.OUTPUT_FILENAME);
        
        
        ifstream sequenceStringFile ("seq.usttemp");
        
        //for(int si = 0; si<sorter.size(); si++){
        unitig_t startWalkIndex = 0;
        unitig_t prevwalk = -1;
        while(true){
            
            assert(startWalkIndex<sorter.size());
            fourtuple &n = sorter[startWalkIndex];
            unitig_t finalWalkId = get<1>(n);
            
            unitig_t uid = get<0>(n);
            int isTip = get<3>(n);
            unitig_t pos_in_walk =get<2>(n);
            //cout<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
            
            //            string unitigString;
            //            if(nodeSign[uid] == false){
            //                unitigString =  reverseComplement(unitigs.at(uid).sequence);
            //            }else{
            //                unitigString =  (unitigs.at(uid).sequence);
            //            }
            
            
            //
            string unitigString;
            string sequenceFromFile = "";//getline
            getline (sequenceStringFile,sequenceFromFile);
            if(nodeSign[uid] == false){
                unitigString =  reverseComplement(sequenceFromFile);
            }else{
                unitigString =  sequenceFromFile;
            }
            
            
            
            if(prevwalk!=finalWalkId){
                tipFile<<">\n";
            }
            prevwalk = finalWalkId;
            
            //            if( MODE_ABSORPTION_TIP && isTip == 0){//
            //                walkString = plus_strings(walkString, unitigString, K);
            //
            if(isTip == 0 ){
                
                if(startWalkIndex==0){
                    tipFile<<unitigString;
                    stat.C_esstip += unitigString.length();
                }else if(finalWalkId != get<1>(sorter[startWalkIndex-1])){
                    tipFile<<unitigString;
                    stat.C_esstip += unitigString.length();
                }else{
                    tipFile<<unitigString.substr(K - 1, unitigString.length() - (K - 1));
                    stat.C_esstip += unitigString.substr(K - 1, unitigString.length() - (K - 1)).length();
                //tipFile<<cutPref(unitigString, K);
                }
                
            }else if(isTip==1){ //right R   R    ]   ]   ]   ]
                //cut prefix: correct
                if(0==0){
                    unitigString = unitigString.substr(K - 1, unitigString.length() - (K - 1));
                    //                    if(walkString.length()<K){
                    //                        cout<<"pos: "<<walkString.length()<<endl;
                    //                    }
                    tipFile<<"(";
                    tipFile<<unitigString;
                    tipFile<<")";
                    stat.C_esstip += unitigString.length() + 2;
                    stat.C_nondna_esstip += 2;
                }
                if(1==0){
                    tipFile<<">pref\n"<<unitigString<<endl;
                }
                
            }else if(isTip==2){ //left L   L    [ [ [
                //cut suffix: correct
                if(0==0){
                    unitigString = unitigString.substr(0, unitigString.length() - (K - 1));
                    //                    if(walkString.length()<K){
                    //                        cout<<"pos: "<<walkString.length()<<endl;
                    //                    }
                    tipFile<<"{";
                    tipFile<<unitigString;
                    tipFile<<"}";
                    stat.C_esstip += unitigString.length() + 2;
                    stat.C_nondna_esstip += 2;
                }
                if(1==0){
                    tipFile<<">suf\n"<<unitigString<<endl;
                }
            }
            
            if(startWalkIndex+1 == sorter.size()) {
                //                int brackets1 = std::count(walkString.begin(), walkString.end(), '(');
                //                int brackets2 = std::count(walkString.begin(), walkString.end(), ')');
                //                int stringPlus = std::count(walkString.begin(), walkString.end(), '{');
                //                int stringMinus = std::count(walkString.begin(), walkString.end(), '}');
                //
                //                C_tip_special += brackets1- brackets2 -stringPlus-stringMinus;
                //                C_tip_ustitch += walkString.length();
                //
                //tipFile<<">\n"<<walkString<<endl;
                tipFile<<endl;
                stat.V_esstip++;
                break;
            }else if(get<1>(sorter[startWalkIndex+1]) != finalWalkId){
                //tipFile<<">\n"<<walkString<<endl;
                tipFile<<endl;
                
                stat.V_esstip++;
                //walkString = "";
                //walkString.shrink_to_fit();
                //break;
            }
            startWalkIndex++;
        }
        //}
        
        if(system("rm -rf tmpUSTdir *.usttemp")!=0) exit(3);
        
        sequenceStringFile.close();
        tipFile.close();
    }
    
    
    vector<ESSTip::fourtuple> sorterMaker() {
        //OVERLOAD
        bool BRACKETCOMP = true;
        
        
        vector<ESSTip::fourtuple> sorter;
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
        
        //cout<<"Basic V loop time: "<<readTimer() - time_a<<" sec"<<endl;
        time_a = readTimer();
        
        
        for (unitig_t j = 0; j < V_bcalm; j++) {
            unitig_t u;
            
            if(0==1){
                u = sortStruct[j].node;
            }else{
                u = j;
            }
            
            if( BRACKETCOMP){
                if(global_issinksource[u]==1){
                    continue;
                }
            }
            
            
            if (color[u] == 'w') {  //DFS_visit(u)
                
                unordered_map<unitig_t, vector<edge_t> > sinkSrcEdges; //int is the unitig id (old id)
                
                if( BRACKETCOMP){
                    if(global_issinksource[u]==1){
//                        vector<edge_t> adju = adjList.at(u);
//                        vector<edge_t> myvector;
//                        for (edge_t e : adju) {
//                            myvector.push_back(e);
//                        }
//                        sinkSrcEdges[u] = myvector;
                        continue;
                    }
                }
                
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
                        unitig_t u = unitigs.at(x).ln; //unitig length
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
                            
                            if( BRACKETCOMP){
                                if(global_issinksource[y] == true){
                                    continue;
                                }
                            }
                            
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
        
        delete [] p_dfs;
        delete [] saturated;
        if(param.VERBOSE_MODE) cout<<"Done. DFS time: "<<readTimer() - time_a<<" sec"<<endl;
        
        
        /***MERGE START***/
        bool* merged = new bool[countNewNode];  // now we will do union-find with path compresison for both way merge
        bool* obsoleteWalkId = new bool[countNewNode];
        
        for (unitig_t i = 0; i<countNewNode; i++) {
            merged[i] = false;
            obsoleteWalkId[i] = false;
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
                        
                        //POTENTIAL BUG @BRACKETCOMP
                        if(i!=commonWalkId){
                            obsoleteWalkId[i] = true;
                        }
                        
                        // travesing the walk list of walk ID i
                        for(unitig_t uid: newToOld[i]){
                            oldToNew[uid].finalWalkId = commonWalkId;
                            oldToNew[uid].pos_in_walk = posOffset++;
                        }
                    }
                    oldToNew[newToOld[lst.back()].back()].isWalkEnd = true;
                    stat.V_nontip ++;
                }
            }
            
            for (unitig_t newNodeNum = 0; newNodeNum<countNewNode; newNodeNum++){
                if(merged[newNodeNum] == false){
                    oldToNew[newToOld[newNodeNum].back()].isWalkEnd = true;
                    stat.V_nontip++;
                }
            }
            delete [] merged;
            vector<list<unitig_t> >().swap(newToOld);
        }
        
        
        //@@@@@ BRACKETED
        reassign(obsoleteWalkId);
        delete [] obsoleteWalkId;
        
        bool* hasStartTip = new bool[V_bcalm];
        bool* hasEndTip = new bool[V_bcalm];
        for (unitig_t i = 0; i<V_bcalm; i++) {
            hasStartTip[i] = false;
            hasEndTip[i] = false;
        }
        
        if( BRACKETCOMP){
            for (unitig_t sinksrc = 0; sinksrc<V_bcalm; sinksrc++){
                
                if(global_issinksource[sinksrc]==0){
                    continue;
                }
                
                for(edge_t e: adjList[sinksrc]){
                    
                    // when can this occur? it does occur
                    if(color[sinksrc] != 'w'){
                        break;
                    }
                    
                    //there are 3 cases
                    //if consistent this way [[[if(nodeSign[e.toNode] == e.right)]]]
                    //case fwd1: sinksrc -> contig start
                    //case fwd2. sinksrc -> contig middle/end -> ... (say that sinksrc is LEFT)
                    //case fwd3. sinksrc -> sinksrc_other (i'd say ignore this for now)
                    //
                    
                    //case bwd1. contig end -> sinksrc
                    //case bwd2. .... -> contig middle/start -> sinksrc (say that sinksrc is RIGHT)
                    //case bwd3. sinksrc_other -> sinksrc  (i'd say ignore this for now)
                    
                    // 3 fwd cases
                    if(nodeSign[e.toNode] == e.right){  //ensure it is a fwd case
                        if(color[e.toNode]!='w' && color[e.toNode]!='r' && color[e.toNode]!='l'){//  this ensures that this to vertex is NOT sinksrc_other
                            //case 1, 2
                            unitig_t whichwalk = oldToNew[e.toNode].finalWalkId;
                            //*** case fwd1 : sinksrc -> contigStart
                            //case fwd2. sinksrc -> contig middle/end -> ... (say that sinksrc is LEFT)
                            //let's merge case fwd1 & fwd2
                            //color[sinksrc] = 'b';
                            
                            nodeSign[sinksrc] = e.left;
                            color[sinksrc] = 'l';
                            oldToNew[sinksrc].finalWalkId = whichwalk;
                            oldToNew[sinksrc].pos_in_walk = oldToNew[e.toNode].pos_in_walk - 1;
                            assert(oldToNew[e.toNode].pos_in_walk != -1);
                            
                            // @@DBG_BLOCK int k = oldToNew[e.toNode].pos_in_walk;
                            // @@DBG_BLOCK bool jjjj = hasStartTip[e.toNode];
                            if(oldToNew[e.toNode].pos_in_walk == 1 && hasStartTip[e.toNode] == false ){
                                oldToNew[sinksrc].isTip = 0;
                                hasStartTip[e.toNode] = true;
                                oldToNew[sinksrc].pos_in_walk = oldToNew[e.toNode].pos_in_walk - 1 ;
                                unitig_t j = oldToNew[sinksrc].pos_in_walk;
                                
                            }else{
                                oldToNew[sinksrc].isTip = 2;
                                //assert(false);
                            }
                            
                        }
                        
                    }else{
                        // 3 bwd cases
                        
                        if((color[e.toNode]!='w' && color[e.toNode]!='r' && color[e.toNode]!='l')){
                            unitig_t whichwalk = oldToNew[e.toNode].finalWalkId;
                            
                            //*** case bwd1: contigend --> sinksrc
                            //*** case bwd2: contigmiddle--> sinksrc
                            
                            
                            nodeSign[sinksrc] = !e.left;
                            //color[sinksrc] = 'b';
                            color[sinksrc] = 'r';
                            oldToNew[sinksrc].finalWalkId = whichwalk;
                            oldToNew[sinksrc].pos_in_walk = oldToNew[e.toNode].pos_in_walk ;
                            assert(oldToNew[e.toNode].pos_in_walk != -1);
                            
                            if(oldToNew[e.toNode].isWalkEnd == true && !hasEndTip[e.toNode] ){
                                oldToNew[sinksrc].isTip = 0;
                                hasEndTip[e.toNode] = true;
                                oldToNew[sinksrc].pos_in_walk = oldToNew[e.toNode].pos_in_walk + 1;
                            }else{
                                oldToNew[sinksrc].isTip = 1;
                                //assert(false);
                            }
                        }
                    }
                    
                }
            }
            delete [] hasStartTip;
            delete [] hasEndTip;
            // now take care of all the remaining edges
            //            for (auto const& x : sinkSrcEdges)
            //            {
            //                int sinksrc = x.first;
            //                if(color[sinksrc] == 'w'){  //still white, that means it goes isolated now
            //                    list<int> xxx;
            //                    xxx.push_back(sinksrc);
            //                    newToOld.push_back(xxx);
            //                    oldToNew[sinksrc].serial = countNewNode++;
            //                    oldToNew[sinksrc].finalWalkId = oldToNew[sinksrc].serial;
            //                    oldToNew[sinksrc].pos_in_walk = 1;
            //                    oldToNew[sinksrc].isTip = 0;
            //                    // error resolved in sept 14
            //                    color[sinksrc] = 'b';
            //                }
            //            }
            
            for (unitig_t sinksrc = 0; sinksrc<V_bcalm; sinksrc++) {
                bool istipit = (oldToNew[sinksrc].isTip==1 || oldToNew[sinksrc].isTip==2);
                if(global_issinksource[sinksrc] == 1 && color[sinksrc] == 'w'){
                    list<unitig_t> xxx;
                    xxx.push_back(sinksrc);
                    newToOld.push_back(xxx);
                    
                    if(istipit)  {
                        oldToNew[sinksrc].finalWalkId = countNewNode++;
                        oldToNew[sinksrc].pos_in_walk = 1;
                    }else if(color[sinksrc] == 'r' || color[sinksrc] == 'l'  ){
                        oldToNew[sinksrc].isTip = 0;
                    }else{
                        oldToNew[sinksrc].finalWalkId = countNewNode++;
                        oldToNew[sinksrc].pos_in_walk = 1;
                        oldToNew[sinksrc].isTip = 0;
                    }
                    // error resolved in sept 14, 2019
                    color[sinksrc] = 'b';
                }
            }
            delete [] color;
            
            // make the absorb graph
            //BRACKETCOMP encoder and printer::::
            //vector<fourtuple> sorter;
            for(unitig_t uid = 0 ; uid< V_bcalm; uid++){
                ESSTip::new_node_info_t nd = oldToNew[uid];
                sorter.push_back(make_tuple(uid, nd.finalWalkId, nd.pos_in_walk, nd.isTip));
            }
            
            stable_sort(sorter.begin(),sorter.end(),sort_by_tipstatus);
            stable_sort(sorter.begin(),sorter.end(),sort_by_pos);
            stable_sort(sorter.begin(),sorter.end(),sort_by_walkId);
            
            
            
        }
        assert(sorter.size()!=0);
        
        delete [] global_issinksource;
        return sorter;
    }
    
    
    void degreePopulate(){
        unitig_t* global_indegree;
        unitig_t* global_outdegree;
        unitig_t* global_plusindegree;
        unitig_t* global_plusoutdegree;
        //bool* global_issinksource;
        
        global_indegree = new unitig_t[V_bcalm];
        global_outdegree = new unitig_t[V_bcalm];
        global_plusindegree = new unitig_t[V_bcalm];
        global_plusoutdegree = new unitig_t[V_bcalm];
        global_issinksource = new bool[V_bcalm];
        
        for (unitig_t i = 0; i < V_bcalm; i++) {
            global_indegree[i] = 0;
            global_outdegree[i] = 0;
            global_plusindegree[i] = 0;
            global_plusoutdegree[i] = 0;
            global_issinksource[i] = false;
        }
        unitig_t xc = 0;
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
        
        for(unitig_t i = 0; i<V_bcalm; i++){
            if(global_plusoutdegree[i] == 0 && global_plusindegree[i] != 0){
                global_issinksource[i] = 1;
            }
            if(global_plusindegree[i] == 0 && global_plusoutdegree[i] != 0){
                global_issinksource[i] = 1;
            }
            if(global_indegree[i] == 0){
                global_issinksource[i] = 1;
            }
        }
        
        delete [] global_indegree;
        delete [] global_outdegree;
        delete [] global_plusindegree;
        delete [] global_plusoutdegree;
    }
    
    

    void readUnitigFile(const string& unitigFileName, vector<unitig_struct_t>& unitigs, vector<vector<edge_t> >& adjList)
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
        
        
        // if(smallestK > K ){
        //     cout<<"\n :::: :::: :::: :::: !!!!!!!!! WARNING !!!!!!!!!!! :::: :::: :::: ::::"<<endl;
        //     cout<<"The length of the smallest string we found was " << smallestK << ". Please make sure you are using the correct value of 'k' to ensure correctness of output."<<endl;
        //     cout << "------------------------------------------------------"<<endl;
        // }

        //cout << "Complete reading input unitig file (bcalm2 file)." << endl;
    }

    
    void run(string graph_file_name, int K, bool runFlag = 0){
        param.VERBOSE_MODE = runFlag;
        //OVERRIDE
        if(param.VERBOSE_MODE) cout<<"Running ESS-Compress (TIP-mode)"<<endl;
        this->K = K;
        param.UNITIG_FILE = graph_file_name;
        //collectInput(argc, argv, graph_file_name, K, FLG_ABUNDANCE); //input: argc, argv
        
        if(param.VERBOSE_MODE) cout << "[1] Please wait input file is being loaded... " << graph_file_name << ": k = "<<K<<endl;
        double startTime = readTimer();
        
        readUnitigFile(graph_file_name, unitigs, adjList); //input: graph_filename, output: last 2
        
        double TIME_READ_SEC = readTimer() - startTime;
        if(param.VERBOSE_MODE) cout<<"Done. TIME to read file "<<TIME_READ_SEC<<" sec."<<endl;
        
        //initialization phase
        param.OUTPUT_FILENAME = "kmers.esstip";
        //param.OUTPUT_FILENAME = param.UNITIG_FILE+".essc";
        
//        globalStatFile.open(("stat_esstip_"+getFileName(param.UNITIG_FILE).substr(0, getFileName(param.UNITIG_FILE).length()-11)+".txt").c_str(), std::fstream::out);
        globalStatFile.open("stat_esstip.txt", std::fstream::out);

        
        V_bcalm = adjList.size();
        nodeSign = new bool[V_bcalm];
        oldToNew = new ESSTip::new_node_info_t[V_bcalm];
        
        sortStruct = new struct node_sorter[V_bcalm];
        
        for (unitig_t i = 0; i < V_bcalm; i++) {
            nodeSign[i] = false;
            disSet.make_set(i);
            oldToNew[i].finalWalkId = -1;
            sortStruct[i].sortkey = 0;
            sortStruct[i].node = i;
        }
        
        if(param.PROFILE_AND_STAT){
            //stat collector
            stat.E_bcalm = 0;
            for (unitig_t i = 0; i < V_bcalm; i++) { //count total number of edges
                stat.E_bcalm += adjList[i].size();
            }
            stat.statPrinter(globalStatFile, "E_BCALM", stat.E_bcalm);
        }
        
        stat.nKmers = 0;
        stat.C_bcalm = 0;
        for (unitig_struct_t unitig : unitigs) {    //count total number of kmers and characters
            stat.C_bcalm += unitig.ln;
            stat.nKmers +=  unitig.ln - K + 1;
        }
        
        
        if(param.VERBOSE_MODE) cout<<"[2] Please wait while number of dead-ends are being counted... "<<endl;
        double time_a = readTimer();
        degreePopulate();
        if(param.VERBOSE_MODE) cout<<"Done. TIME to count the nodes in unitig graph: "<<readTimer() - time_a<<" sec."<<endl;


         stat.statPrinter(globalStatFile, "K", K);
         stat.statPrinter(globalStatFile, "N_KMER", stat.nKmers);
         stat.statPrinter(globalStatFile, "V_BCALM", stat.V_bcalm);
         stat.statPrinter(globalStatFile, "C_BCALM", stat.C_bcalm);
     
         
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
         //##################################################//
         
         //DFS
         if(param.VERBOSE_MODE) cout<<"[3] DFS algorithm for stitching unitigs is running... "<<endl;
         sorter = this->sorterMaker();
        
         
         time_a = readTimer();
         
        //  if(runFlag==1){
        //      return;
        //  }
         
         if(param.VERBOSE_MODE) cout<<"[4] Writing strings to disk.... "<<endl;
         this->ustOutputToDisk(sorter);
        
        //VALIDATION
       // decodeTip(K, "fa.esstip", "fa.spss.fa");
         if(param.VALIDATE){
             cout<<"## [4] Decoding ESS...\n";
             double decodetime = readTimer();
             decodeTip(K,param.OUTPUT_FILENAME);
             decodetime = readTimer() - decodetime;
             stat.statPrinter(globalStatFile, "TIME_DECODE_SEC", decodetime, true);
             cout<<"[COMPLETE][4] Decoding done. TIME: "<<decodetime<<"sec. \n";
             cout << "------------------------------------------------------"<<endl;
             
             //validate();
             cout<<"## [5] Validating decoded ESS...\n";
             if(system((param.DSK_PATH +"dsk -file "+param.UNITIG_FILE+" -kmer-size "+to_string(K)+" -out list_reads.unitigs.h5 -abundance-min 1  -verbose 0").c_str())!=0) exit(3);
             if(system((param.DSK_PATH + "dsk -file absorbDecompressed.fa -kmer-size "+to_string(K)+" -abundance-min 1  -verbose 0").c_str())!=0) exit(3);
             if(system((param.DSK_PATH+"dsk2ascii -file list_reads.unitigs.h5 -out output-bcalm.txt  -verbose 0").c_str())!=0) exit(3);
             if(system((param.DSK_PATH + "dsk2ascii -file absorbDecompressed.h5 -out output-my.txt   -verbose 0").c_str())!=0) exit(3);
             //cout<<"doing highly  accurate validation................"<<endl;
             if(system("sort -k 1 -n output-bcalm.txt -o a.txt; sort -k 1 -n output-my.txt -o b.txt")!=0) exit(3);
             if(system("cmp a.txt b.txt && echo '### SUCCESS: Files Are Identical! ###' || echo '### WARNING: Files Are Different! ###'")!=0) exit(3);
             if(system("rm -rf a.txt b.txt output-bcalm.txt output-my.txt list_reads.unitigs.h5 absorbDecompressed.h5")!=0) exit(3);
         }
        
        
         if(param.VERBOSE_MODE) cout<<"Done. TIME to output: "<<readTimer() - time_a<<" sec."<<endl;
         double TIME_TOTAL_SEC = readTimer() - startTime;
         
         
         
         float percent_saved_c = (1-(stat.C_esstip*1.0/stat.C_bcalm))*1.0;
         float esstipCharPerKmer = stat.C_esstip*1.0/stat.nKmers;
        
         
         stat.statPrinter(globalStatFile, "V_ESSTIP", stat.V_esstip);
         stat.statPrinter(globalStatFile, "C_ESSTIP", stat.C_esstip);
        stat.statPrinter(globalStatFile, "C_NONDNA_ESSTIP", stat.C_nondna_esstip);

        stat.statPrinter(globalStatFile, "CHAR_PER_KMER_ESSTIP", esstipCharPerKmer);
        stat.statPrinter(globalStatFile, "CHAR_PER_KMER_NONDNA_ESSTIP", stat.C_nondna_esstip*1.0/stat.nKmers);
        
        
        stat.statPrinter(globalStatFile, "TIME_TOTAL_SEC", TIME_TOTAL_SEC);
        
         uint64_t maxulen = maximumUnitigLength();
         stat.statPrinter(globalStatFile, "MAX_UNITIG_LEN_MB", maxulen*1.0/1024.0/1024.0);
         stat.statPrinter(globalStatFile, "E_MB", stat.E_bcalm*8.0/1024.0/1024.0);
         stat.statPrinter(globalStatFile, "V_MB", stat.V_bcalm*8.0/1024.0/1024.0);
         
         globalStatFile.close();
         
        if(param.VERBOSE_MODE) cout << "------ End of ESS-Compress (TIP mode) core : Success! ------"<<endl;
        //Output is in file "<<param.OUTPUT_FILENAME <<
         cout << "Total number of unique "<<K<<"-mers " <<  "= " << stat.nKmers << endl;
        
         cout << "Size of ESS-Tip-Compress representation" <<  "= " <<esstipCharPerKmer << " char/k-mer"<< endl;

    }
};
    
