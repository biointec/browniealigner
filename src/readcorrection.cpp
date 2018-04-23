 
/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
 *   This file is part of Brownie                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <thread>
#include <string>
#include <iomanip>

#include "library.h"
#include "readcorrection.h"
#include <algorithm>

using namespace std;

// ============================================================================
// ALIGNMENT METRICS CLASS
// ============================================================================

void AlignmentMetrics::addMetrics(const AlignmentMetrics& rhs)
{
        lock_guard<mutex> lock(metricMutex);
        
        numCorrReads += rhs.numCorrReads;
        numCorrByMEM += rhs.numCorrByMEM;
        numSubstitutions += rhs.numSubstitutions;
        
        if (rhs.maxReadLen > maxReadLen)
                maxReadLen = rhs.maxReadLen;
        if (numReads + rhs.numReads > 0)
                avgReadLen = (double) ( avgReadLen * numReads + rhs.numReads * rhs.avgReadLen)/(double)(numReads+ rhs.numReads);
        numReads += rhs.numReads;
        
        
}
void AlignmentMetrics::updateMetrics(size_t &tootalNumOfReads, size_t &numOfSub, size_t &numCorrectedByMEM ,size_t & numCorrectedReads)
{
        lock_guard<mutex> lock(metricMutex);
        numCorrReads += numCorrectedReads;
        numCorrByMEM += numCorrectedByMEM;
        numSubstitutions += numOfSub;
        numReads = tootalNumOfReads;
        
        
}

void AlignmentMetrics::printStatistics() const
{
        size_t numCorrByKmer = numCorrReads - numCorrByMEM;
        size_t numUncorrected = numReads - numCorrReads;

        cout << "\nRead alignment report:\n";
        cout << "\tNumber of reads handled: " << numReads << endl;
        cout << "\tNumber of aligned reads: " << numCorrReads
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numCorrReads, numReads) << "%)" << endl;
        cout << "\tNumber of reads aligned by kmer seeds: " << numCorrByKmer
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numCorrByKmer, numReads) << "%)" << endl;
        cout << "\tNumber of reads aligned by MEM seeds: " << numCorrByMEM
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numCorrByMEM, numReads) << "%)" << endl;
        cout << "\tNumber of substitutions in reads: " << numSubstitutions
             << fixed << setprecision(2) << " (avg of "
             << double(numSubstitutions)/double(numReads) << " per read)" << endl;
        cout << "\tNumber of unaligned reads: " << numUncorrected
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numUncorrected, numReads) << "%)" << endl;
 
}

// ============================================================================
// SEED CLASS
// ============================================================================

void Seed::createNodePosition(const DBGraph& dbg, vector<NodePosPair>& npp) const
{
        size_t currReadPos = readFirst;
        for (vector<NodeID>::const_iterator it = nodeID.begin(); it != nodeID.end(); it++) {
                SSNode node = dbg.getSSNode(*it);
                size_t nodeOffset = (it == nodeID.begin()) ? nodeFirst : 0;
                size_t nodeOL = min(node.getMarginalLength() - nodeOffset, readEnd - currReadPos);
                for (size_t i = 0; i < nodeOL; i++)
                        npp[currReadPos+i] = NodePosPair(*it, i + nodeOffset);
                currReadPos += nodeOL;
        }
}

void Seed::mergeSeeds(const vector<Seed>& seeds,
                      vector<Seed>& mergedSeeds)
{
        if (seeds.empty())
                return;

        mergedSeeds.reserve(seeds.size());
        mergedSeeds.push_back(seeds.front());

        for (size_t i = 1; i < seeds.size(); i++) {
                const Seed& left = seeds[i-1];
                const Seed& right = seeds[i];

                bool consistent = false;
                if (left.nodeID.size() == 1 && right.nodeID.size() == 1)
                        if (left.nodeID.front() == right.nodeID.front())
                                if ((right.nodeFirst - left.nodeFirst) == (right.readFirst - left.readFirst))
                                        consistent = true;

                if (consistent)
                        mergedSeeds.back().readEnd = right.readEnd;
                else
                        mergedSeeds.push_back(right);
        }
}

// ============================================================================
// READ CORRECTION CLASS
// ============================================================================

void ReadCorrection::revCompl(vector< NodePosPair >& npp)
{
        reverse(npp.begin(), npp.end());

        for (size_t i = 0; i < npp.size(); i++) {
                if (!npp[i].isValid())
                        continue;

                NodeID nodeID = npp[i].getNodeID();
                size_t pos = npp[i].getPosition();
                const SSNode node = dbg.getSSNode(nodeID);
                npp[i] = NodePosPair(-nodeID, node.getMarginalLength() - 1 - pos);
        }
}

void ReadCorrection::findNPPSlow(const string& read, vector<NodePosPair>& npp)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair result = dbg.findNPP(kmer);
                npp[it.getOffset()] = result;
        }
}


void ReadCorrection::findNPPFast(const string& read, vector<NodePosPair>& nppv)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair npp = dbg.findNPP(kmer);
                nppv[it.getOffset()] = npp;

                if (!npp.isValid())
                        continue;

                NodeID nodeID = npp.getNodeID();
                const SSNode node = dbg.getSSNode(nodeID);

                size_t readPos = it.getOffset() + Kmer::getK();
                size_t nodePos = npp.getPosition() + Kmer::getK();

                while ((readPos < read.size()) && (nodePos < node.getLength())) {
                        if (read[readPos] == node.getNucleotide(nodePos)) {
                                it++;
                                nppv[it.getOffset()] = NodePosPair(nodeID, nodePos - Kmer::getK() + 1);
                                nodePos++;
                                readPos++;
                        } else
                                break;

                }
        }
}

void ReadCorrection::extractSeeds(const vector<NodePosPair>& nppv,
                                     vector<Seed>& seeds)
{
        size_t prev = nppv.size();
        for (size_t i = 0; i < nppv.size(); i++) {
                if (!nppv[i].isValid())
                        continue;

                // is it the first time we encounter a valid npp?
                if (prev == nppv.size()) {
                        prev = i;
                        seeds.push_back(Seed(nppv[i].getNodeID(), nppv[i].getPosition(), i, i + 1));
                        continue;
                }

                bool consistent = false;
                const SSNode thisNode = dbg.getSSNode(nppv[i].getNodeID());
                const SSNode prevNode = dbg.getSSNode(nppv[prev].getNodeID());

                // no, check for check for consistency
                if (thisNode.getNodeID() == prevNode.getNodeID()) {
                        if ((nppv[i].getPosition() - nppv[prev].getPosition()) == (i - prev))
                                consistent = true;
                } else {                // we're in different nodes
                        if (prevNode.getRightArc(thisNode.getNodeID()) != NULL) {
                                size_t thisPos = prevNode.getMarginalLength() + nppv[i].getPosition();
                                if ((thisPos - nppv[prev].getPosition()) == (i - prev)) {
                                        seeds.back().nodeID.push_back(thisNode.getNodeID());
                                        consistent = true;
                                }
                        }
                }

                prev = i;
                if (!consistent) {
                        seeds.push_back(Seed(nppv[i].getNodeID(), nppv[i].getPosition(), i, i + 1));
                } else {
                        seeds.back().readEnd = i + 1;
                }
        }
}



void ReadCorrection::recSearch(NodeID curr, string& read, vector<NodePosPair>& npp,
                                  size_t currReadPos, size_t& counter,
                                  int currScore, int& bestScore, size_t& seedLast, bool &fullyCorrected,string currentNode)
{
        const SSNode node = dbg.getSSNode(curr);
        
        counter++;
        if (counter > (size_t)settings.getReadCorrDFSNodeLimit())
        {
                fullyCorrected = false;
                return ;
        }
        vector<DFSNode> dfsNode;

        size_t readCharLeft = getMarginalLength(read) - currReadPos;

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                NodeID nextID = it->getNodeID();
                const SSNode nextNode = dbg.getSSNode(nextID);

                size_t OLSize = min(nextNode.getMarginalLength(), readCharLeft);
                size_t nextReadPos = currReadPos + OLSize;

                string nodeOL = nextNode.substr(Kmer::getK()-1, OLSize);
                string readOL = read.substr(currReadPos + Kmer::getK() - 1, OLSize);
                string discoveredNode = currentNode + nodeOL;

                /*cout << nextID <<endl;
                string discoveredRead = read.substr(currReadPos + readOL.size()+ Kmer::getK() - 1 - discoveredNode.size(), discoveredNode.size());
                int fullScore = alignment.align(discoveredNode, discoveredRead);
                alignment.printAlignment(discoveredNode, discoveredRead);
                cout <<fullScore <<endl;*/
                //int thisScore = align(readOL.c_str(), nodeOL.c_str());
                int thisScore = alignment.alignBanded(readOL, nodeOL);
                
                int nextScore = currScore + thisScore;
                float nextRelScore = (float)thisScore / (float)OLSize;

                dfsNode.push_back(DFSNode(nextID, nextReadPos, nextScore, nextRelScore, discoveredNode));
                // =====================
                /*string str = nextNode.substr(Kmer::getK()-1, OLSize);
                string readSubStr = read.substr(currReadPos + Kmer::getK() -1, OLSize);
                for (size_t i = 0; i < currReadPos + Kmer::getK() - 1; i++)
                        cout << " ";
                cout << str << " (Node: " << nextID << ", curr: " << nextScore << ", best: " << bestScore << ", rel. score.:" << nextRelScore << ")" << endl;*/
                // ======================
        }
        sort(dfsNode.begin(), dfsNode.end());
        
               for (auto& it : dfsNode) {
                NodeID nextID = it.nodeID;
                int nextScore = it.score;
                size_t nextReadPos = it.readPos;
                size_t OLSize = nextReadPos - currReadPos;
                string discoveredNode = it.discoveredNode;

                int maxAttainScore = nextScore + getMarginalLength(read) - nextReadPos;

                // save and, if necessary, update the best score
                int prevBestScore = bestScore;
                if ((nextScore > bestScore) /*&& (currReadPos + OLSize > seedLast)*/) {
                        bestScore = nextScore;
                        seedLast = currReadPos + OLSize;
                }

                // =====================
               /* SSNode nextNode = dbg.getSSNode(nextID);
                string str = nextNode.substr(Kmer::getK()-1, OLSize);
                string readSubStr = read.substr(currReadPos + Kmer::getK() -1, OLSize);
                for (size_t i = 0; i < currReadPos + Kmer::getK() - 1; i++)
                        cout << " ";
                cout << str << " (Node: " << nextID << ", curr: " << nextScore << ", best: " << bestScore << ", max: " << maxAttainScore << ")" << endl;*/
                // ======================

                // descend in a child node only if there is chance this will improve the best score
                if (maxAttainScore > bestScore)
                        recSearch(nextID, read, npp, nextReadPos, counter, nextScore, bestScore, seedLast,fullyCorrected, discoveredNode);

                // if the best score has been updated in this branch save the npp
                if (bestScore > prevBestScore)
                        for (size_t i = 0; i < OLSize; i++)
                                npp[currReadPos+i] = NodePosPair(nextID, i);
        }
}



SearchRes ReadCorrection::extendSeedSlow (string& read, vector<NodePosPair>& npp,
                                   size_t& seedFirst, size_t& seedLast, bool& branch){
        branch = false;
          // read is fully aligned
        if (seedLast >= getMarginalLength(read))
                return OPTIMAL;
        // remove at most k nucleotides from the seed
        for (size_t i = 0; i < Kmer::getK(); i++) {
                if (seedLast == seedFirst + 1)
                        break;

                npp[seedLast - 1] = NodePosPair(0, 0);
                seedLast--;
        }

        // try to extend to the seed to the right within the current node
        const SSNode node = dbg.getSSNode(npp[seedLast-1].getNodeID());
        while ((seedLast < getMarginalLength(read)) &&
               (npp[seedLast - 1].getPosition() + 1 < node.getMarginalLength())) {
                npp[seedLast] = NodePosPair(node.getNodeID(), npp[seedLast-1].getPosition() + 1);
                seedLast++;
        }

        if (seedLast >= getMarginalLength(read))
                return OPTIMAL;
        
        branch = true;
        size_t readStart = seedLast + Kmer::getK() - 1;
        NodeID nodeID = npp[seedLast-1].getNodeID();
        size_t nodeStart = npp[seedLast-1].getPosition() + Kmer::getK();

        size_t curr = seedFirst;
        vector<NodeID> nodeChain;
        while (curr < seedLast) {
                NodeID currID = npp[curr].getNodeID();
                nodeChain.push_back(currID);
                SSNode curNode = dbg.getSSNode(currID);
                size_t strLen = min(curNode.getLength() - npp[curr].getPosition(), read.size() - curr);
                curr += (strLen - Kmer::getK() + 1);
        }
        

        GraphAlignment ga;
        ga.updatePath(NodeAlignment(0, 0, 0, 0, 0, -100000, 0), 1);     // make the worst possible start node
         
        SearchRes res = graphaln.DFSAln_guided(read.substr(readStart), NodePosPair(nodeID, nodeStart), ga, markov_filter, branchAndBound, nodeChain);

                
        if (res == EXHAUSTED || ga.getScore() == -100000)
                return EXHAUSTED;

        // transform the results to k-mer space
        vector<NodePosPair> alignedKmerNPP;
        ga.getKmerNPP(alignedKmerNPP);
        size_t elementsToCopy = min(npp.size() - seedLast, alignedKmerNPP.size());
        copy(alignedKmerNPP.begin(), alignedKmerNPP.begin() + elementsToCopy, npp.begin() + seedLast);

        seedLast += elementsToCopy;

        return res;
        
                                           
}
bool ReadCorrection::extendSeedFast(string& read, vector<NodePosPair>& npp,
                                   size_t& seedFirst, size_t& seedLast, bool& branch)
{
        // read is fully aligned
        if (seedLast >= getMarginalLength(read))
                return true;
        // remove at most k nucleotides from the seed
        for (size_t i = 0; i < Kmer::getK(); i++) {
                if (seedLast == seedFirst + 1)
                        break;

                npp[seedLast - 1] = NodePosPair(0, 0);
                seedLast--;
        }

        // try to extend to the seed to the right within the current node
        const SSNode node = dbg.getSSNode(npp[seedLast-1].getNodeID());
        while ((seedLast < getMarginalLength(read)) &&
               (npp[seedLast - 1].getPosition() + 1 < node.getMarginalLength())) {
                npp[seedLast] = NodePosPair(node.getNodeID(), npp[seedLast-1].getPosition() + 1);
                seedLast++;
        }

        // read is fully aligned
        if (seedLast >= getMarginalLength(read))
                return true;
        
        branch = true;
        // try to find the best right path
        size_t counter = 0; int bestScore = -(getMarginalLength(read) - seedLast);
        bool fullyCorrected = true;
        string discoveredNode = read.substr(seedFirst,seedLast-seedFirst+Kmer::getK()-1);
        if (seedLast < getMarginalLength(read))
                 recSearch(node.getNodeID(), read, npp, seedLast, counter, 0, bestScore, seedLast, fullyCorrected,discoveredNode);
        return fullyCorrected;
}

void ReadCorrection::applyReadCorrection(string& read,
                                            const vector<NodePosPair>& npp,
                                            size_t seedFirst, size_t seedLast, vector<NodeID>& nodeChain, pair<size_t, size_t> &nodePos)
{

        // correct the seed
        size_t curr = seedFirst;
        size_t nodeFirstPos = npp[curr].getPosition();
        while (curr < seedLast) {
                NodeID currID = npp[curr].getNodeID();
                 nodeChain.push_back(currID);
                
                SSNode node = dbg.getSSNode(npp[curr].getNodeID());
                size_t strLen = min(node.getLength() - npp[curr].getPosition(), read.size() - curr);
                string str = node.substr(npp[curr].getPosition(), strLen);

                read.replace(curr, strLen, str);
                curr += (strLen - Kmer::getK() + 1);
        }
        size_t nodeLastPos  = npp[seedLast-1].getPosition() + Kmer::getK();
        nodeLastPos = dbg.getSSNode( npp[seedLast-1].getNodeID()).getLength() -  nodeLastPos ;
        nodePos.first = nodeFirstPos;
        nodePos.second = nodeLastPos;
}

SearchRes ReadCorrection::correctRead(string& read, vector<NodePosPair> npp,
                                    size_t& first, size_t& last, vector<NodeID>& nodeChain, bool & branch, pair<size_t, size_t> &nodePos)
{
    
        // extend to the right
        bool branchRight = false;
        SearchRes resRight = extendSeedSlow (read, npp, first, last, branchRight);
        //bool right = extendSeedFast (read, npp, first, last, branchRight);
        // reverse complement all data
        Nucleotide::revCompl(read);
        first = getMarginalLength(read) - first;
        last = getMarginalLength(read) - last;
        swap(first, last);
        revCompl(npp);

        // extend to the right (which used to be left)
        bool branchLeft = false;
        SearchRes resLeft = extendSeedSlow (read, npp, first, last, branchLeft);
        //bool left = extendSeedFast (read, npp, first, last, branchLeft);
        
        Nucleotide::revCompl(read);
        first = getMarginalLength(read) - first;
        last = getMarginalLength(read) - last;
        swap(first, last);
        revCompl(npp);

        branch = branchLeft || branchRight ;
        
        // reverse complement all data
        string original = read;

        applyReadCorrection(read, npp, first, last, nodeChain,  nodePos);
        if (resLeft ==OPTIMAL && resRight== OPTIMAL)
                return OPTIMAL;
        else{
                if (resLeft ==EXHAUSTED || resRight== EXHAUSTED)
                        return EXHAUSTED;
                else
                        return MULTIPLE;
        }
}

void ReadCorrection::findSeedKmer(const std::string& read,
                                     vector<Seed>& mergedSeeds)
{
        vector<NodePosPair> nppv(read.length() + 1 - Kmer::getK());

        // find the node position pairs using the kmer lookup table
        findNPPFast(read, nppv);

        // transform consistent npps to seeds
        vector<Seed> seeds;
        extractSeeds(nppv, seeds);

        // sort seeds according to nodeID
        sort(seeds.begin(), seeds.end());

        // merge seeds
        Seed::mergeSeeds(seeds, mergedSeeds);

        /*#ifdef DEBUG
        // ----------- OUTPUT ------------
        cout << "NPP after kmer lookup search" << endl;
        cout << read << endl;
        for (size_t i = 0; i < Kmer::getK() - 1; i++)
                cout << " ";
        for (size_t i = 0; i < nppv.size(); i++) {
                if (nppv[i].isValid())
                        //cout << i << ":" << npp[i].getNodeID() << ":" << npp[i].getOffset() << " ";
                        cout << "|";
                else
                        cout << "*";
        }
        cout << endl;

        for (size_t i = 0; i < mergedSeeds.size(); i++) {
                for (size_t j = 0; j < mergedSeeds[i].readFirst + Kmer::getK() - 1; j++)
                        cout << " ";
                cout << "[";
                for (size_t j = mergedSeeds[i].readFirst + 1; j < mergedSeeds[i].readEnd; j++)
                        cout << "-";
                cout << "[";
                cout << " (" << mergedSeeds[i].readFirst << " to " << mergedSeeds[i].readEnd << ")\t";
                for (size_t j = 0; j < mergedSeeds[i].nodeID.size(); j++)
                        cout << mergedSeeds[i].nodeID[j] << "\t";
                cout << endl;
        }
#endif*/
}

bool sortByLength(const Seed& a, const Seed& b) {
        return ((a.readEnd-a.readFirst) > (b.readEnd-b.readFirst));
}

void ReadCorrection::findSeedMEM(const string& read,
                                    vector<Seed>& mergedSeeds)
{
        vector<match_t> matches;

        int memSize = Kmer::getK() - 1;
        while (matches.size() < 100 && memSize>15) {
                matches.clear();
                //cout << "Find MEM: " << memSize << endl;
                sa.findMEM(0l, read, matches, memSize, false);
                //cout << "Number of matches for size " << memSize << ": " << matches.size() << endl;
                memSize--;
        }


        vector<Seed> seeds;
        seeds.reserve(matches.size());
        for (match_t it: matches) {

                vector<long>::const_iterator e = upper_bound(startpos.begin(), startpos.end(), it.ref);
                e--;
                NodeID nodeID = distance(startpos.begin(), e) + 1;
                size_t nodeFirst = it.ref - *e;

                SSNode node = dbg.getSSNode(nodeID);
                if (nodeFirst > node.getLength()) {
                        nodeID = -nodeID;
                        nodeFirst = nodeFirst - node.getLength() - 1;
                }

                // don't select MEMs that are closer than k nucleotides to the
                // right edge of a node as this MEM is also in its right neighbor

                if (nodeFirst >= node.getMarginalLength())
                        continue;

                size_t readFirst = it.query;

                // don't select MEMs that are closer than k nucleotides to the
                // right edge of the read: NECESSARY ?!?

                if (readFirst >= getMarginalLength(read))
                        continue;

                size_t readEnd = readFirst + it.len; //(it.len > Kmer::getK() ? it.len +1 - Kmer::getK() : 1);

                seeds.push_back(Seed(nodeID, nodeFirst, readFirst, readEnd));
        }

        // sort seeds according to nodeID
        sort(seeds.begin(), seeds.end());
        // merge seeds
        Seed::mergeSeeds(seeds, mergedSeeds);
        // sort seeds according to length
        sort(mergedSeeds.begin(), mergedSeeds.end(), sortByLength);

        // retain only 10 largest seeds
        if (mergedSeeds.size() > 10)
                mergedSeeds.resize(10);

        for (size_t i = 0; i < mergedSeeds.size(); i++)
                mergedSeeds[i].readEnd -= min(Kmer::getK() - 1, mergedSeeds[i].readEnd - mergedSeeds[i].readFirst - 1);


        // ----------- OUTPUT ------------
        /*cout << "NPP after kmer lookup search" << endl;
        cout << read << endl;

        for (size_t i = 0; i < mergedSeeds.size(); i++) {
                for (size_t j = 0; j < mergedSeeds[i].readFirst; j++)
                        cout << " ";
                cout << "[";
                for (size_t j = mergedSeeds[i].readFirst + 1; j < mergedSeeds[i].readEnd; j++)
                        cout << "-";
                cout << "[";
                cout << " (" << mergedSeeds[i].readFirst << " to " << mergedSeeds[i].readEnd << ")\t";
                for (size_t j = 0; j < mergedSeeds[i].nodeID.size(); j++)
                        cout << mergedSeeds[i].nodeID[j] << "\t";
                cout << endl;
        }*/
        // ----------- OUTPUT ------------
}

SearchRes ReadCorrection::correctRead(const string& read,
                                 string& bestCorrectedRead,
                                 const vector<Seed>& seeds,int &bestScore, vector<NodeID>& bestNodeChain, bool &branch,  pair<size_t, size_t> &nodePos )
{
        SearchRes CorrectionStatus = EXHAUSTED;
        pair<size_t, size_t> nodePosLocal;
        for (auto& it : seeds) {
                size_t first = it.readFirst;
                size_t last = it.readEnd;
                // create a npp vector
                vector<NodePosPair> npp(read.length() + 1 - Kmer::getK());
                it.createNodePosition(dbg, npp);
                string correctedRead = read;
                vector<NodeID> nodeChain;
                bool branchLocal = false;
                SearchRes res = correctRead(correctedRead, npp, first, last, nodeChain, branchLocal, nodePosLocal);
                int score = alignment.alignBanded(read, correctedRead);
                size_t correctedLength = Kmer::getK() - 1 + last - first;
                size_t uncorrectedLength = read.size() - correctedLength;
                score = score - uncorrectedLength;
                if (score > bestScore ) {
                        CorrectionStatus = res;
                        branch = branchLocal;
                        bestScore = score ;
                        bestCorrectedRead = correctedRead;
                        bestNodeChain = nodeChain;
                        nodePos = nodePosLocal;
                        
                }
        }
        return (CorrectionStatus );
}



SearchRes ReadCorrection::correctRead(const string &read, string &bestCorrectedRead , bool& correctedByMEM, int minSim, int & bestScore, vector<NodeID> &bestNodeChain, bool &branch , pair<size_t, size_t>& nodePos){
        vector<Seed> seeds;
        findSeedKmer(read, seeds);
        bool branchLocal = false;
        SearchRes CorrectionStatus = correctRead(read, bestCorrectedRead, seeds, bestScore, bestNodeChain, branchLocal, nodePos);
        branch = branchLocal;
        if (seeds.empty() || (bestScore <= minSim && CorrectionStatus != EXHAUSTED) ) {
                findSeedMEM(read, seeds);
                bestNodeChain.clear();
                int bestScoreMem = - read.size();
                SearchRes CorrectionStatusMem = correctRead(read, bestCorrectedRead, seeds, bestScoreMem, bestNodeChain, branchLocal, nodePos);
                if (bestScoreMem > minSim && CorrectionStatusMem != EXHAUSTED){
                        branch = branchLocal;
                        correctedByMEM = true;
                        bestScore = bestScoreMem ;
                        CorrectionStatus = CorrectionStatusMem;
                }
        }
        return CorrectionStatus;
}

void ReadCorrection::correctReadFirstAttempt(ReadRecord& record,
                                    AlignmentMetrics& metrics)
{
        
        bool correctedByMEM = false, readCorrected = false ,  branch = false;
        vector<NodeID> bestNodeChain;
        string& read = record.getRead(), bestCorrectedRead = "";
        int minSim = ((int)read.size() / 2), bestScore = -read.size();
        size_t numSubstitutions = 0;
        
        // the offset positions in the first node that read aligns (normal and RC)
        pair<size_t, size_t> nodePos;
        string alignmentInfo = record.preRead.substr(0,record.preRead.length()-1);
        
        // if the read is too short, get out
        if (read.length() < Kmer::getK()){
                record.alignmentInfo  = alignmentInfo+ "\tA(0)A\tS(0)S\tB(0)B\tP(0,0)P\tC(0)C";
                record.preRead.insert(record.preRead.length()-1,"\tA(0)A\tS(0)S\tB(0)B\tP(0,0)P\tC(0)C");
                return;
        }
        

        SearchRes correctionStatus = correctRead(read, bestCorrectedRead, correctedByMEM, minSim, bestScore , bestNodeChain, branch, nodePos);

        
        //Don't use the alignment info for the MM if it is exhusted or the similarity score is low
        if (correctionStatus == EXHAUSTED || bestScore < minSim){
                bestScore = -read.size();
                bestNodeChain.clear();
        }
        // if the read is needed to be aligned accros the nodes with dfs algorithm, 
        // then branch is set to 1 and kept as insertedStr in the header of read for the next round of alignment
        
        string insertedStr = makeAlignmentInfo(branch, bestScore, bestNodeChain, nodePos);
    
        record.nodeChain = bestNodeChain;
        
        // In the first attemp we don't wirte back reads that are aligned accros multiple nodes with dfs algorithm
        // This can change to CorrectionStatus != EXHAUSTED if you want to consider also MULTIPLE
        
        if ( !branch  && bestScore > minSim && correctionStatus == OPTIMAL && dbg.validateChain(bestNodeChain)) {
                insertedStr =  "\tA(1)A\t" + insertedStr;
                read = bestCorrectedRead;
                readCorrected = true;
                numSubstitutions = (read.length() - bestScore)/2;
        }
        else{
                if (!branch)
                        insertedStr =  "\tA(0)A\t" + insertedStr;
                else
                        insertedStr =  "\tA(1)A\t" + insertedStr;
        }
        record.preRead.insert(record.preRead.length()-1,insertedStr);
        record.alignmentInfo =  alignmentInfo + insertedStr ;
        metrics.addObservation(readCorrected, correctedByMEM, numSubstitutions, read.length());
}

void ReadCorrection::correctReadSecondAttempt(ReadRecord& record,
                                    AlignmentMetrics& metrics)
{
        
        bool correctedByMEM = false, readCorrected = false, branch = false;;
        vector<NodeID> bestNodeChain;
        string& read = record.getRead();
        int minSim = ((int)read.size() / 2);
        string bestCorrectedRead;
        int bestScore = -read.size();
        pair<size_t, size_t> nodePos;
        SearchRes correctionStatus = correctRead(read, bestCorrectedRead, correctedByMEM, minSim, bestScore , bestNodeChain, branch, nodePos);
        
        size_t last = record.preRead.find("\tA(");
        string alignmentInfo = record.preRead.substr (0,last) ;
        
        
        string insertedStr = makeAlignmentInfo(branch, bestScore, bestNodeChain, nodePos );
      
        size_t numSubstitutions = 0;
        if (bestScore > minSim && correctionStatus == OPTIMAL && dbg.validateChain(bestNodeChain)) {
                read = bestCorrectedRead;
                readCorrected = true;
                numSubstitutions = (read.length() - bestScore)/2;

                insertedStr = "\tA(1)A\t"+ insertedStr ;
        }
        else
                insertedStr = "\tA(0)A\t" +insertedStr ;
        
     
        record.alignmentInfo =  alignmentInfo + insertedStr;
        record.preRead = record.alignmentInfo +"\n";
        metrics.addObservation(readCorrected, correctedByMEM, numSubstitutions, read.length());
}
string ReadCorrection::makeAlignmentInfo(bool branch, int bestScore, vector<NodeID> bestNodeChain, pair<size_t, size_t> nodePos){
        string branchStr = branch ? "1" : "0";
        
        //save the original read beside the corrected one, later in refinement we need to correct it again.
        string insertedStr  = "S("+to_string( bestScore)+")S\tB("+ branchStr +")B";
        insertedStr = insertedStr+ "\tP(" +to_string( nodePos.first) + "," +to_string( nodePos.second) + ")P";

        string nodeCahinStr = "\tC(";
        if (bestNodeChain.empty())
                nodeCahinStr = nodeCahinStr+ "0" ;
        else{
                for (size_t j = 0; j < bestNodeChain.size() - 1; j++)
                        nodeCahinStr = nodeCahinStr +to_string( bestNodeChain[j]) + " ";
                nodeCahinStr = nodeCahinStr + to_string(  bestNodeChain.back()) ;
        }
        insertedStr = insertedStr + nodeCahinStr+")C" ;
        return insertedStr;
}

void ReadCorrection::correctChunk(vector<ReadRecord>& readChunk,
                                     AlignmentMetrics& metrics)
{
        for (auto& it : readChunk)
                correctReadFirstAttempt(it, metrics);
}


void ReadCorrection::refineChunk(vector<ReadRecord>& readChunk,
                                 AlignmentMetrics& metrics)
{

        for (auto& it : readChunk){
                size_t end = it.preRead.find(")A\t");
                string alignmentInfo = it.preRead.substr (0,end+2) ;
                if (it.read.length() < Kmer::getK()){
                        it.alignmentInfo = it.preRead.substr(0,it.preRead.length()-1);
                        it.preRead = alignmentInfo + "\n";
                        continue;
                }
                size_t first = it.preRead.find("B(") + 2;
                size_t last = it.preRead.find(")B");
                bool branch = it.preRead.substr (first,last-first) == "0" ? false : true;
                if (branch)
                        correctReadSecondAttempt(it, metrics);
                else
                        it.alignmentInfo = it.preRead.substr(0,it.preRead.length()-1);
                end = it.preRead.find(")A\t");
                alignmentInfo = it.preRead.substr (0,end+2) ;
                it.preRead = alignmentInfo + "\n";
                

        }
        
}
void ReadCorrectionHandler::workerThread(size_t myID, LibraryContainer& libraries,
                                         AlignmentMetrics& metrics)
{
        ReadCorrection readCorrection(dbg, settings, *sa, startpos, mch, false, branchAndBound);

        // local storage of reads
        vector<ReadRecord> myReadBuf;

        // performance counters per thread
        AlignmentMetrics threadMetrics;

        while (true) {
                size_t blockID, recordID;
                bool result = libraries.getRecordChunk(myReadBuf, blockID, recordID);

                readCorrection.correctChunk(myReadBuf, threadMetrics);

                if (result)
                        libraries.commitRecordChunk(myReadBuf, blockID, recordID);
                else
                        break;
        }
        
        // update the global metrics with the thread info (thread-safe)
        metrics.addMetrics(threadMetrics);
}

void ReadCorrectionHandler::initEssaMEM()
{
        size_t length = 0;
        startpos.clear();
        startpos.reserve(dbg.getNumNodes());
        for (NodeID nodeID = 1; nodeID <= dbg.getNumNodes(); nodeID++) {
                SSNode node = dbg.getSSNode(nodeID);
                assert(node.isValid());
                startpos.push_back(length);
                length += node.getLength() * 2 + 2;
        }
        reference.clear();
        reference.reserve(length);
        for (NodeID nodeID = 1; nodeID < dbg.getNumNodes(); nodeID++) {
                SSNode node = dbg.getSSNode(nodeID);
                if (!node.isValid())
                        continue;

                string thisSequence = node.getSequence();
                reference.append(thisSequence);
                reference.append(">");

                Nucleotide::revCompl(thisSequence);
                reference.append(thisSequence);
                reference.append(">");
        }

        std::vector<std::string> refdescr;
        refdescr.push_back("");

        bool printSubstring = false;
        bool printRevCompForw = false;

        sa = new sparseSA(reference,                            // reference
                          refdescr,
                          startpos,                             // start index for each string
                          false,                                // 4 column format or not
                          settings.getEssaMEMSparsenessFactor(),// ESSA sparseness factor
                          true,                                 // suffixlinks
                          true,                                 // child arrays,
                          true,                                 // kmertable
                          1,             
                          // skip parameter
                          10,                                   // kmer size
                          printSubstring,
                          printRevCompForw,
                          false                                 // nucleotides only
                          );
        sa->construct();
}
void ReadCorrectionHandler::saveMetricsInFile(AlignmentMetrics &metrics){
        ofstream myfile;
        myfile.open (settings.getTempDirectory()+"sta.txt");
        myfile <<"NumOfReads:"       << metrics.getNumReads()         <<endl;
        myfile <<"AvgReadLen:"       << metrics.getAvgReadLen()       <<endl;
        myfile <<"MaxReadLen:"       << metrics.getMaxReadLen()       <<endl;
        myfile <<"NumSubstitutions:" << metrics.getNumSubstitutions() <<endl;
        myfile <<"NumCorrReads:"     << metrics.getNumCorrReads()     <<endl;
        myfile <<"NumCorrByMEM:"     << metrics.getNumCorrByMEM()     <<endl;
        myfile.close();
}
vector<string> splitString(const string str,const char delimiter)
{

        vector<string> parsed;
        size_t pos = 0;
        while (true) {
                size_t colon = str.find(delimiter, pos);
                if (colon == ::std::string::npos) {
                        parsed.push_back(str.substr(pos));
                        break;
                } else {
                        parsed.push_back(str.substr(pos, colon - pos));
                        pos = colon + 1;
                }
        }
        return parsed;
}
void  ReadCorrectionHandler::loadMetrics(size_t& tootalNumOfReads, size_t &numOfSub, size_t &numCorrectedByMEM ,size_t & numCorrReads )
{
        std::ifstream ifs(settings.getTempDirectory()+"sta.txt");
        std::set<string> ids;
        while(ifs.is_open() && !ifs.eof())
        {
                string line = "";
                std::getline(ifs, line);
                vector<string> items = splitString(line,':');
                if(items[0] ==  "NumSubstitutions"){
                         numOfSub= atoi (items[1].c_str());
                }
                if(items[0] ==  "NumCorrReads"){
                        numCorrReads = atoi (items[1].c_str());
                }
                if(items[0] ==  "NumCorrByMEM"){
                        numCorrectedByMEM = atoi (items[1].c_str());
                }
                if (items[0] == "NumOfReads"){
                        tootalNumOfReads = atoi (items[1].c_str());
                } 
        }
        
        ifs.close();
}


void ReadCorrectionHandler::doErrorCorrection(LibraryContainer& libraries)
{
        const unsigned int& numThreads = settings.getNumThreads();;//
        cout << "Number of threads: " << numThreads << endl;

        libraries.startIOThreads(settings.getThreadWorkSize(),
                                 10 * settings.getThreadWorkSize() * numThreads,
                                 true);

        // start worker threads
        AlignmentMetrics metrics;
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&ReadCorrectionHandler::workerThread,
                                          this, i, ref(libraries), ref(metrics));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
        
        libraries.joinIOThreads();
        //metrics.printStatistics();
        saveMetricsInFile(metrics);
}

ReadCorrectionHandler::ReadCorrectionHandler(DBGraph& g, const Settings& s,bool markov_filter_) :
        dbg(g), settings(s), sa(NULL),mch(dbg,settings), markov_filter(markov_filter_), branchAndBound(settings.getBranchAndBound())
{
        Util::startChrono();
        cout << "Creating kmer lookup table... "; cout.flush();
        dbg.buildKmerNPPTable();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        Util::startChrono();
        cout << "Building suffix array (sparseness factor: "
             << settings.getEssaMEMSparsenessFactor() << ")..."; cout.flush();
        initEssaMEM();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
}

void ReadCorrectionHandler::workerThreadRefine(size_t myID, LibraryContainer& libraries,
                                         AlignmentMetrics& metrics)
{
        
        ReadCorrection readCorrection(dbg, settings, *sa, startpos, mch, markov_filter, branchAndBound);

        // local storage of reads
        vector<ReadRecord> myReadBuf;

        // performance counters per thread
        AlignmentMetrics threadMetrics;

        while (true) {
                size_t blockID, recordID;
                bool result = libraries.getRecordChunk(myReadBuf, blockID, recordID);
                readCorrection.refineChunk(myReadBuf, threadMetrics);
                
                
                if (result)
                        libraries.commitRecordChunk(myReadBuf, blockID, recordID);
                else
                        break;
        }
        
        // update the global metrics with the thread info (thread-safe)
        metrics.addMetrics(threadMetrics);
}
void ReadCorrectionHandler::doReadRefinement(LibraryContainer& libraries)
{
        // first load the Markov Table into mamory, if it dosn't exist populate table 
        if (markov_filter){
                Util::startChrono();
                bool fileLoad = false;
                if (mch.markovTableFileExist()){
                        cout << "Markov table appear to be present on the disk, loading from file..." << endl ;
                        if (!mch.loadFromFile())
                                cout << "Loading from file was unsuccessful" <<endl;
                        else
                                fileLoad = true;
                }
                if (!fileLoad){
                        cout << "Creating Markov Chain Transition table... "; cout.flush();
                        mch.extractTransitionMatrixFromFile(libraries);
                        mch.writeInFile();
                }
                cout << "done (" << Util::stopChronoStr() << ")" << endl;
        }
        // the result of the first attempt error correction is saved in temp/*.final.*
        LibraryContainer libContErroneous;
        for (size_t i = 0; i < libraries.getSize(); i++) {
                const ReadLibrary &input = libraries.getInput(i);
                string outputFileName = input.getOutputFileName();
                std::ostringstream oss;
                oss << settings.getTempDirectory() <<input.getBaseFilename() << + ".final." << input.getFileType() ;
                std::string tempFileStr = oss.str();
                libContErroneous.insert(ReadLibrary(outputFileName,tempFileStr,settings.getTempDirectory()));
                
        }
        const unsigned int& numThreads = settings.getNumThreads();;//
        cout << "Number of threads: " << numThreads << endl;

        libContErroneous.startIOThreads(settings.getThreadWorkSize(),
                                 10 * settings.getThreadWorkSize() * numThreads,
                                 true);

        // start worker threads
        AlignmentMetrics metrics;
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&ReadCorrectionHandler::workerThreadRefine,
                                          this, i, ref(libContErroneous), ref(metrics));
        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
     
        libContErroneous.joinIOThreads();
        size_t tootalNumOfReads = 0, numOfSub = 0, numCorrectedByMEM =0, numCorrReads =0;
        
        loadMetrics(tootalNumOfReads, numOfSub, numCorrectedByMEM , numCorrReads );
        metrics.updateMetrics(tootalNumOfReads, numOfSub, numCorrectedByMEM , numCorrReads);
        metrics.printStatistics();

        // keep a copy of file to extract the branching reads for further analises, this is temporary 
       /* for (size_t i = 0; i < libraries.getSize(); i++) {
                const ReadLibrary &input = libraries.getInput(i);
              
                std::ostringstream oss;
                oss << input.getBaseFilename() << + ".firstAttempt." << input.getFileType() ;
                std::string tempFileStr = oss.str();  
                string outputFileName = input.getOutputFileName();
                std::rename(  outputFileName.c_str(),  tempFileStr.c_str());
        }*/
        
        // rename the output file to the name that user already gave and remove temporary *.fnal.* file
        for (size_t i = 0; i < libraries.getSize(); i++) {
                const ReadLibrary &input = libraries.getInput(i);
                ostringstream inputFileName ;
                std::ostringstream oss;
                oss << settings.getTempDirectory() <<input.getBaseFilename() << + ".final." << input.getFileType() ;
                std::string tempFileStr = oss.str();  
                string outputFileName = input.getOutputFileName();
                std::rename( tempFileStr.c_str(), outputFileName.c_str());
        }
}
ReadCorrectionHandler::~ReadCorrectionHandler()
{
        delete sa;
        dbg.destroyKmerNPPTable();
}
