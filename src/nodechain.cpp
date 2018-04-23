/***************************************************************************
 *   Copyright (C) 2014 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
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

#include "nodechain.h"

#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <list>

using namespace std;

// ============================================================================
// OPERATOR <<
// ============================================================================

std::ostream &operator<<(std::ostream &out, const vector<NodeID> &path)
{
        if (!path.empty())
                out << path.front();

        for (size_t i = 1; i < path.size(); i++)
                out << " " << path[i];

        return out;
}

std::ostream &operator<<(std::ostream &out, const NodeChain &nc)
{
        for (auto& it : nc)
                out << it << " ";
        out << "(" << nc.count << ")";

        return out;
}

// ============================================================================
// SORT DECLARATIONS
// ============================================================================

bool sortChainByOcc(const NodeChain& left, const NodeChain& right)
{
        return left.getCount() > right.getCount();
}

bool sortChainByNodeID(const NodeChain& left, const NodeChain& right)
{
        for (size_t i = 0; i < left.size() && i < right.size(); i++)
                if (left[i] != right[i])
                        return left[i] < right[i];

        return left.size() > right.size();
}

bool sortByPosition(const NodeChainPos& left, const NodeChainPos& right)
{
        if (left.getChainID() != right.getChainID())
                return left.getChainID() < right.getChainID();
        return left.getChainPos() < right.getChainPos();
}

// ============================================================================
// NODE CHAIN CLASS
// ============================================================================

NodeChain NodeChain::getReverseComplement() const
{
        NodeChain copy = *this;

        reverse(copy.begin(), copy.end());
        for (size_t i = 0; i < copy.size(); i++)
                copy[i] = -copy[i];

        return copy;
}

bool NodeChain::operator<(const NodeChain& rhs) const
{
        if (size() != rhs.size())
                return size() < rhs.size();
        for (size_t i = 0; i < size(); i++)
                if ((*this)[i] != rhs[i])
                        return (*this)[i] < rhs[i];
        return false;
}

bool NodeChain::operator==(const NodeChain& rhs) const
{
        if (size() != rhs.size())
                return false;
        for (size_t i = 0; i < size(); i++)
                if ((*this)[i] != rhs[i])
                        return false;
        return true;
}

// ============================================================================
// NODE CHAIN CONTAINER CLASS
// ============================================================================

ostream &operator<<(ostream &out, const NodeChainContainer& ncc)
{
        for (const auto& it : ncc)
                out << it << endl;

        return out;
}

void NodeChainContainer::buildIndex()
{
        for (size_t i = 0; i < size(); i++) {
                const NodeChain& nc = (*this)[i];

                for (size_t j = 0; j < nc.size(); j++) {
                        NodeChainPos ncp(i, j);
                        insertIndex(nc[j], ncp);
                }
        }
}

NodeChainContainer::NodeChainContainer(const vector<NodeChain>& input)
{
        // create a temporary set
        std::set<NodeChain> tempSet;

        for (const auto& nc : input) {
                NodeChain repr = nc.getRepresentative();

                auto it = tempSet.insert(repr);
                bool inserted = it.second;

                // handle duplicates: sum the counts
                if (!inserted) {
                        size_t count = nc.getCount() + it.first->getCount();
                        const_cast<NodeChain&>(*it.first).setCount(count);
                }
        }

        // convert the set to a vector
        reserve(tempSet.size());
        for (const auto& nc : tempSet)
                push_back(nc);
        tempSet.clear();

        // sort chain by occurence
        sort(begin(), end(), sortChainByOcc);

        // create the index
        buildIndex();
}

void NodeChainContainer::findOcc(const NodeChain& pattern,
                                 vector<NodeChainPos>& occ) const
{
        // clear the output
        occ.clear();

        // get out early if the search pattern is empty
        if (pattern.empty())
                return;

        // use the index to quickly find good candidates
        NodeID seedID = pattern.front();
        for (auto it = index.find(abs(seedID)); it != index.end(); it++) {
                if (it->first != abs(seedID))
                        break;

                const size_t& id = it->second.getChainID();
                const size_t& pos = it->second.getChainPos();

                const NodeChain& nc = (*this)[id];
                assert(abs(nc[pos]) == abs(seedID));

                if (nc[pos] == seedID) {
                        if (pos + pattern.size() > nc.size())
                                goto BreakLoop;

                        for (size_t i = 0, j = pos; i < pattern.size(); i++, j++)
                                if (nc[j] != pattern[i])
                                        goto BreakLoop;

                        occ.push_back(NodeChainPos(id, pos, false));
                } else {
                        if (pos + 1 < pattern.size())
                                goto BreakLoop;

                        for (size_t i = 0, j = pos; i < pattern.size(); i++, j--)
                                if (nc[j] != -pattern[i])
                                        goto BreakLoop;

                        occ.push_back(NodeChainPos(id, pos, true));
                }

                BreakLoop:

                continue;
        }
}

size_t NodeChainContainer::getChainIdx(const NodeChain& nc) const
{
        NodeID leftID = nc.front();

        for (multimap<NodeID, NodeChainPos>::const_iterator it = index.find(abs(leftID)); it != index.end(); it++) {
                if (it->first != abs(leftID))
                        break;

                const NodeChain& target = at(it->second.getChainID());

                if (nc == target)
                        return it->second.getChainID();

        }

        // we should never reach this position
        assert(false);
}

void NodeChainContainer::addContainers(const vector<string>& filenames)
{
        
        //The format of this file changed , 
        //Now per line you have:  \tS(score)S\tB(branch)B\t\P(nodePos, nodePosRC)P\tC(node chains)C
        //score is the alignment score, branch is 1 if read needs to be aligned with the dfs aglorithm
        //The pair of position in P tag, indicates the position of alignment on the node for the read and reverse complement of that read
        //node chain is the chain of nodes that read maps on the graph separated by space like before
        //to get the previous content which was only the node chain you can use :
        //size_t first = line.find("C(") + 2;
        //size_t last = line.find(")C");
        //string chain = line.substr (first,last-first);
        // create a temporary set
        std::set<NodeChain> tempSet;

        // load all node chains
        for (const string& filename : filenames) {
                ifstream ifs(filename.c_str());
                cout << "Reading file " << filename << endl;
                size_t count = 0;

                while (true) {
                        // get a line from the input file
                        string line;
                        getline(ifs, line);
                        if (!ifs.good())
                                break;

                        if (count++ % OUTPUT_FREQUENCY == 0) {
                                cout << "Processing read " << count << "\r";
                                cout.flush();
                        }
                        
                        // this code added by mahdi
                        size_t first = line.find("C(") + 2;
                        size_t last = line.find(")C");
                        line = line.substr (first,last-first);
                        //end of modification
                 
                        // convert line to a node chain
                        vector<NodeID> inputVector;
                        istringstream iss(line);

                        NodeID nodeID;
                        while (iss >> nodeID)
                                inputVector.push_back(nodeID);

                        //if (inputVector.size() < 3)
                        //        continue;

                        NodeChain nc(inputVector);
                        NodeChain repr = nc.getRepresentative();

                        auto it = tempSet.insert(repr);
                        bool inserted = it.second;

                        // handle duplicates: sum the counts
                        if (!inserted) {
                                size_t count = nc.getCount() + it.first->getCount();
                                const_cast<NodeChain&>(*it.first).setCount(count);
                        }
                }

                cout << "Processing read " << count << endl;

                ifs.close();
        }

        cout << "Converting set to vector..." << endl;
        // convert the set to a vector
        clear();
        reserve(tempSet.size());
        for (const auto& nc : tempSet)
                push_back(nc);
        tempSet.clear();

        cout << "Building index..." << endl;
        // create the index
        buildIndex();

        cout << "Make paired-end reads index..." << endl;

        // read files again to figure out which nodechains are paired
        ifstream ifs[2];
        for (size_t i = 0; i < 2; i++)
                ifs[i].open(filenames[i].c_str());

        while (true) {
                // get a line from the input file
                string line1;
                getline(ifs[0], line1);
                if (!ifs[0].good())
                        break;
                // this code added by mahdi
                size_t first = line1.find("C(") + 2;
                size_t last = line1.find(")C");
                line1 = line1.substr (first,last-first);
                //end of modification
                
                string line2;
                getline(ifs[1], line2);
                if (!ifs[1].good())
                        break;

                // this code added by mahdi
                first = line2.find("C(") + 2;
                last = line2.find(")C");
                line2 = line2.substr (first,last-first);
                //end of modification
                
                // convert line to a node chain
                vector<NodeID> inputVector1;
                istringstream iss(line1);

                NodeID nodeID;
                while (iss >> nodeID)
                        inputVector1.push_back(nodeID);

                NodeChain nc1(inputVector1);
                NodeChain repr1 = nc1.getRepresentative();

                vector<NodeID> inputVector2;
                istringstream iss2(line2);

                while (iss2 >> nodeID)
                        inputVector2.push_back(nodeID);

                NodeChain nc2(inputVector2);
                NodeChain repr2 = nc2.getRepresentative();

                int64_t idx1 = getChainIdx(repr1);
                int64_t idx2 = getChainIdx(repr2);

                if (nc1 != repr1)
                        idx1 = -idx1;
                if (nc2 == repr2)               // TODO: Change innies to outies : auto-detect
                        idx2 = -idx2;

                PER[pair<size_t, size_t>(idx1, idx2)]++;
        }

        for (size_t i = 0; i < 2; i++)
                ifs[i].close();
}

NodeChain NodeChainContainer::getNodeChainSection(NodeID nodeID,
                                                  size_t id, size_t pos) const
{
        const NodeChain& nc = (*this)[id];
        bool forward = (nc[pos] == nodeID);

        NodeChain result;
        result.setCount(nc.getCount());

        if (forward) {
                result.push_back(nc[pos]);
                for (pos++; pos < nc.size(); pos++) {
                        result.push_back(nc[pos]);
                        if (nc[pos] == nodeID)  // break at loops
                                break;
                }
        } else {
                result.push_back(-nc[pos]);
                for (ssize_t cnt = pos - 1; cnt >= 0; cnt--) {
                        result.push_back(-nc[cnt]);
                        if (-nc[cnt] == nodeID) // break at loops
                                break;
                }
        }

        return result;
}

vector<NodeChain> NodeChainContainer::getNodeChainTree(NodeID leftID, NodeID rightID) const
{
        vector<NodeChain> result;

        for (multimap<NodeID, NodeChainPos>::const_iterator it = index.find(abs(leftID)); it != index.end(); it++) {
                if (it->first != abs(leftID))
                        break;

                size_t id = it->second.getChainID();
                size_t pos = it->second.getChainPos();

                const NodeChain& nc = (*this)[id];
                assert(abs(nc[pos]) == abs(leftID));

                bool forward = (nc[pos] == leftID);

                if (forward) {
                        if ((pos+1 >= nc.size()) || (nc[pos+1] != rightID))
                                continue;
                } else {
                        if ((pos == 0) || (nc[pos-1] != -rightID))
                                continue;
                }

                NodeChain nct = getNodeChainSection(leftID, id, pos);
                if (nct.size() >= 3)
                        result.push_back(nct);
        }

        // sort the paths according to content and length
        sort(result.begin(), result.end(), sortChainByNodeID);

        return result;
}

size_t NodeChainContainer::haveConsensus(vector<pair<NodeID, size_t> >& consensus) const
{
        if (consensus.size() == 1)
                return 0;

        return consensus.size();
}

vector<NodeChain> NodeChainContainer::getNonBranchingPath(NodeID leftID, NodeID rightID) const
{
        vector<NodeChain> result;

        // create a node chain tree
        vector<NodeChain> nodeChainTree = getNodeChainTree(leftID, rightID);

        if (nodeChainTree.empty())
                return result;

        /*for (auto it : nodeChainTree)
                cout << it << endl;
        cout << " === " << endl;*/

        //exit(0);

        // find all indices in the node chain tree
        set<size_t> indices;
        for (size_t i = 0; i < nodeChainTree.size(); i++)
                indices.insert(i);

        NodeChain consensus;

        for (size_t pos = 0; !indices.empty(); pos++) {
                vector<pair<NodeID, size_t> > nodeCons;

                set<size_t> toDelete;
                for (size_t id : indices) {
                        if (pos >= nodeChainTree[id].size()) {
                                toDelete.insert(id);
                                continue;
                        }

                        // can we simply update the counter?
                        const NodeID& thisID = nodeChainTree[id][pos];
                        const NodeID& thisCount = nodeChainTree[id].getCount();

                        if (!nodeCons.empty())
                                if (nodeCons.back().first == thisID) {
                                        nodeCons.back().second += thisCount;
                                        continue;
                                }

                        // if not, add a new pair
                        nodeCons.push_back(pair<NodeID, size_t>(thisID, thisCount));
                }

                size_t consensusID = haveConsensus(nodeCons);
                if (consensusID == nodeCons.size())
                        break;

                // delete indices for which the node chain is too short
                for (size_t id : toDelete)
                        indices.erase(id);

                consensus.push_back(nodeCons[consensusID].first);
                consensus.setCount(nodeCons[consensusID].second);

                if (pos >= 2)
                        result.push_back(consensus);
        }

        return result;
}

bool NodeChainContainer::isNonBranchingPath(const NodeChain& nbPath) const
{
        if (nbPath.empty())
                return false;

        // create a node chain tree
        vector<NodeChain> nodeChainTree = getNodeChainTree(nbPath[0], nbPath[1]);

        // find all indices in the node chain tree
        set<size_t> indices;
        for (size_t i = 0; i < nodeChainTree.size(); i++)
                indices.insert(i);

        for (size_t pos = 0; pos < nbPath.size(); pos++) {

                // if reduction is longer than longest chainTree...
                if (indices.empty())
                        return false;

                // check the contents of the chainTrees
                set<size_t> toDelete;
                for (size_t id : indices) {
                        if (pos >= nodeChainTree[id].size()) {
                                toDelete.insert(id);
                                continue;
                        }

                        if (nodeChainTree[id][pos] != nbPath[pos])
                                return false;
                }

                // delete indices for which the node chain is too short
                for (size_t id : toDelete)
                        indices.erase(id);
        }

        return true;
}

void NodeChainContainer::insertIndex(NodeID nodeID, NodeChainPos pos)
{
        index.insert(pair<NodeID, NodeChainPos>(abs(nodeID), pos));
}

void NodeChainContainer::deleteIndex(NodeID nodeID, NodeChainPos oldPos)
{
        for (multimap<NodeID, NodeChainPos>::iterator it = index.find(abs(nodeID)); it != index.end(); it++) {
                if (it->first != abs(nodeID))
                        return;
                if (oldPos.getChainID() != it->second.getChainID())
                        continue;
                if (oldPos.getChainPos() != it->second.getChainPos())
                        continue;

                index.erase(it);
                return;
        }
}

void NodeChainContainer::updateIndex(NodeID nodeID, NodeChainPos oldPos, NodeChainPos newPos)
{
        for (multimap<NodeID, NodeChainPos>::iterator it = index.find(abs(nodeID)); it != index.end(); it++) {
                if (it->first != abs(nodeID))
                        break;
                if (oldPos.getChainID() != it->second.getChainID())
                        continue;
                if (oldPos.getChainPos() != it->second.getChainPos())
                        continue;

                it->second = newPos;
        }
}

void NodeChainContainer::printIndex() const
{
        for (auto it : index)
                cout << it.first <<  " " << it.second.getChainID() << " " << it.second.getChainPos() << endl;
}

size_t NodeChainContainer::updateNodeChains(const NodeChain& seedPattern,
                                            const NodeChain& oldPattern,
                                            const NodeChain& newPattern)
{
        // make sure the seed pattern is the initial part of the old pattern
        assert(seedPattern.size() <= oldPattern.size());
        for (size_t i = 0; i < seedPattern.size(); i++)
                assert(seedPattern[i] == oldPattern[i]);

        size_t numChanged = 0;

        vector<NodeChainPos> occ;
        findOcc(seedPattern, occ);

        // sort occurrences by chain identifier and position
        sort(occ.begin(), occ.end(), sortByPosition);

        const NodeChain newPatternRC = newPattern.getReverseComplement();

        // replace the patterns with the new patterns while updating the index
        for (size_t i = 0; i < occ.size(); i++) {
                size_t id = occ[i].getChainID();
                NodeChain& nc = (*this)[id];

                NodeChain copy;
                copy.reserve(nc.size());
                copy.setCount(nc.getCount());
                numChanged += nc.getCount();

                size_t srcPos = 0, dstPos = 0;

                do {
                        ssize_t pos = occ[i].getChainPos();
                        size_t patBegin = (occ[i].getReverse()) ? max(pos + 1 - (ssize_t)oldPattern.size(), (ssize_t)0) : pos;
                        size_t patEnd = (occ[i].getReverse()) ? pos + 1 : min(pos + oldPattern.size(), nc.size());

                        const NodeChain& newRep = occ[i].getReverse() ? newPatternRC : newPattern;

                        // copy the piece before the pattern occurrence
                        for ( ; srcPos < patBegin; srcPos++, dstPos++) {
                                copy.push_back(nc[srcPos]);
                                updateIndex(nc[srcPos], NodeChainPos(id, srcPos), NodeChainPos(id, dstPos));
                        }

                        // delete old pattern index
                        for (size_t i = 0; srcPos < patEnd; i++, srcPos++)
                                deleteIndex(nc[srcPos], NodeChainPos(id, srcPos));

                        // insert the new pattern (or its reverse complement)
                        for (size_t i = 0; i < newRep.size(); i++, dstPos++) {
                                copy.push_back(newRep[i]);
                                insertIndex(newRep[i], NodeChainPos(id, dstPos));
                        }

                        // check whether we're finished with this chain
                        if ( ((i+1) >= occ.size()) || (occ[i+1].getChainID() != id) )
                                break;

                        i++;

                } while (true);

                // copy the remainder past the last pattern occurrence
                for ( ; srcPos < nc.size(); srcPos++, dstPos++) {
                        copy.push_back(nc[srcPos]);
                        updateIndex(nc[srcPos], NodeChainPos(id, srcPos), NodeChainPos(id, dstPos));
                }

                nc = copy;
        }

        return numChanged;
}

size_t NodeChainContainer::processReduction(const NodeChain& reduction)
{
        size_t result = 0;

        // Assume a valid reduction
        assert(reduction.size() >= 3);

        // Assume the transformation A - B - ... - X - Y to A - Y
        NodeChain searchPattern = reduction;
        searchPattern.pop_back();

        NodeChain newPattern;
        newPattern.push_back(reduction.front());

        NodeChain seed(vector<NodeID>(searchPattern.begin(), searchPattern.begin() + 2));

        //cout << "Replacing: " << searchPattern << " by " << newPattern << endl;

        // Find all instances of "A - B - ..." and replace by "A"
        result += updateNodeChains(seed, searchPattern, newPattern);

        searchPattern = reduction.getReverseComplement();
        newPattern.clear();
        newPattern.push_back(searchPattern.front());
        newPattern.push_back(searchPattern.back());

        seed.clear();
        seed.push_back(searchPattern[0]);
        seed.push_back(searchPattern[1]);

        //cout << "Replacing: " << searchPattern << " by " << newPattern << endl;

        // Find all instances of "Y~ - X~ - ..." and replace by "Y~ - A~"
        result += updateNodeChains(seed, searchPattern, newPattern);

        return result;
}

size_t NodeChainContainer::processConcatentation(const NodeChain& concatenation)
{
        size_t result = 0;

        // Assume a valid concatenation
        if (concatenation.size() < 2)
                return result;

        NodeChain newPattern;
        newPattern.push_back(concatenation.front());

        for (size_t i = 0; i < concatenation.size(); i++) {
                NodeChain searchPattern(vector<NodeID>(concatenation.begin() + i, concatenation.end()));
                NodeChain seedPattern;
                seedPattern.push_back(searchPattern.front());
                result += updateNodeChains(seedPattern, searchPattern, newPattern);
        }

        return result;
}


void NodeChainContainer::removeSubChain(const NodeChainPos& pos)
{
        NodeChain& nc = at(pos.getChainID());

        for (size_t i = 0; i < nc.size(); i++)
                deleteIndex(nc[i], NodeChainPos(pos.getChainID(), i));
        nc.clear();
}

void NodeChainContainer::smoothPath(NodeID nodeID)
{
        // find all occurences of nodeID
        vector<NodeChainPos> occ;

        for (auto it = index.find(abs(nodeID)); it != index.end(); it++) {
                if (it->first != abs(nodeID))
                        break;

                const size_t& id = it->second.getChainID();
                const size_t& pos = it->second.getChainPos();
                bool reverse = at(id)[pos] == nodeID;

                occ.push_back(NodeChainPos(id, pos, reverse));
        }

        /*cout << "Before: " << endl;
        for (auto pos : occ) {
                for ( ; validNCP(pos); pos++)
                        cout << getNodeID(pos) << " ";
                cout << endl;
        }*/

        list<NodeChainPos> fwdIt(occ.begin(), occ.end());
        while(!fwdIt.empty()) {

                // check for consensus at the current position
                map<NodeID, size_t> consensus;
                for (auto it = fwdIt.begin(); it != fwdIt.end(); ++it)
                        consensus[getNodeID(*it)] += getCount(*it);

                // find the nodeID with the highest count
                NodeID consensusID = 0; size_t maxCount = 0;
                for (auto it : consensus)
                        if (it.second > maxCount) {
                                consensusID = it.first;
                                maxCount = it.second;
                        }

                // modify the nodeChains with alternate content
                for (auto it = fwdIt.begin(); it != fwdIt.end(); ) {
                        if (getNodeID(*it) != consensusID) {
                                removeSubChain(*it);
                                it = fwdIt.erase(it);
                        } else {
                                ++(*it);        // move to the next position

                                if (!validNCP(*it))
                                        it = fwdIt.erase(it);
                                else
                                        ++it;
                        }
                }
        }

        /*cout << "After: " << endl;
        for (auto pos : occ) {
                for ( ; validNCP(pos); pos++)
                        cout << getNodeID(pos) << " ";
                cout << endl;
        }*/
}

void NodeChainContainer::printPER() const
{
        cout << PER.size() << endl;
        for (auto it : PER) {
                int64_t idx1 = it.first.first, idx2 = it.first.second;

                if (idx1 == idx2)
                        continue;

                NodeChain nc1 = at(abs(idx1));
                if (idx1 < 0)
                        nc1 = nc1.getReverseComplement();

                NodeChain nc2 = at(abs(idx2));
                if (idx2 < 0)
                        nc2 = nc2.getReverseComplement();

                if (nc1 == nc2)
                        continue;
                if (nc1.empty() || nc2.empty())
                        continue;

                cout << nc1 << " --> " << nc2 << " (" << it.second << ")" << endl;
        }

        cout << " === " << endl;

        // create a node chain tree
        vector<NodeChain> nodeChainTree = getNodeChainTree(1279, 1323);

        for (auto it : nodeChainTree)
                cout << it << endl;
        cout << " === " << endl;
}
