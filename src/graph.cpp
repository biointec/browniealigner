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

#include "graph.h"
#include "graphaln.h"
#include "kmeroverlap.h"
#include "settings.h"
#include "nodeendstable.h"
#include "library.h"
#include <cmath>

#include <mutex>
#include <queue>
#include <iomanip>

using namespace std;

DSNode* SSNode::nodes = NULL;
const DBGraph* DBGraph::graph = NULL;

// ============================================================================
// GRAPH STATISTICS
// ============================================================================

std::ostream &operator<<(std::ostream &out, const GraphStats &stats)
{
        out << "Number of nodes: " << stats.numNodes << "\n";
        out << "Number of arcs: " << stats.numArcs << "\n";
        out << "N50: " << stats.N50 << "\n";
        out << "Total size (kmers): " << stats.totMargLength;

        return out;
}

bool sortNodeByLength(const NodeID& left, const NodeID& right)
{
        return DBGraph::graph->getDSNode(left).getMarginalLength() <
               DBGraph::graph->getDSNode(right).getMarginalLength();
}

// ============================================================================
// PARALLEL GRAPH CLASSES
// ============================================================================

void ParGraph::getNodeChunk(size_t& chunkOffset, size_t& chunkSize)
{
        // lock the mutex
        std::unique_lock<std::mutex> lock(mutex);

        // no more work: send termination signal (chunkSize == 0)
        if (currOffset > numNodes) {
                chunkOffset = chunkSize = 0;
                lock.unlock();
                return;
        }

        // assign a workload
        chunkOffset = currOffset;
        chunkSize = min(targetChunkSize, numNodes + 1 - currOffset);

        currOffset += chunkSize;

        double perc = 100.0 * (double)chunkOffset / (double)numNodes;
        cout << std::fixed << std::setprecision(1) << "\tProcessing graph (" << perc << "%)\r";
        cout.flush();

        lock.unlock();
}

// ============================================================================
// GRAPH CLASS
// ============================================================================

DBGraph::DBGraph(const Settings& settings) : settings(settings),
nodes(NULL), arcs(NULL), numNodes(0), numArcs(0), alignment(2, 1, -1, -3) {
        DBGraph::graph = this;
}

NodeID DBGraph::getFirstValidNode(NodeID seed)
{
        for (NodeID id = seed; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                #ifdef DEBUG
           /*     if (trueMult[id] == 0)
                        continue;*/
                #endif
                if (getSSNode(id).isValid())
                        return id;
        }
        return 0;
}

void DBGraph::sanityCheck()
{
        // sanity check
        for (NodeID i = -numNodes; i <= numNodes; i++) {
                if (i == 0)     // skip node 0, doesn't exist
                        continue;

                // shortcuts
                SSNode n = getSSNode(i);

                string sequence = n.getSequence();
                if (!n.isValid()) {
                        if ((n.getNumLeftArcs() != 0) || (n.getNumRightArcs() != 0))
                                cerr << "\t\tNode " << n.getNodeID()
                                << " is invalid but has arcs." << endl;
                }
                // check the continuity of the kmers
                Kmer firstKmer(sequence);
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        SSNode t = getSSNode(it->getNodeID());
                        if (!t.isValid())
                                cerr << "\t\tNode " << n.getNodeID()
                                << " has a left arc to an invalid node." << endl;

                        assert(t.getRightArc(i)->getCoverage() == it->getCoverage());
                        //  assert(it->getCoverage() > 0);
                        Kmer tKmer = t.getRightKmer();
                        tKmer.pushNucleotideRight(firstKmer.peekNucleotideRight());
                        if (firstKmer != tKmer)
                                cerr << "\tError in the continuity of the nodes" << endl;
                        assert(firstKmer == tKmer);
                }
                Kmer lastKmer(sequence, sequence.size()-Kmer::getK());
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        SSNode t = getSSNode(it->getNodeID());
                        if (!t.isValid())
                                cerr << "\t\tNode " << n.getNodeID()
                                << " has a right arc to an invalid node." << endl;

                        if (t.getLeftArc(i)->getCoverage() != it->getCoverage())
                                cerr << "\t\tNode " << n.getNodeID()
                                << " has a right arc with cov " << (int)it->getCoverage()
                                << " to node " << t.getNodeID() << " whereas that node"
                                << " has a left arc with cov " << (int)t.getLeftArc(i)->getCoverage() << endl;

                        assert(t.getLeftArc(i)->getCoverage() == it->getCoverage());

                        Kmer tKmer = t.getLeftKmer();
                        tKmer.pushNucleotideLeft(lastKmer.peekNucleotideLeft());

                        if (lastKmer != tKmer)
                                cerr << "\tError in the continuity of the nodes" << endl;
                        assert (lastKmer == tKmer);
                }
        }
}

DBGraph::~DBGraph()
{
        delete [] arcs;
        delete [] nodes;
}

bool DBGraph::getLeftUniqueSSNode(const SSNode &node, SSNode &leftNode) const
{
        if (node.getNumLeftArcs() != 1)
                return false;

        leftNode = getSSNode(node.leftBegin()->getNodeID());

        if (leftNode.getNumRightArcs() > 1)
                return false;

        return true;
}

bool DBGraph::getRightUniqueSSNode(const SSNode &node, SSNode &rightNode) const
{
        if (node.getNumRightArcs() != 1)
                return false;

        rightNode = getSSNode(node.rightBegin()->getNodeID());

        if (rightNode.getNumLeftArcs() > 1)
                return false;

        return true;
}

void DBGraph::increaseCoverage(const NodeEndRef &left, const NodeEndRef &right)
{
        SSNode lNode = getSSNode(left.getNodeID());
        SSNode rNode = getSSNode(right.getNodeID());

        // Avoid increasing the coverage of a non-existing intra-node arc.
        // This might arise when the node consists of two overlapping kmers
        // or when (the first part of) the node is a palindrome
        if (abs(left.getNodeID()) == abs(right.getNodeID()))
                if (lNode.getLength() > Kmer::getK())
                        if (lNode.getLeftKmer() == left.getKmer())
                                return;

        lNode.getRightArc(right.getNodeID())->incReadCov();
        rNode.getLeftArc(left.getNodeID())->incReadCov();
}

void DBGraph::convertNodesToString(const vector<NodeID> &nodeSeq,
                                   string &output)
{
        output.clear();

        if (nodeSeq.empty())
                return;

        size_t size = 0;
        for (size_t i = 0; i < nodeSeq.size(); i++)
                size += getSSNode(nodeSeq[i]).getLength();
        size -= (nodeSeq.size() - 1) * (Kmer::getK() - 1);

        output = getSSNode(nodeSeq[0]).getSequence();
        for (size_t i = 1; i < nodeSeq.size(); i++)
                output.append(getSSNode(nodeSeq[i]).getSequence().substr(Kmer::getK() - 1));

        assert(size == output.size());
}

GraphStats DBGraph::getGraphStats()
{
        vector<size_t> nodeLength;

        size_t numNodes = 0;            // number of (valid) nodes in the graph
        size_t numArcs = 0;             // number of (valid) arcs in the graph
        size_t N50 = 0;                 // N50 of the nodes
        size_t totMargLength = 0;       // total marginal length of all nodes
        size_t totLength = 0;           // total length of all nodes

        for (NodeID id = 1; id <= this->numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                numNodes++;

                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++)
                        numArcs++;

                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++)
                        numArcs++;

                totMargLength += node.getMarginalLength();
                totLength += node.getLength();
                nodeLength.push_back(node.getLength());
        }

        sort(nodeLength.begin(), nodeLength.end());

        size_t currLength = 0;
        for (size_t i = 0; i < nodeLength.size(); i++) {
                currLength += nodeLength[i];
                if (currLength >= 0.5*totLength) {
                        N50 = nodeLength[i];
                        break;
                }
        }

        GraphStats stats;
        stats.setMetrics(numNodes, numArcs, N50, totMargLength);
        return stats;
}

void DBGraph::writeCytoscapeGraph(const std::string& filename,
                                  NodeID seedNodeID, size_t maxDepth) const
{
        cout << "Writing cytoscape graph of all the nodes..." <<endl;
        // a map containing nodeIDs to handle + their depth
        priority_queue<PathDFS, vector<PathDFS>, PathDFSComp> nodeDepth;
        // set of nodes that were handled
        set<NodeID> nodesHandled;

        // if default input values are provided, write the entire graph
        if (seedNodeID == 0) {
                for (NodeID id = -numNodes; id <= numNodes; id++) {
                        if (id == 0)
                                continue;
                        if (!getSSNode(id).isValid())
                                continue;
                        nodeDepth.push(PathDFS(id, 0, 0));
                }
        } else {        // else check if seedNode is valid
                if ((abs(seedNodeID) > numNodes) || !getSSNode(seedNodeID).isValid()) {
                        cerr << "WARNING: trying to use an invalid node as a "
                                "seed in writeCytoscapeGraph!" << endl;
                        return;
                }
                nodeDepth.push(PathDFS(seedNodeID, 0, 0));
        }

        // map of all nodes in the local graph and its depth
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs" << " for writing" << endl;
        ofs << "Source node\tTarget node\tArc coverage" << endl;

        // A) write all arcs
        while (!nodeDepth.empty()) {
                // get and erase the current node
                PathDFS currTop = nodeDepth.top();
                nodeDepth.pop();
                NodeID thisID = currTop.nodeID;
                size_t thisDepth = currTop.length;

                // if the node was already handled, skip
                if (nodesHandled.find(thisID) != nodesHandled.end())
                        continue;

                // if we're too far in the graph, stop
                if (thisDepth > maxDepth) {
                        nodesHandled.insert(thisID);
                        continue;
                }

                // write the right arcs
                SSNode thisNode = getSSNode(thisID);
                for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
                        SSNode rNode = getSSNode(it->getNodeID());
                        if (!rNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        ofs << thisID << "\t" << it->getNodeID() << "\t" << it->getCoverage() << "\n";
                        PathDFS nextTop(it->getNodeID(), 0, thisDepth + 1);
                        nodeDepth.push(nextTop);
                }

                // write the left arcs
                for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
                        SSNode lNode = getSSNode(it->getNodeID());
                        if (!lNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        if (it->getNodeID() != thisID)  // avoid writing this arc twice
                                ofs << it->getNodeID() << "\t" << thisID << "\t" << it->getCoverage() << "\n";
                        PathDFS nextTop(it->getNodeID(), 0, thisDepth + 1);
                        nodeDepth.push(nextTop);
                }
                if (thisNode.getNumLeftArcs() ==0 && thisNode.getNumRightArcs()==0){
                     ofs << thisID << "\t" << ""<< "\t" << "" << "\n";

                }

                // mark this node as handled
                nodesHandled.insert(thisID);

        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes" << " for writing" << endl;
        ofs << "Node ID\tMarginal length\tNum left arcs\tNum right arcs\tTrue multiplicity\tEstimated multiplicity"
               "\tKmer coverage\tRead start coverage\tSequence" << "\n";
        for (set<NodeID>::iterator it = nodesHandled.begin(); it != nodesHandled.end(); it++) {
                SSNode node = getSSNode(*it);

#ifdef DEBUG
                int thisTrueMult = (trueMult.empty()) ? 0 : trueMult[abs(*it)];
#else
                int thisTrueMult = 0;
#endif


                ofs << *it << "\t" << node.getMarginalLength() << "\t"
                    << (int)node.getNumLeftArcs() << "\t" << (int)node.getNumRightArcs() << "\t"
                    << thisTrueMult << "\t" << "0" << "\t"
                    << node.getAvgKmerCov()
                    << "\t" << double(node.getReadStartCov()/node.getMarginalLength())
                    << "\t" << node.getSequence() << "\n";
        }
        ofs.close();
}
void DBGraph::writeCytoscapeGraph(const std::string& filename,
                                  std::set<int> &nodeSet, size_t maxDepth) const
{
        // a map containing nodeIDs to handle + their depth
        priority_queue<PathDFS, vector<PathDFS>, PathDFSComp> nodeDepth;
        // set of nodes that were handled
        set<NodeID> nodesHandled;

        // if default input values are provided, write the entire graph
        set <int> ::iterator it;
        for (it =nodeSet.begin();it!=nodeSet.end();it++)
        {        // else check if seedNode is valid
                NodeID seedNodeID = *it;
                if ((abs(seedNodeID) > numNodes) || !getSSNode(seedNodeID).isValid()) {
                        cerr << "WARNING: trying to use an invalid node as a "
                                "seed in writeCytoscapeGraph!" << endl;
                        return;
                }
                nodeDepth.push(PathDFS(seedNodeID, 0, 0));
        }

        // map of all nodes in the local graph and its depth
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs" << " for writing" << endl;
        ofs << "Source node\tTarget node\tArc coverage" << endl;

        // A) write all arcs
        while (!nodeDepth.empty()) {
                // get and erase the current node
                PathDFS currTop = nodeDepth.top();
                nodeDepth.pop();
                NodeID thisID = currTop.nodeID;
                size_t thisDepth = currTop.length;

                // if the node was already handled, skip
                if (nodesHandled.find(thisID) != nodesHandled.end())
                        continue;

                // if we're too far in the graph, stop
                if (thisDepth > maxDepth) {
                        nodesHandled.insert(thisID);
                        continue;
                }

                // write the right arcs
                SSNode thisNode = getSSNode(thisID);
                for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
                        SSNode rNode = getSSNode(it->getNodeID());
                        if (!rNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        ofs << thisID << "\t" << it->getNodeID() << "\t" << it->getCoverage() << "\n";
                        PathDFS nextTop(it->getNodeID(), 0, thisDepth + 1);
                        nodeDepth.push(nextTop);
                }

                // write the left arcs
                for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
                        SSNode lNode = getSSNode(it->getNodeID());
                        if (!lNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        if (it->getNodeID() != thisID)  // avoid writing this arc twice
                                ofs << it->getNodeID() << "\t" << thisID << "\t" << it->getCoverage() << "\n";
                        PathDFS nextTop(it->getNodeID(), 0, thisDepth + 1);
                        nodeDepth.push(nextTop);
                }

                // mark this node as handled
                nodesHandled.insert(thisID);

        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes" << " for writing" << endl;
        ofs << "Node ID\tMarginal length\tNum left arcs\tNum right arcs\tTrue multiplicity\tEstimated multiplicity"
               "\tKmer coverage\tRead start coverage\tSequence\tinPath" << "\n";
        for (set<NodeID>::iterator it = nodesHandled.begin(); it != nodesHandled.end(); it++) {
                SSNode node = getSSNode(*it);

#ifdef DEBUG
                int thisTrueMult = (trueMult.empty()) ? 0 : trueMult[abs(*it)];
#else
                int thisTrueMult = 0;
#endif
                const bool is_in = nodeSet.find(*it) != nodeSet.end();
                ofs << *it << "\t" << node.getMarginalLength() << "\t"
                    << (int)node.getNumLeftArcs() << "\t" << (int)node.getNumRightArcs() << "\t"
                    << thisTrueMult << "\t" << "0" << "\t"
                    << node.getAvgKmerCov()
                    << "\t" << double(node.getReadStartCov()/node.getMarginalLength())
                    << "\t" << node.getSequence()
                    <<"\t"<<is_in <<"\n";
        }
        ofs.close();
}

void DBGraph::loadGraph(const string& nodeFilename,
                        const string& arcFilename,
                        const string& metaDataFilename)
{
        // auxiliary variables
        string dS, descriptor;
        int dI, length;
        Coverage kmerCov, readStartCov;

        // read the metadata
        ifstream metaDataFile(metaDataFilename.c_str());
        if (!metaDataFile)
                throw ios_base::failure("Can't open " + metaDataFilename);
        metaDataFile >> numNodes >> numArcs;
        numValidNodes = numNodes;
        numValidArcs = numArcs;
        metaDataFile.close();

        NodeEndTable table(settings.isDoubleStranded(), 2*numNodes);

        // A) create the nodes
        ifstream nodeFile(nodeFilename.c_str());
        if (!nodeFile)
                throw ios_base::failure("Can't open " + nodeFilename);

        nodes = new DSNode[numNodes+1];
        SSNode::setNodePointer(nodes);
        for (NodeID id = 1; id <= numNodes; id++) {
                // read the node info
                nodeFile >> dS >> dI >> length >> kmerCov >> readStartCov >> descriptor;

                DSNode& node = getDSNode(id);
                node.setSequence(descriptor);

                node.setKmerCov(kmerCov);
                node.setReadStartCov(readStartCov);

                Kmer firstKmer = node.getLeftKmer();
                Kmer finalKmer = node.getRightKmer();

                if (node.getMarginalLength() > 1) {
                        if (!table.insert(firstKmer, id))
                                cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
                        if (!table.insert(finalKmer, id))
                                cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
                } else if (!table.insert(firstKmer, id))
                        cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
        }

        nodeFile.close();

        // B) create the arcs
        ifstream arcFile(arcFilename.c_str());
        if (!arcFile)
                throw ios_base::failure("Can't open " + arcFilename);

        // +2 because index 0 isn't used, final index denotes 'end'.
        arcs = new Arc[numArcs+2];
        DSNode::setArcsPointer(arcs);

        ArcID arcOffset = 1;
        for (NodeID i = 1; i <= numNodes; i++) {
                int dI, bfLeft, bfRight;

                // read the arc info
                arcFile >> dI >> bfLeft >> bfRight;
                KmerOverlap overlap((bfLeft << 4) + bfRight);

                int numLeftArcs = overlap.getNumLeftOverlap();
                int numRightArcs = overlap.getNumRightOverlap();

                DSNode& node = getDSNode(i);

                node.setNumLeftArcs(numLeftArcs);
                node.setNumRightArcs(numRightArcs);

                node.setFirstLeftArcID(arcOffset);
                node.setFirstRightArcID(arcOffset+numLeftArcs);

                // connect the left arcs alphabetically (ACGT)
                Kmer firstKmer = node.getLeftKmer();
                for (NucleotideID j = 0; j < 4; j++) {
                        char n = Nucleotide::nucleotideToChar(j);
                        if (!overlap.hasLeftOverlap(n))
                                continue;

                        Kmer kmer = firstKmer;
                        kmer.pushNucleotideLeft(n);

                        NodeEndRef ref = table.find(kmer);
                        if (ref.first == table.end())
                                throw ios_base::failure("Mismatch between nodes"
                                "and arc file.");
                        int arcCov;
                        arcFile >> arcCov;
                        arcs[arcOffset].setCoverage(arcCov);
                        arcs[arcOffset++].setNodeID(ref.getNodeID());
                }

                // connect the right arcs alphabetically (ACGT)
                Kmer finalKmer = node.getRightKmer();
                for (NucleotideID j = 0; j < 4; j++) {
                        char n = Nucleotide::nucleotideToChar(j);
                        if (!overlap.hasRightOverlap(n))
                                continue;

                        Kmer kmer = finalKmer;
                        kmer.pushNucleotideRight(n);

                        NodeEndRef ref = table.find(kmer);
                        if (ref.first == table.end())
                                throw ios_base::failure("Mismatch between nodes"
                                "and arc file.");
                        int arcCov;
                        arcFile >> arcCov;
                        arcs[arcOffset].setCoverage(arcCov);
                        arcs[arcOffset++].setNodeID(ref.getNodeID());
                }
        }

        arcFile.close();

        if(arcOffset != numArcs+1)
                throw ios_base::failure("Mismatch between nodes and arc file.");
}

void DBGraph::writeGraph(const std::string& nodeFilename,
                         const std::string& arcFilename,
                         const std::string& metaDataFilename) const
{
        ofstream nodeFile(nodeFilename.c_str());
        ofstream arcFile(arcFilename.c_str());

        size_t numExtractedNodes = 0, numExtractedArcs = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                numExtractedNodes++;

                nodeFile << ">NODE" << "\t" << id << "\t"
                << node.getLength() << "\t" << node.getKmerCov()
                << "\t" << node.getReadStartCov() << "\n"
                << node.getSequence() << "\n";

                KmerOverlap ol;
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                        char c = getSSNode(it->getNodeID()).getRightKmer().peekNucleotideLeft();
                        ol.markLeftOverlap(c);
                }

                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        char c = getSSNode(it->getNodeID()).getLeftKmer().peekNucleotideRight();
                        ol.markRightOverlap(c);
                }

                arcFile << numExtractedNodes << "\t" << (int)ol.getLeftOverlap()
                        << "\t" << (int)ol.getRightOverlap();

                // write the left arcs alphabetically (ACGT)
                for (NucleotideID j = 0; j < 4; j++) {
                        char n = Nucleotide::nucleotideToChar(j);
                        if (!ol.hasLeftOverlap(n))
                                continue;

                        NodeID leftID = node.getLeftArcNodeID(n);
                        arcFile << "\t" << (int)node.getLeftArc(leftID)->getCoverage();
                }

                // write the right arcs alphabetically (ACGT)
                for (NucleotideID j = 0; j < 4; j++) {
                        char n = Nucleotide::nucleotideToChar(j);
                        if (!ol.hasRightOverlap(n))
                                continue;

                        NodeID rightID = node.getRightArcNodeID(n);
                        arcFile << "\t" << (int)node.getRightArc(rightID)->getCoverage();
                }

                arcFile << endl;

                numExtractedArcs += ol.getNumLeftOverlap() + ol.getNumRightOverlap();
        }

        nodeFile.close();
        arcFile.close();

        ofstream metadataFile(metaDataFilename.c_str());
        metadataFile << numExtractedNodes << "\t" << numExtractedArcs << endl;
        metadataFile.close();

        cout << "Extracted " << numExtractedNodes << " nodes and "
        << numExtractedArcs << " arcs." << endl;
}

void DBGraph::loadGraphBin(const std::string& nodeFilename,
                           const std::string& arcFilename,
                           const std::string& metaDataFilename)
{
        // auxiliary variables
        string dS, descriptor;

        // read the metadata
        ifstream metaDataFile(metaDataFilename.c_str());
        if (!metaDataFile)
                throw ios_base::failure("Can't open " + metaDataFilename);
        metaDataFile >> numNodes >> numArcs;
        numValidNodes = numNodes;
        numValidArcs = numArcs;
        metaDataFile.close();

        // A) create the nodes
        ifstream nodeFile(nodeFilename.c_str(), ios::binary);
        if (!nodeFile)
                throw ios_base::failure("Can't open " + nodeFilename);

        nodes = new DSNode[numNodes+1];
        SSNode::setNodePointer(nodes);
        for (NodeID id = 1; id <= numNodes; id++) {
                // read the node info

                DSNode& node = getDSNode(id);
                node.read(nodeFile);
        }
        nodeFile.close();

        // B) create the arcs
        ifstream arcFile(arcFilename.c_str());
        if (!arcFile)
                throw ios_base::failure("Can't open " + arcFilename);

        // +2 because index 0 isn't used, final index denotes 'end'.
        arcs = new Arc[numArcs+2];
        DSNode::setArcsPointer(arcs);
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].read(arcFile);
        arcFile.close();
}

void DBGraph::writeGraphBin(const std::string& nodeFilename,
                            const std::string& arcFilename,
                            const std::string& metaDataFilename) const
{
        // A) Write node file
        ofstream nodeFile(nodeFilename.c_str(), ios::binary);
        for (NodeID id = 1; id <= numNodes; id++) {
                DSNode& node = getDSNode(id);
                // write the node contents
                node.write(nodeFile);
        }
        nodeFile.close();

        // B) Write arc file
        ofstream arcFile(arcFilename.c_str(), ios::binary);
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].write(arcFile);
        arcFile.close();

        // C) Write metadata file
        ofstream metadataFile(metaDataFilename.c_str());
        metadataFile << numNodes << "\t" << numArcs << endl;
        metadataFile.close();

        cout << "Wrote " << numNodes << " nodes and "
             << numArcs << " arcs" << endl;
}

void DBGraph::writeGraphFasta() const
{
        vector<size_t> nodeLengths;
        string nodeFileName = settings.getTempDirectory() + "/DBGraph.fasta";
        ofstream nodeFile(nodeFileName);

        size_t numExtractedNodes = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                numExtractedNodes++;

                nodeFile << ">NODE" << "\t" << numExtractedNodes << "\t"
                << node.getLength();


                KmerOverlap ol;
                nodeFile << "\t" << (int)node.getNumLeftArcs();
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                        nodeFile << "\t" << it->getNodeID();
                }

                nodeFile << "\t" << (int)node.getNumRightArcs();
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        nodeFile << "\t" << it->getNodeID();
                }
                nodeFile << "\n" << node.getSequence() << "\n";;
        }

        nodeFile.close();
}

bool DBGraph::dijkstra(NodeID srcID, NodeID dstID, size_t maxLen,
                       vector<size_t>& dist_v, vector<bool>& visited_v) const
{
        vector<NodeID> modifiedNodes_v;        // which entries are modified?
        bool returnValue = false;

        // seed priority queue with source node
        priority_queue<PathDFS, vector<PathDFS>, PathDFSComp> heap;
        heap.push( PathDFS(srcID, 0, 0) );
        modifiedNodes_v.push_back(srcID);

        while(!heap.empty()) {
                // get the top node
                PathDFS currTop = heap.top();
                heap.pop();
                NodeID currID = currTop.nodeID;
                size_t currLength = currTop.length;
                SSNode curr = getSSNode(currID);

                // if node was visited before, get out: this is possible because
                // we don't update distances in the PQ but rather insert doubles
                if (visited_v[currID + numNodes])
                        continue;
                visited_v[currID + numNodes] = true;

                for (ArcIt it = curr.rightBegin(); it != curr.rightEnd(); it++) {
                        NodeID nextID = it->getNodeID();
                        SSNode next = getSSNode(nextID);
                        size_t nextLength = currLength + next.getMarginalLength();

                        // we found a path to the destination
                        if (nextID == dstID) {
                                returnValue = true;
                                goto exitRoutine;
                        }

                        // if the node was visited before OR path is too long
                        if (visited_v[nextID + numNodes] || nextLength > maxLen)
                                continue;

                        // if we can reach the node in a shorter distance
                        if (nextLength < dist_v[nextID + numNodes]) {
                                dist_v[nextID + numNodes] = nextLength;
                                heap.push( PathDFS(nextID, 0, nextLength) );
                                modifiedNodes_v.push_back(nextID);
                        }
                }
        }

        // label definition to break out of nested loops
        exitRoutine:

        // undo all modifications to dist_v and visited_v
        for (auto it : modifiedNodes_v) {
                dist_v[it + numNodes] = numeric_limits<size_t>::max();
                visited_v[it + numNodes] = false;
        }

        return returnValue;
}

bool DBGraph::findPath(const NodePosPair& srcNpp, const NodePosPair& dstNpp,
                       size_t maxLen, vector<size_t>& dist_v,
                       vector<bool>& visited_v) const
{
        // shortcut notation
        NodeID srcID = srcNpp.getNodeID();
        NodeID dstID = dstNpp.getNodeID();

        // assert valid positions in the graph
        assert(getSSNode(srcID).isValid());
        assert(getSSNode(dstID).isValid());
        assert(srcNpp.getPosition() < getSSNode(srcID).getMarginalLength());
        assert(dstNpp.getPosition() < getSSNode(dstID).getMarginalLength());

        // handle case in which path exists within a node
        if ((srcID == dstID) && (srcNpp.getPosition() <= dstNpp.getPosition()))
                return (srcNpp.getPosition() + maxLen >= dstNpp.getPosition());

        // correct maxLen for source node
        size_t correction = getSSNode(srcID).getMarginalLength() -
                srcNpp.getPosition() + dstNpp.getPosition();
        if (correction > maxLen)
                return false;
        maxLen -= correction;

        // use the dijkstra algorithm to figure out if there exists a path
        return dijkstra(srcID, dstID, maxLen, dist_v, visited_v);
}

void DBGraph::performReduction(const NodeChain& reduction)
{
        NodeID firstID = reduction.front();
        NodeID nextID = reduction[1];
        NodeID lastID = reduction.back();
        NodeID prevID = reduction[reduction.size() - 2];

        SSNode first = getSSNode(firstID);
        SSNode next = getSSNode(nextID);
        SSNode last = getSSNode(lastID);
        SSNode prev = getSSNode(prevID);

        // rewire the arcs
        first.replaceRightArc(nextID, lastID);
        next.deleteLeftArc(firstID);
        last.replaceLeftArc(prevID, firstID);
        prev.deleteRightArc(lastID);

        first.getRightArc(lastID)->setCoverage(last.getLeftArc(firstID)->getCoverage());

        // FIXME : update coverage !!

        /*size_t newKmerCov = first.getKmerCov();
        size_t newReadStartCov = first.getReadStartCov();
        for (size_t i = 1; i < reduction.size() - 1; i++) {
                newKmerCov += getSSNode(reduction[i]).getKmerCov();
                newReadStartCov += getSSNode(reduction[i]).getReadStartCov();
                getSSNode(reduction[i]).deleteAllLeftArcs();
                getSSNode(reduction[i]).deleteAllRightArcs();
                getSSNode(reduction[i]).invalidate();
        }

        first.setKmerCov(newKmerCov);
        first.setReadStartCov(newReadStartCov);*/

        string str;
        convertNodesToString(vector<NodeID>(reduction.begin(), reduction.begin() + reduction.size() - 1), str);

#ifdef DEBUG
        if (!trueMult.empty()) {
                for (size_t i = 1; i < reduction.size() - 1; i++)
                        trueMult[abs(reduction[i])] -= trueMult[abs(reduction[0])];
        }
#endif

        first.setSequence(str);
}

bool DBGraph::validateChain(const NodeChain& nc) const
{
        for (size_t i = 0; i < nc.size(); i++) {
                if (nc[i] == 0)
                        continue;

                SSNode node = getSSNode(nc[i]);

                if (!node.isValid()) {
                        cout << "Not valid" << endl;
                        return false;
                }
                if (i+1 == nc.size()) {
                        break;
                }
                if (node.getRightArc(nc[i+1]) == NULL) {
                        cout << "Not connected" << endl;
                        return false;
                }
        }

        return true;
}

void DBGraph::validateChainContainer(const NodeChainContainer& ncc)
{
        for (const NodeChain& nc : ncc)
                if (!validateChain(nc))
                        cout << nc << " is no longer valid in the graph!" << endl;
}

void DBGraph::findReductions(vector<NodeChain>& reductionv)
{
        reductionv.clear();

        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;

                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                if (node.getNumRightArcs() != 1)
                        continue;
                NodeID rightID = node.rightBegin()->getNodeID();

                vector<NodeChain> nbPaths = ncc.getNonBranchingPath(id, rightID);

                for (const auto& nbPath : nbPaths) {
                        // check whether the final node in the nbPath is part of a loop structure
                        bool isLoop = false;
                        for (size_t i = 0; i < nbPath.size() - 1; i++)
                                if (nbPath[i] == nbPath.back())
                                        isLoop = true;
                        if (isLoop)
                                continue;

                        // TODO : what about hairpins: A - B - A~ ??

                        // check whether the reverse-complement is OK
                        NodeChain nbPathRC = nbPath.getReverseComplement();

                        if (ncc.isNonBranchingPath(nbPathRC)) {
                                bool typeA = getSSNode(nbPath.back()).getNumLeftArcs() >= 2;
                                bool repr = nbPath < nbPathRC;

                                // save the reduction when it is
                                if (repr || typeA) {
                                        reductionv.push_back(nbPath);
                                }
                                break;
                        }
                }
        }
}

