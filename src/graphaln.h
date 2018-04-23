/***************************************************************************
 *   Copyright (C) 2017 Jan Fostier (jan.fostier@ugent.be)                 *
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

#ifndef GRAPHALN_H
#define GRAPHALN_H

#include <vector>

#include "alignment.h"
#include "markovChain.h"
// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class DBGraph;
class NodeAlignment;

typedef std::pair<NodeAlignment, size_t> QueueEl;   // <NodeAlignment, Depth>

enum SearchRes { OPTIMAL, EXHAUSTED, MULTIPLE };

// ============================================================================
// DIJKSTRA AUXILIARY CLASSES
// ============================================================================

class PathDFS {

public:
        NodeID nodeID;          // current node identifier
        NodeID prevID;          // previous node identifier
        size_t length;          // length to current node

        /**
         * Default constructor
         * @param nodeID node identifier
         * @param prevID previous node identifier
         * @param length length to the current node
         */
        PathDFS(NodeID nodeID, NodeID prevID, size_t length) :
                nodeID(nodeID), prevID(prevID), length(length) {};
};

struct PathDFSComp {

        /**
         * Compare two paths (by length) (strict weak ordering (irreflexive))
         */
        bool operator()(const PathDFS& f, const PathDFS& s) {
                return f.length > s.length;
        }
};

// ============================================================================
// NODE PROPERTY CLASS
// ============================================================================

// This class behaves as a vector where you can use negative indices. This is
// useful when using node identifiers as key as they can become negative.

template<class T>
class NodeProp : private std::vector<T>
{
private:
        size_t numNodes;        // number of nodes

public:
        /**
         * Default constructor
         * @param numNodes Number of nodes
         */
        NodeProp(size_t numNodes, const T& val) : numNodes(numNodes) {
                vector<T>::resize(2*numNodes + 1, val);
        }

        T& operator[](ssize_t n) {
                return vector<T>::operator[](n + numNodes);
        }

        const T& operator[](ssize_t n) const {
                return vector<T>::operator[](n + numNodes);
        }
};

// ============================================================================
// NODE ALIGNMENT CLASS
// ============================================================================

class NodeAlignment
{
public:
        NodeID nodeID;          // node identifier to which is aligned
        size_t nodeBegin;       // node begin position to which is aligned
        size_t nodeEnd;         // node end position to which is aligned

        size_t pattBegin;       // pattern begin position to which is aligned
        size_t pattEnd;         // pattern end position to which is aligned

        int score;              // alignment score
        float relScore;         // alignment score / maximum alignment score

        /**
         * Default constructor
         */
        NodeAlignment() : nodeID(0), nodeBegin(0), nodeEnd(0), pattBegin(0),
                pattEnd(0), score(0), relScore(0.0f) {}

        /**
         * Default constructor with initialization
         */
        NodeAlignment(NodeID nodeID, size_t nodeBegin, size_t nodeEnd,
                      size_t pattBegin, int pattEnd, int score, float relScore) :
                nodeID(nodeID), nodeBegin(nodeBegin), nodeEnd(nodeEnd),
                pattBegin(pattBegin), pattEnd(pattEnd), score(score),
                relScore(relScore) {}

        /**
         * Get the node position pair of the end of the alignment
         * @return The node position pair of the end of the alignment
         */
        NodePosPair getNodeEnd() {
                return NodePosPair(nodeID, nodeEnd);
        }

        /**
         * Operator < overloading
         * @param rhs Right hand size operand
         * @return True if this object has a lower relative score than rhs
         */
        bool operator< (const NodeAlignment& rhs) const {
                if (relScore != rhs.relScore)
                        return relScore < rhs.relScore;
                return score < rhs.score;
        }
};

// ============================================================================
// GRAPH ALIGNMENT CLASS
// ============================================================================

class GraphAlignment
{
private:
        std::vector<NodeAlignment> path;        // sequence of node alignments
        std::vector<int> pathScore;             // path score

public:
        /**
         * Default constructor
         */
        GraphAlignment() {}

        /**
         * Clear the graph alignment
         */
        void clear() {
                path.clear();
        }

        /**
         * Get the graph alignment score
         * @return The graph alignment score
         */
        int getScore() const {
                return pathScore.empty() ? 0 : pathScore.back();
        }

        /**
         * Get the depth of the path
         * @return The depth of the path (= number of nodes in the alignment)
         */
        size_t getDepth() const {
                return path.size();
        }

        /**
         * Get the length of the sequence that is actually aligned
         * @return Length of the sequence that is aligned
         */
        size_t getAlignedLength() const {
                return path.empty() ? 0 : path.back().pattEnd;
        }

        /**
         * Update graph alignment with new node alignment
         * @param nodeAln Node alignment
         */
        void updatePath(NodeAlignment nodeAln, size_t depth) {

                // update the path vector
                path.resize(depth);
                path[depth-1] = nodeAln;

                // update the path score vector
                pathScore.resize(depth);
                int prevPathScore = (depth > 1) ?
                        pathScore[depth-2] : 0;
                pathScore[depth-1] = prevPathScore + nodeAln.score;
        }

        /**
         * Get the sequence implied by the alignment
         * @param dbg Const-ref to the matching de Bruijn graph
         * @return Sequence string
         */
        std::string getSequence(const DBGraph& dbg) const;

        /**
         * Get the average kmer coverage for a given graph alignment
         * @param dbg Const-ref to the matching de Bruijn graph
         * @return The average kmer coverage
         */
        double getAvgKmerCov(const DBGraph& dbg) const;

        /**
         * Get a const-reference to the node alignments
         * @return A const-reference to the node alignments
         */
        const std::vector<NodeAlignment>& getPath() const {
                return path;
        }
        /**
        * returns the current visited nodes in a chain
        * @return A vector contains all the nodes in the current alignment path
        */
        void getNodeChain(vector<NodeID> &visitedNodes) ;


        void getKmerNPP(vector<NodePosPair>& kmerNPP) const {
                kmerNPP.clear();

                if (path.empty())
                        return;

                kmerNPP.reserve(path.back().pattEnd);

                // extract the k-mers from the first node
                const NodeAlignment& na = path[0];
                for (size_t i = na.nodeBegin; i < na.nodeEnd + 1 - Kmer::getK(); i++)
                        kmerNPP.push_back(NodePosPair(na.nodeID, i));

                // extract the k-mers from the other nodes
                for (size_t j = 1; j < path.size(); j++) {
                        const NodeAlignment& na = path[j];
                        for (size_t i = na.nodeBegin; i < na.nodeEnd; i++)
                                kmerNPP.push_back(NodePosPair(na.nodeID, i + 1 - Kmer::getK()));
                }
        }

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ga GraphAlignment object to display
         * @return Output file stream
         */
        friend std::ostream &operator<<(std::ostream &out,
                                        const GraphAlignment& ga);
};

// ============================================================================
// GRAPH ALIGNER CLASS
// ============================================================================

class GraphAligner {

private:
        const DBGraph& dbg;     // const-ref to the underlying dbg
        size_t maxVisits;       // maximum number of nodes visited
        NWAligner aligner;      // alignment class
        MarkovChainHandler &mch; // markov chain transition matrix
         /**
         * Visualize a node alignment
         */
        void visualizeNodeAln(const NodeAlignment& nodeAln,
                              const GraphAlignment& ga, const string& P) const;

        /**
         * Given a pattern and a partial graph alignment, return the maximum
         * score that could possible be obtained by completing the alignment
         * @param P pattern
         * @param ga Partial graph alignment
         * @param shortestPath Shortest path to destination (if applicable)
         */
        int getMaxAttainableScore(const std::string& P,
                                  const GraphAlignment& ga,
                                  size_t shortestPath = 0) const;

        /**
         * Get the length of a node
         * @param nodeID Node identifier
         * @return The length of a node
         */
        size_t nodeLength(NodeID nodeID) {
                return dbg.getSSNode(nodeID).getLength();
        }

        /**
         * Get the node position pair pointing to the begin of a node
         * @param nodeID Node identifier
         * @return Node position pair pointing to the begin of a node
         */
        NodePosPair getNodeBegin(NodeID nodeID) {
                return NodePosPair(nodeID, Kmer::getK() - 1);
        }

        /**
         * Get the node position pair pointing to the end of a node
         * @param nodeID Node identifier
         * @return Node position pair pointing to the end of a node
         */
        NodePosPair getNodeEnd(NodeID nodeID) {
                return NodePosPair(nodeID, dbg.getSSNode(nodeID).getLength());
        }

        /**
         * Best alignment of P to N allowing unaligned suffix in either P or N
         * "Unaligned Trailing Pattern - Unaligned Trailing Node"
         * @param patt Pattern P to align
         * @param pattBegin Begin position of the pattern P
         * @param nodeBegin Begin position of the node N
         * @return Node alignment object
         */
        NodeAlignment alnToNodeOpenEnd(const string& patt, size_t pattBegin,
                                      const NodePosPair& nodeBegin);

        /**
         * Best alignment of P to N allowing unaligned suffix in P
         * "Unaligned Trailing Pattern"
         * @param patt Pattern P to align
         * @param pattBegin Begin position of the pattern P
         * @param nodeBegin Begin position of the node N
         * @param nodeEnd End position of the node N
         * @return Node alignment object
         */
        NodeAlignment nodeAlnTo(const string& patt, size_t pattBegin,
                                   const NodePosPair& nodeBegin,
                                   const NodePosPair& nodeEnd);

        /**
         * Best alignment of P to a N
         * @param patt Pattern P to align
         * @param pattBegin Begin position of the pattern P
         * @param pattEnd End position of the pattern P
         * @param nodeBegin Begin position of the node N
         * @param nodeEnd End position of the node N
         * @return Node alignment object
         */
        NodeAlignment nodeAln(const string& patt,
                                size_t pattBegin, size_t pattEnd,
                                const NodePosPair& nodeBegin,
                                const NodePosPair& nodeEnd);

public:
        /**
         * Default constructor
         * @param dbg de Bruijn graph structure
         * @param readCorrDFSNodeLimit max number of visited nodes before give up in dfs search
         */
        GraphAligner(const DBGraph& dbg, MarkovChainHandler& mch_, size_t readCorrDFSNodeLimit) : dbg(dbg), maxVisits(readCorrDFSNodeLimit),
                aligner(2, 1, -1, -3),mch(mch_) {}

        /**
         * Set the maximum number of visited nodes in the search procedure
         * @param maxVisits Maximum number of visited nodes
         */
        void setMaxVisits(size_t maxVisits) {
                this->maxVisits = maxVisits;
        }

        /**
         * Find the best graph alignment given a pattern and a start NPP
         * @param P Pattern to align
         * @param srcNPP Source node-position pair
         * @param graphAln Best graph alignment (= highest alignment score)
         * @param markovFilter using markov framework to filter somt nodes
         * @param branchAndBoundFileter using branch and bound pruning technique to reduce the search space 
         * @param seedNodeChain already visited nodes in seed
         * @return OPTIMAL or EXHAUSTED depending search result
         */

        SearchRes DFSAln_guided(const string& P, const NodePosPair& srcNPP,
                                 GraphAlignment& bestGA, bool markovFilter, bool branchAndBoundFileter, vector<NodeID>seedNodeChain);


        /**
         * Find the best graph alignment given a pattern and a start NPP
         * @param P Pattern to align
         * @param srcNPP Source node-position pair
         * @param graphAln Best graph alignment (= highest alignment score)
         * @return OPTIMAL or EXHAUSTED depending search result
         */
        SearchRes DFSAln(const std::string& P, const NodePosPair& srcNPP,
                         GraphAlignment& graphAln);


        /**
         * Find the best graph alignment given a pattern and a src/dst NPP
         * @param P Pattern to align
         * @param srcBegin Source node-position pair
         * @param dstEnd Destination node-position pair
         * @param lenToDst Length from a node to the destination point
         * @param bestGA Best graph alignment (input/output)
         * @return OPTIMAL or EXHAUSTED depending on search result
         */
        SearchRes DFSAlnTo(const string& P, const NodePosPair& srcBegin,
                           const NodePosPair& dstEnd,
                           const NodeProp<size_t>& lenToDst,
                           GraphAlignment& bestGA);

        /**
         * Given dstEnd, compute the shortest path smaller than maxLen to
         * dstEnd. By definition, the shortest path from dst to dstEnd is the
         * shortest loop that first contains the entire dst node.
         * @param dstEnd Desination node-position pair
         * @param maxLen Maximum length to consider
         * @param lenToDst Shortest path to dstEnd (initialized to max)
         * @param modNodes List of nodes that were modified
         */
        void findShortestPathToNPP(const NodePosPair& dstEnd, size_t maxLen,
                                   NodeProp<size_t>& lenToDst,
                                   std::vector<NodeID>& modNodes) const;

        /**
         * Depth-first search for a parallel path that is parallel to and
         * most similar to the sequence in node with nodeID identifier
         * @param nodeID node identifier
         * @param graphAln Alignment result
         */
        void findParallelPath(NodeID nodeID, GraphAlignment& graphAln);
};

#endif
