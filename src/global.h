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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <cassert>
#ifdef _MSC_VER
        #include "msstdint.h"     // uint64_t types, etc.
#else
        #include <stdint.h>
#endif
#include "util.h"       // timing routines
#include <thread>
#include <cmath>

// ============================================================================
// DEFINITIONS
// ============================================================================

#ifndef MAXKMERLENGTH
        #define MAXKMERLENGTH 63        // this is the maximum k-mer length possible
#endif

#define KMERBYTESIZE ((MAXKMERLENGTH+3)/4)
#define KMERBYTEREDUCTION 2
typedef uint32_t KmerLSB;       // set to uint16_t and all hell will break loose

#ifdef __GNUC__
        #define ATTRIBUTE_PACKED __attribute__ ((packed))
#else
        #define ATTRIBUTE_PACKED
#endif

#define MAXGAPS 3

#define MAX_COVERAGE 65535
#define MAX_MULTIPLICITY 255
#define OUTPUT_FREQUENCY 32768

#define MAX_PATH_SEARCHES 100000
#define MAX_PATH_SOLUTIONS 20

#define MAX_NUM_OBSERVATIONS 1000

#define ZERO_SCORE_WINDOW 5
#define MAX_THRESHOLD 1
#define MULT_SIGN_STD 3

#define FLOAT_SMALL 1E-6
#define DOUBLE_SMALL 1E-15

#define SCAFFOLD_SNAP_LINK 5

// The number of recordBlocks simulataneously held in memory. Higher values
// result in more parallel chunks at the cost of increased memory use.
#define NUM_RECORD_BLOCKS 2

// ============================================================================
// TYPEDEFS
// ============================================================================

typedef int32_t NodeID; // max 2 billion nodes
const int32_t MAX_NODE_ID = 2147000000;

typedef int32_t ArcID;  // max 2 billion arcs
typedef uint32_t NodeLength;    // max length is 4 billion
typedef NodeLength NodePosition;    // position in a contig or read
typedef uint64_t NucleotideID;
typedef uint32_t Coverage;       // coverage of an arc
typedef uint8_t Multiplicity;   // multiplicity of a node, arc, etc
typedef int32_t ReadID; // max 2 billion reads
typedef float Time;    // time, as defined by D.Z.
typedef int64_t ssize_t;
typedef int32_t sPositionID;

template<size_t numBytes>
class TKmer;

template<size_t numBytes>
struct TKmerHash;

typedef TKmer<KMERBYTESIZE> Kmer;      // full blown k-mer
typedef TKmer<KMERBYTESIZE - KMERBYTEREDUCTION> RKmer; // reduced k-mer
typedef TKmerHash<KMERBYTESIZE> KmerHash;
typedef TKmerHash<KMERBYTESIZE - KMERBYTEREDUCTION> RKmerHash;

typedef std::pair<NodeID, NodeID> NodePair;

std::ostream &operator<<(std::ostream &out, const NodePair& np);
typedef std::pair<NodeID, size_t> NodeSpec;

// ============================================================================
// CONSTANTS
// ============================================================================

static const int indelScore = 0;        // negative == penalty
static const int matchScore = 1;        // match score (> 0 !!)
static const int scoreMatrix[4][4] = {  // make sure the diagonal elements equal
        {matchScore, 0, 0, 0},
        {0, matchScore, 0, 0},
        {0, 0, matchScore, 0},
        {0, 0, 0, matchScore}
};

// ============================================================================
// ENUMS
// ============================================================================

enum MapResult { AmbiguousMap, NoMap, MultSolutions, MapToSame, MapToConnected,
                 UniqueMap, IncompleteMap };


#endif
