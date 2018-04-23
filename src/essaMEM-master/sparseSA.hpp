#ifndef __sparseSA_hpp__
#define __sparseSA_hpp__

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>
#include <limits.h>


using namespace std;

static const unsigned int BITADD[256] =
{
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //0-9
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //10-19
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //20-29
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //30-39
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //40-49
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //50-59
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 0,        UINT_MAX, 1,        UINT_MAX, UINT_MAX,        //60-69                65:A        67:C
        UINT_MAX, 2,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //70-79                71:G
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 3,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //80-89                84:T
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 0,        UINT_MAX, 1,                //90-99                97:a        99: c
        UINT_MAX, UINT_MAX, UINT_MAX, 2,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //100-109        103:g
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 3,        UINT_MAX, UINT_MAX, UINT_MAX,        //110-119        116:t
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //120-129
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //130-139
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //140-149
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //150-159
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //160-169
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //170-179
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //180-189
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //190-199
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //200-209
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //210-219
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //220-229
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //230-239
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,        //240-249
        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX                                                 //250-255
};

// Stores the LCP array in an unsigned char (0-255). Values larger
// than or equal to 255 are stored in a sorted array.
// Simulates a vector<int> LCP;
struct vec_uchar {
        struct item_t{
                item_t(){}
                item_t(size_t i, int v) { idx = i; val = v; }
                size_t idx; int val;
                bool operator < (item_t t) const { return idx < t.idx; }
        };
        vector<unsigned char> vec; // LCP values from 0-65534
        vector<item_t> M;
        void resize(size_t N) { vec.resize(N); }
        // Vector X[i] notation to get LCP values.
        int operator[] (size_t idx) const {
                if(vec[idx] == numeric_limits<unsigned char>::max())
                        return lower_bound(M.begin(), M.end(), item_t(idx,0))->val;
                else
                        return vec[idx];
        }
        // Actually set LCP values, distingushes large and small LCP
        // values.
        void set(size_t idx, int v) {
                if(v >= numeric_limits<unsigned char>::max()) {
                        vec.at(idx) = numeric_limits<unsigned char>::max();
                        M.push_back(item_t(idx, v));
                }
                else { vec.at(idx) = (unsigned char)v; }
        }
        // Once all the values are set, call init. This will assure the
        // values >= 255 are sorted by index for fast retrieval.
        void init() { sort(M.begin(), M.end()); /*cout << "M.size()=" << M.size() << endl;*/ std::vector<item_t>(M).swap(M);}

        long index_size_in_bytes() const {
                long indexSize = 0L;
                indexSize += sizeof(vec) + vec.capacity()*sizeof(unsigned char);
                indexSize += sizeof(M) + M.capacity()*(sizeof(size_t)+sizeof(int));
                return indexSize;
        }
};

// Match find by findMEM.
struct match_t {
        match_t() { ref = 0; query = 0, len = 0; }
        match_t(long r, long q, long l) { ref = r; query = q; len = l; }
        long ref; // position in reference sequence
        long query; // position in query
        long len; // length of match
};

struct saTuple_t {
        saTuple_t(): left(0), right(0) {}
        saTuple_t(unsigned int l, unsigned int r): left(l), right(r){}
        unsigned int left;
        unsigned int right;
};

// depth : [start...end]
struct interval_t {
        interval_t() { start = 1; end = 0; depth = -1; }
        interval_t(long s, long e, long d) { start = s; end = e; depth = d; }
        void reset(long e) { start = 0; end = e; depth = 0; }
        long depth, start, end;
        long size() { return end - start + 1; }
};

struct sparseSA {
        vector<string> const &descr; // Descriptions of concatenated sequences.
        vector<long> &startpos; // Lengths of concatenated sequences.
        long maxdescrlen; // Maximum length of the sequence description, used for formatting.
        bool _4column; // Use 4 column output format.

        long N; //!< Length of the sequence.
        long logN; // ceil(log(N))
        long NKm1; // N/K - 1
        string &S; //!< Reference to sequence data.
        vector<unsigned int> SA; // Suffix array.
        vector<int> ISA; // Inverse suffix array.
        vec_uchar LCP; // Simulates a vector<int> LCP.
        vector<int> CHILD; //child table
        vector<saTuple_t> KMR;

        long K; // suffix sampling, K = 1 every suffix, K = 2 every other suffix, K = 3, every 3rd sffix
        bool hasChild;
        bool hasSufLink;
        //fields for lookup table of sa intervals to a certain small depth
        bool hasKmer;
        long kMerSize;
        long kMerTableSize;
        int sparseMult;
        bool printSubstring;
        bool printRevCompForw;
        bool forward;
        bool nucleotidesOnly;

        long index_size_in_bytes() {
                long indexSize = 0L;
                indexSize += sizeof(forward);
                indexSize += sizeof(printRevCompForw);
                indexSize += sizeof(printSubstring);
                indexSize += sizeof(sparseMult);
                indexSize += sizeof(hasSufLink);
                indexSize += sizeof(hasChild);
                indexSize += sizeof(K);
                indexSize += sizeof(NKm1);
                indexSize += sizeof(logN);
                indexSize += sizeof(N);
                indexSize += sizeof(_4column);
                indexSize += sizeof(maxdescrlen);
                indexSize += sizeof(descr);
                indexSize += sizeof(hasKmer);
                indexSize += sizeof(kMerSize);
                indexSize += sizeof(kMerTableSize);
                indexSize += sizeof(nucleotidesOnly);
                for(size_t i = 0; i < descr.size(); i++){
                        indexSize += descr[i].capacity();
                }
                indexSize += sizeof(startpos) + startpos.capacity()*sizeof(long);
                indexSize += S.capacity();
                indexSize += sizeof(SA) + SA.capacity()*sizeof(unsigned int);
                indexSize += sizeof(ISA) + ISA.capacity()*sizeof(int);
                indexSize += sizeof(CHILD) + CHILD.capacity()*sizeof(int);
                indexSize += sizeof(KMR) + KMR.capacity()*(2*sizeof(unsigned int));
                indexSize += LCP.index_size_in_bytes();
                return indexSize;
        }

        // Maps a hit in the concatenated sequence set to a position in that sequence.
        void from_set(long hit, long &seq, long &seqpos) const {
                // Use binary search to locate index of sequence and position
                // within sequence.
                vector<long>::iterator it = upper_bound(startpos.begin(), startpos.end(), hit); // SG: should use vector<long>::const_iterator
                seq = distance(startpos.begin(), it) - 1;
                it--;
                seqpos = hit - *it;
        }

        // Constructor builds sparse suffix array.
        sparseSA(string &S_, vector<string> const &descr_, vector<long> &startpos_,
        bool __4column, long K_, bool suflink_, bool child_, bool kmer_, int sparseMult_,
        int kMerSize_, bool printSubstring_, bool printRevCompForw_, bool nucleotidesOnly_);

        // Modified Kasai et all for LCP computation.
        void computeLCP();
        //Modified Abouelhoda et all for CHILD Computation.
        void computeChild();
        //build look-up table for sa intervals of kmers up to some depth
        void computeKmer();

        // Radix sort required to construct transformed text for sparse SA construction.
        void radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h);

        // Prints match to cout.
        void print_match(match_t m) const;
        void print_match(match_t m, vector<match_t> &buf) const; // buffered version
        void print_match(string meta, vector<match_t> &buf, bool rc) const; // buffered version

        //Check if the matches are correct
        void checkMatches(std::string const &P, std::vector<match_t> const &matches, int const min_len) const;

        //find the first l index of an lcp-interval
        int get_first_l(int const start, int const end) const;

        // Binary search for left boundry of interval.
        inline long bsearch_left(char c, long i, long s, long e) const;
        // Binary search for right boundry of interval.
        inline long bsearch_right(char c, long i, long s, long e) const;

        // Simple suffix array search.
        inline bool search(string const &P, long &start, long &end) const;

        // Simple top down traversal of a suffix array.
        inline bool top_down(char c, long i, long &start, long &end) const;
        inline bool top_down_faster(char c, long i, long &start, long &end) const;
        inline bool top_down_child(char c, interval_t &cur) const;

        // Traverse pattern P starting from a given prefix and interval
        // until mismatch or min_len characters reached.
        inline void traverse(string const &P, long prefix, interval_t &cur, int min_len) const;
        inline void traverse_faster(const string &P, const long prefix, interval_t &cur, int min_len) const;

        // Simulate a suffix link.
        inline bool suffixlink(interval_t &m) const;

        // Expand ISA/LCP interval. Used to simulate suffix links.
        inline bool expand_link(interval_t &link) const {
                long thresh = 2 * link.depth * logN, exp = 0; // Threshold link expansion.
                long start = link.start;
                long end = link.end;
                while(LCP[start] >= link.depth) {
                        exp++;
                        if(exp >= thresh) return false;
                        start--;
                }
                while(end < NKm1 && LCP[end+1] >= link.depth) {
                        exp++;
                        if(exp >= thresh) return false;
                        end++;
                }
                link.start = start; link.end = end;
                return true;
        }

        // Given a position i in S, finds a left maximal match of minimum
        // length within K steps.
        inline void find_Lmaximal(string const &P, long prefix, long i, long len, vector<match_t> &matches, int min_len, bool print) const;

        // Given an interval where the given prefix is matched up to a
        // mismatch, find all MEMs up to a minimum match depth.
        void collectMEMs(string const &P, long prefix, interval_t mli, interval_t xmi, vector<match_t> &matches, int min_len, bool print) const;

        // Find all MEMs given a prefix pattern offset k.
        void findMEM(long k, string const &P, vector<match_t> &matches, int min_len, bool print) const;

        // NOTE: min_len must be > 1
        void findMAM(string const &P, vector<match_t> &matches, int min_len, long& memCount, bool print) const;
        inline bool is_leftmaximal(string const &P, long p1, long p2) const;

        // Maximal Almost-Unique Match (MAM). Match is unique in the indexed
        // sequence S. as computed by MUMmer version 2 by Salzberg
        // et. al. Note this is a "one-sided" query. It "streams" the query
        // P throught he index. Consequently, repeats can occur in the
        // pattern P.
        void MAM(string const &P, vector<match_t> &matches, int min_len, long& memCount, bool forward_, bool print) {
                forward = forward_;
                if(K != 1) return; // Only valid for full suffix array.
                findMAM(P, matches, min_len, memCount, print);
        }

        // Find Maximal Exact Matches (MEMs)
        void MEM(string &P, vector<match_t> &matches, int min_len, bool print, long& memCount, bool forward_, int num_threads = 1);

        // Maximal Unique Match (MUM)
        void MUM(string const &P, vector<match_t> &unique, int min_len, long& memCount, bool forward_, bool print);

        //save index to files
        void save(const string &prefix);

        //load index from file
        bool load(const string &prefix);

        //construct
        void construct();
};


#endif // __sparseSA_hpp__

