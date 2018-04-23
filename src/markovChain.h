
#include "settings.h"
#include "graph.h"
#include <mutex>

#include "library.h"
//the fixed size of the maximum order of markov chain
const size_t maxMCOrder = 10;
const size_t minMCOrder = 2;
const float illuminaErrorRate =.005 ;
const size_t minimumExpectedCov = 10;
const size_t minOddratioLog = 5;
const float errorMean = 1;

const double inf = std::numeric_limits<double>::infinity();
class KmerSeed
{
public:
        std::vector<NodeID> nodeID;     // chain of node IDs
        size_t nodeFirst;               // offset within first node
        size_t readFirst;               // start position in read
        size_t readEnd;                 // end position in read


        KmerSeed() : nodeID(0), nodeFirst(0), readFirst(0), readEnd(0) {}

        KmerSeed(NodeID nodeID_, size_t nodeFirst_, size_t readFirst_, size_t readEnd_) :
             nodeFirst(nodeFirst_), readFirst(readFirst_), readEnd(readEnd_) {
                nodeID.push_back(nodeID_);
        }

        bool operator< (const KmerSeed& rhs) const {
                if (nodeID.front() != rhs.nodeID.front())
                        return nodeID.front() < rhs.nodeID.front();
                return nodeFirst < rhs.nodeFirst;
        }

        static void mergeSeeds(const std::vector<KmerSeed>& seeds,
                               std::vector<KmerSeed>& mergedSeeds);


};


class ParGraphMCH
{
private:
        size_t numNodes;                // number of nodes in the graph
        size_t targetChunkSize;         // target chunk size
        size_t currOffset;              // current node offset

        std::mutex mutex;               // mutex

public:
        /**
         * Get a chunk of nodes (thread-safe)
         * @param chunkOffset First nodeID to handle
         * @param chunkSize Number of nodes to handle
         */
        void getNodeChunk(size_t& chunkOffset, size_t &chunkSize);

        /**
         * Default constructor
         * @param numNodes_ Number of nodes in the graph
         * @param targetChunkSize_ Target size per chunk
         */
        ParGraphMCH(size_t numNodes_, size_t targetChunkSize_) :
                numNodes(numNodes_), targetChunkSize(targetChunkSize_),
                currOffset(1)  {}
};
//                              Decision
// ============================================================================
//      each chain of nodes, has a vector of Decision objects.
//      Decision object contains the next node at the end of that chain.
// ============================================================================

class Decision{
        NodeID nextNode;                        // Next node at the end of the chain
        size_t frequency;                       // Frequency of visiting nextNode at the end of the chain
        bool accepted;                          // Its true if this node is the winner among all possible next nodes.
        double oddRatioOfAcceptance;            // The odd ratio of correctness of choosing this node as the winner among all other next nodes.


public :
        /**
         * Constructor
         * @param nodeID_ next node ID
         * @param freq_ the frequency of visiting this node.
         * by default this node is not the winner node and the oddratio set to zero
         */
        Decision(NodeID nodeID_, size_t freq_):accepted(true), oddRatioOfAcceptance(0){
                nextNode = nodeID_;
                frequency = freq_;
        }


        /**
         * Constructor for already decided one
         * @param nodeID_ next node ID
         * @param freq_ the frequency of visiting this node.
         * by default this node is not the winner node and the oddratio set to zero
         */
        Decision(NodeID nodeID_, size_t freq_, bool accepted_, double oddratio_):nextNode(nodeID_),
        frequency(freq_),  accepted(accepted_),oddRatioOfAcceptance(oddratio_){}

        /**
         * Defult constructor
         *
         * */
        Decision(){};
        /**
         * returns the id of the next node
         */
        NodeID getNodeID(){
                return nextNode;
        }
        /**
         * returns the frequency of the next node
         * @return the frequency value
         */
        size_t getFrequency(){
                return frequency;
        }
        /**
         * set the frequency of the next node
         * @param value the value of frequency variable
         */
        void setFrequency( size_t value){
                frequency = value;

        }
        /**
         * set accepted variable to true if the node is the winner, and false if its not
         * @param value the value of accepted variable
         */
        void setDecision(bool value)
        {
                accepted = value;
        }
        bool getDecision(){
                return accepted;
        }
        double getOddRatioOfAcceptance(){
                return oddRatioOfAcceptance;
        }
        void setOddRatioOfAcceptance(double value){
                oddRatioOfAcceptance = value;
        }
};

class ChainRow{
private:

        size_t totalSampleSize ;         // This keeps the total number of smaples from the same order and same condition
        float length_avg;                // the avg length of seed
   
        
  

private :
        size_t order;                    // order of this row
        size_t multiplicity;             // the multiplicity of the chain
        bool singleOutgoint;
public :
        vector<NodeID> nodeChain;        // The list of observed node,  leftmost node is pushed into 0

        vector<Decision> decisionVector; // contains all the decisions which can be made with this chain of node.
        /**
         * Default constructor
         * 
         */
        ChainRow(){};
        /**
         * constructor
         * @param nodeChain_ The visited node chain
         * @param freq_ the frequency of visiting this chain, defult valu is 1
         */
        ChainRow(std::vector< NodeID > nodeChain_, size_t len, size_t freq);

        /**
         * Constructor for already decided rows
         **/
        ChainRow(vector<NodeID> nodeChain_,size_t length_,
                 size_t sampleSize ,size_t multiplicity,  vector<Decision> decisionVector_);


        void updateAVGlength(float len, size_t freq_ );

        /**
         * for all the possible decisions, if the coverage is appropriate, then calculates the updateOddRatio then make decision
         * @return returns true if any decision made
         **/
        bool makeDecision();
        /**
         * @param expectedCov expected num of read to support based on the multiplicities and length
         * @param mul multiplicity of the node
         */
        bool makeDecision( float expectedCov);


        void increaseFreqOfNode(NodeID nodeID, size_t len, size_t value);
        /**
         * compare two chins of nodes and returns true if they are exactly similar
         * @param nodeList chain of node that thould be compared
         */
        bool equal (vector<NodeID> nodeList);
        /**
         * compare two chins of nodes and returns true if the smaller chain is simlar to the end of bigger one, so they don't need to have same size
         * @param nodeList chain of node that thould be compared
         */

        bool equalEnd(vector<NodeID> nodeList);

        /**
         * print this chain
         */
        size_t getTotalSampleSize(){
                return totalSampleSize;
        }
        void setTotalSampleSize(size_t value){
                totalSampleSize = value;
        }
   
        bool getSingleValidOut (){
                return singleOutgoint;
        }
        void setMultiplicity(size_t _mul){
                multiplicity = _mul;
        }
        size_t getMultiplicity(){
                return multiplicity;
        }


        /**
         * returns the order of this chain.
         *
         * */
        size_t getOrder(){
                return order;
        }
        void printRow ();
        size_t getLength(){
                return length_avg;
        }
        void setLength(size_t value)
        {
                length_avg = value;
        }
      
        struct sortDecisions {
                bool operator()( Decision left, Decision right) {

                        if ( left.getFrequency() == right.getFrequency())
                                return (left.getNodeID() < right.getNodeID());
                        else
                                return (left.getFrequency() <right.getFrequency());
                }
        };
};


// ============================================================================
// contains different chains with the same orde at a specific current ndoe
// ============================================================================

class TransitionTable{
private:

        size_t chainLength;
        size_t order ;
public :
        vector <ChainRow> rows;                                                 // transition rows in this order
        bool modified;                                                          // if table has any row the value is true.
        /**
         * default constructor
         */
        TransitionTable():modified(false){}
        /**
         * constructor
         * @param order order of this table (markov chain order)
         */
        TransitionTable(size_t order_):order(order_), modified(false){}
        /**
         * adding a new row to this table
         * @param chainList list of visited nodes in the seed
         * @param len size of the seed
         */
        void addRowSlow(vector<NodeID> chainList, size_t len);
        /**
         * merging tables from different threads
         * @param subTable table need to be merged
         */
        void mergeTables(TransitionTable &subTable);
        /**
         * printing the table information in a formatted way
         *
         */
        void printRows();
        vector <ChainRow> getRows(){
                return rows;
        }

        bool smartSearchRow(vector<NodeID> chainList, ChainRow& row);

        /**
         * removes the chain rows which ends exactly as the given nodeChain
         *
         */
        void forceRemove(vector<NodeID> nodeChain);
        /** 
         * if a decision is made in lower order , it updates the higher orders
         */
        void updateDecisions (ChainRow row);
        /**
         * sort the rows in this table
         **/
        void sort();
        void addNewRow(ChainRow row);
        struct sortRows {
                bool operator()( ChainRow left, ChainRow right) {
                        if (left.nodeChain.size() == right.nodeChain.size()){
                                size_t i = 0;
                                while ( i < left.nodeChain.size()){
                                        if (left.nodeChain[i] != right.nodeChain[i])
                                                return (left.nodeChain[i] < right.nodeChain[i]);
                                        else
                                                i++;
                                }
                        }
                        else
                                return (left.nodeChain.size() < right.nodeChain.size());
                        return true;

                }

        };



};
// ============================================================================
// contains all the transition tables from different order
// ============================================================================

class TableContainer{

private:

public:
        
        bool modified;                                         //if the table has any vlue this variable is true (useful in merging and printing table)
        
        TransitionTable orderTable [maxMCOrder-minMCOrder+1];  //contains all the subTables of transitions from this node in different orders
        /**
         * Default constructor
         * */
        TableContainer ():modified(false){
        }
         void addNewChain(vector<NodeID> chainList, size_t len);
         void addNewRow(ChainRow &row, size_t order);
        /**
         * print the table information
         **/
        void printTable();

        /**
         * merge Table containers after all threads finish
         * @param subTables tables which needs to be merged for this node
         */
        void mergeTables(TableContainer &subTables);

        /**
         * if any decision has been chosen in a lower order, it propagate it to higher order by removing them.
         *
         */
        void propagateDecision( ChainRow row);
        /**
         * for each order, sorts the rows
         *
         */
        void sort();
        
    

};
// ============================================================================
// MARKOV CHAIN CLASS
// ============================================================================

class MarkovChain
{
private:
        const DBGraph &dbg;
        const Settings &settings;
        int numNodes;


public:
         size_t numReads;               // total number of reads in the data set, used to estimate coverage
         double avg_readLength;         // the avg length of reads in the data set, used to estimate coverage
         size_t maxReadSize;            // the maximum length of a read in the data set

         TableContainer  *markovTable;  // a pointer to an array conataining all the transitions from different nodes, in different orders.
                                        // to access the transitions from node nodeID, go to (numNodes + nodeID)th element.
private:
         /**
         * first finds the seeds then fill the transaction matrix
         * @param record Record to extract the transition
         */
        void extractTransitionMatrixFromSeed(const ReadRecord& record);
        
        
        
        /**
         * extract transaction matrix from subchainLists
         * 
         **/
        void extractTransactionFromSubchainLists (vector<vector<vector< NodeID>>> &subchainLists );
        
        void findSeedKmer(const std::string& read,
                          std::vector<KmerSeed>& seeds);
        /**
         * Correct the records in one chunk
         * @param npp Vector of node position pairs
         * @param seeds Vector containining consistent seeds (output)
         **/
        void extractSeeds(const std::vector<NodePosPair>& nppv,
                          std::vector<KmerSeed>& seeds);
        /**
         * Find the node position pairs for a read
         * @param read Reference to the read
         * @param npp  Vector of node position pairs
         */
        void findNPPFast(const std::string& read, std::vector<NodePosPair>& npp);
        /**
         * Extract all the subchains from length =2 to maxMCOrder from the seed and the reverse direction
         * @param nodeIDs  Vector containining consistent seeds (input)
         * @return a list conataining all the subchains of the seed.
         **/
        vector< vector<vector< NodeID >> > extractSubchains(const vector< NodeID> nodeIDs );
        /**
         * return the minimum marginal Length of a node chain.
         * @param nodeChain  list of nodes
         * @return length
         **/
        size_t getChainLength(vector<NodeID> nodeChain);


        void updateStatistic(size_t readLength);
public:

        /**
         * extract transitions of the records in one chunk
         * @param readChunk Chunk of records to correct
         */
        void extractTransitionMatrixInChunk(std::vector<ReadRecord>& readChunk);
        /**
         * fill the transitions table from chains on nodechain
         * @param chainNodes chan of nodechains
         */
        void fillTransitionMatrixFromList(vector < vector< NodeID >> &chainNodes);
        
        /**
         * Default constructor
         * @param dbg_ Reference to the De Bruijn graph
         * @param settings_ Reference to the settings class
         */
         MarkovChain(const DBGraph& dbg_, const Settings& settings_, const size_t maxMCOrder_ ) :
                          dbg(dbg_), settings(settings_), numNodes(dbg_.getNumNodes()), numReads(0),avg_readLength(0), maxReadSize(0){
                                  markovTable = new TableContainer [numNodes * 2 +1];

                        }

         /**
         * Destructor
         */
        ~MarkovChain() {
                delete [] markovTable;
                markovTable = NULL;
        }


};

// ============================================================================
// MARKOV CHAIN  HANDLER CLASS
// ============================================================================



class MarkovChainHandler
{
private:


        DBGraph &dbg;
        const Settings &settings;
        void workerThread(size_t myID, LibraryContainer& libraries);
        void aggregateResult(MarkovChain &markovChain);
        std::mutex mergeMutex;         // mutex for merging tables
        std::mutex searchMutex;
        int  numNodes;
        TableContainer  *markovTable;
        size_t estiGenomeSize;
        size_t estiCoverage;
        size_t numOfReads;
        double avg_readLength;
        size_t maxReadSize;
        string filename;

        double meanInsertSize;
        double stdInsertSize;
        vector<float> expectedCoverageByLen;    // keeps the expected coverage of seeds with different size.
        /**
         * count number of kmers in the graph
         * @return number of kmers in the graph as the size
         * */
        size_t getGraphSize();


        /**
         * based on the number of reads, size of graph and read size estimate the coverage.
         **/
        void updateStatistic();
        /**
         * upload some statistical info about the read file from disk
         * @parm fileName , name of the file that those static are stored.
         **/
        void uploadStatisticFromFile(string fileName);
        /**
         *search specific row in the table
         *
         */ 
        bool searchInMCTable(vector<NodeID> nodeChain,  ChainRow &row);
        /**
         * this function produces some random query from the table to see if the search function works properly.
         */
        void testSearchFunction();


        const unsigned char delimiter = '\t';

        enum command {ESTGENOMESIZE, NUMOFNODES, NUMOFARCS,NODEID,ORDER,NEWROW, NEWDECISION,NUMOFRESOLVEDREPEAT, ESTCOVERAGE , EXPECTEDCOVERAGE, MAXREADSIZE};
        std::map< command, string > IOCommand = {
                {NUMOFNODES, ">numberNodes"},
                {NUMOFARCS, ">numberArcs"},
                {NODEID, ">NODE"},
                {ORDER, ">order"},
                {NEWROW, ">newRow"},
                {NEWDECISION, ">newDecision"},
                {NUMOFRESOLVEDREPEAT, ">numOfSolvedRepeat"},
                {ESTCOVERAGE, ">estiCoverage"},
                {ESTGENOMESIZE, ">estiGenomeSize"},
                {EXPECTEDCOVERAGE,">expectedCovBylen"},
                {MAXREADSIZE,">maxReadSize"}
        };



        float getExpectedCovByLen(size_t len);
        
        /**
         * alternative solution to find the chain of nodes, instead of findign the seeds.
         * @param chainChunk output contains all chains of chain longer than minMCOrder
         * @param filename ncf file name which contan tha path which reads align to.
         */ 
        void readNodeChainFile ( vector < vector< NodeID >>& chainChunk  ,string fileName );

           
        /**
         * estimate the valuse of insertSizeMean and insertSizeSTD
         * @param inserSizeVe a vector contains all the insert size of all the mapped paired reads (both pairs map on the sam node)
         **/
        void estimateInsertSize(const std::vector<size_t> inserSizeVe);

        
public:
        
        /**
         * based on a given correct chain which is extracted from perfect reads, 
         * it checks if Markov chain mistakenly prohibits a true path or not.
         * 
         **/
        void verifyChains(string fileName);
        /**
         * splits the string of multiple words into a vector of strings
         * @param str the given string
         * @param delimiter the seperator charachter
         */
        vector<string> splitString(const string str,const char delimiter);
        /**
         * writes the final trnasactions from the nodes into the filename
         **/
        void writeInFile();
        /**
         * load from the file the trnasactions instead of calculate it from the scrach.
         */
        bool loadFromFile();
         /**
         * check if the file exists in the local disk to avoid populating it from the scrach
         * @return true if it exists
         **/
        bool markovTableFileExist();


        void printTable();
        /**
         * Default constructor
         */
        MarkovChainHandler(DBGraph& g, const Settings& s);

        /**
         * Destructor
         */
        ~MarkovChainHandler(){
                delete [] markovTable;
                markovTable = NULL;
        }

        /**
         *
         * @param libraries Library container with libraries
         */
        void extractTransitionMatrix(LibraryContainer &libraries);


        /**
         * it emiminates the wrong path according to the observed chain of nodes
         * @param nodeChain a vector contains chian of nodes.
         *  @return set<NodeID> the eligible nodes set in the right of last node in the chain
         */
        set<NodeID>  getPotentialNextNodes(vector<NodeID> nodeChain);

        /**
         * alternative way to train the markov model instead of seedng.
         * for each input read file , there is a ncf file contains the path which those read align.
         * this procedure reads that file to train the model.
         */
        void extractTransitionMatrixFromFile(LibraryContainer& libraries);

};

/**
 * ============================================================================
        Repeat Resolver CLASS, based on the probabilistic model it resolves the repeats
 * ============================================================================
**/

class RepeatResolver{
        TableContainer * markovTable;           // Table contains all the information about the transitions of nodes in the graph
        int numNodes;                           // Number of nodes in the grap, in order to access the markovTable
        size_t estiCoverage;                    // The estimated coverage of the data set
        double avg_readLength;                  // The average of read length in the input data set
        size_t maxReadSize;                     // The maximum size of a read in the input data set
        vector<float> expectedCoverageByLen;    // keeps the expected coverage of seeds with different size.
        const Settings &settings;               // the graph and libraries settings and parameters
        size_t minimumSampleSize ;              // the minimum amount of reads that should cover a chain 


       /**
        * it once calculate and keeps the expected coverage of seeds with different size.
        * @parma errorFreeAssumtion it assums that the error free coverage is needed (from seeds)
        *
        **/
        void updateExpectedCovByLen(bool errorFreeAssumtion);

        /**
         * returns the thread work size
         * 
         **/
        size_t getThreadWorkSize() const {
                return 32768;
        }
      
        /**
         * it solves the Repeat in a lower order if that would be possible
         * @param table a table contains all the subTables in different orders
         **/
        void makeDecision(TableContainer &tables);


public :

        /**
         * returns the expected coverage of a seed with a given length
         * @param len length of a seed
         *
         * @return the expected coverage
         **/

        float getExpectedCovByLen(size_t len);
        /**
         * Default constructor
         *
         **/
        RepeatResolver(TableContainer *markovTable_, int numNodes_,size_t estiCoverage_  , size_t avg_readLength_ ,size_t maxReadSize_ ,Settings settings_):markovTable(markovTable_), numNodes(numNodes_),
        estiCoverage(estiCoverage_), avg_readLength(avg_readLength_), maxReadSize(maxReadSize_), settings(settings_), minimumSampleSize(10){
        }
        /**
         * loops over nodes in the grah and resolves the repeats
         *
         **/

        void fitModelIntoTable(TableContainer &tables);
        /**
         *  It finds the repeats which can be solved for nodes between given interval (useful for parallel implementation )
         **/
        void fitModelIntoTableThread(size_t threadID, ParGraph& wlb);
        /**
         * @param errorFreeAssumtion it assumes that the coverage is calculated based on error free parts of reads (seeds)
         */
        void fitModelIntoTable(bool errorFreeAssumtion);
};



