
#include <thread>
#include <string>
#include <iomanip>
#include <map>
#include "markovChain.h"


MarkovChainHandler::MarkovChainHandler(DBGraph& g, const Settings& s) :
dbg(g), settings(s), numNodes(dbg.getNumNodes()), estiGenomeSize(getGraphSize()),numOfReads(0),avg_readLength(0),maxReadSize(0), filename(settings.getTempDirectory()+"markovTable.dat"), meanInsertSize(0), stdInsertSize(0) {
        
        markovTable = new TableContainer[numNodes * 2 +1];
}

bool MarkovChainHandler::markovTableFileExist(){
        std::ifstream file(filename.c_str(), std::ios::in);
        bool OK = file.good();
        file.close();
        return OK;
}
bool MarkovChainHandler::loadFromFile()
{

        ifstream tableFile(filename.c_str());
        if (!tableFile){
                cout <<"Can't open " + filename;
                return false;
        }
        std::string line;
        NodeID currentNode = 0;
        size_t order = 0 ;
        size_t sampleSize =0;
        size_t length = 0;
        size_t multiplicity = 0;
        size_t numOfDecisions = 0;
        size_t decisionVectorSize = 0;
        vector<Decision> decisionVector;
        while (std::getline(tableFile, line))
        {
             
                TableContainer* tc;
                vector< string> tokens = splitString(line,delimiter);

                // all lines has at least two tokens
                if (tokens.size() <2)
                        return false;

                // check number of nodes to be same
                if  (tokens[0].compare(IOCommand[NUMOFNODES]) == 0) {
                        if (dbg.getNumNodes() != stoi( tokens[1])){
                                cout << "Mismatch between nodes" <<endl;
                                return false;
                        }
                        continue;

                }
                // check number of nodes to be same
                if  (tokens[0].compare(IOCommand[MAXREADSIZE]) == 0) {
                        maxReadSize = stoi( tokens[1]);
                        continue;

                }
                
                //check number of arcs to be the same
                if  (tokens[0].compare(IOCommand[NUMOFARCS]) == 0) {
                        if (dbg.getNumArcs() != stoi( tokens[1])){
                                cout << "Mismatch between arcs" <<endl;
                                return false;
                        }
                        continue;
                }
                // load estimated coverage
                if  (tokens[0].compare(IOCommand[ESTCOVERAGE]) == 0) {
                        estiCoverage = stoi( tokens[1]);
                        continue;

                }
                // load estimated coverage
                if  (tokens[0].compare(IOCommand[ESTGENOMESIZE]) == 0) {
                        estiGenomeSize = stoi( tokens[1]);
                        continue;

                }
                //load the expected coverage talbe 
                if  (tokens[0].compare(IOCommand[EXPECTEDCOVERAGE]) == 0) {
                        for ( size_t i = 1 ; i< tokens.size(); i++){
                                expectedCoverageByLen.push_back( stof( tokens[i]));
                        }
                        continue;
                }
                
                if  (tokens[0].compare(IOCommand[NODEID]) == 0) {
                        //if (currentNode !=0)
                        //    markovTable[currentNode + numNodes] = tc;
                        currentNode  = stoi( tokens[1]);
                        tc = &markovTable[currentNode + numNodes];
                        continue;
                }
                if  (tokens[0].compare(IOCommand[ORDER]) == 0) {
                        order  = stoi( tokens[1]);
                        continue;
                }
                if  (tokens[0].compare(IOCommand[NEWROW]) == 0) {
                        if (tokens.size() <4)
                                return false;
                        sampleSize  = stoi( tokens[1]);
                        length = stoi( tokens[2]);
                        decisionVectorSize = stoi(tokens[3]);
                        multiplicity = stoi(tokens[4]);
                        numOfDecisions = 0;
                        decisionVector.clear();;

                        continue;
                }
                // tokenName, nodeChain, nextNode, freq, accepcted, oddratio
                if  (tokens[0].compare(IOCommand[NEWDECISION]) == 0) {
                        if (tokens.size() <1 + order +1 +3)
                                return false;
                        vector<NodeID> nodeChain;
                        NodeID nextNode = stoi( tokens[order + 1]);
                        size_t frequency = stoi( tokens[order + 2]);;
                        bool accepcted = false;
                        if (stoi( tokens[order + 3])==1)
                                accepcted = true;
                        else if (stoi( tokens[order + 3])==0)
                                accepcted = false;
                        else
                                return false;
                        double oddRatio = 0;
                        if (tokens[ order + 4]=="inf")
                                oddRatio = std::numeric_limits<double>::infinity();
                        else
                                oddRatio = stoi(tokens[ order + 4]);

                        for (size_t  i = 1; i < order+1; i++)
                                nodeChain.push_back(stoi( tokens[ i]));
                        
                        Decision decision(nextNode, frequency,accepcted,oddRatio);
                        decisionVector.push_back(decision);
                        numOfDecisions ++;

                        if (numOfDecisions == decisionVectorSize){
                                ChainRow row(nodeChain,length,sampleSize , multiplicity, decisionVector);
                                tc->addNewRow(row, order);
                        }
                }
        }
        return true;
}

vector<string> MarkovChainHandler::splitString(const string str,const char delimiter)
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
void MarkovChainHandler::writeInFile()
{
        ofstream tableFile(filename.c_str());
        tableFile << IOCommand[NUMOFNODES]<<delimiter <<dbg.getNumNodes()<<endl;
        tableFile << IOCommand[NUMOFARCS]<<delimiter<< dbg.getNumArcs()<<endl;
        tableFile << IOCommand[ESTCOVERAGE]<<delimiter<< estiCoverage <<endl;
        tableFile << IOCommand[ESTGENOMESIZE]<<delimiter<< estiGenomeSize <<endl;
        tableFile << IOCommand[MAXREADSIZE]<<delimiter<< maxReadSize <<endl;
        size_t numOfSolvedRepeat = 0;
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id ==0)
                        continue;
                TableContainer tc = markovTable[ id + numNodes];
                if (!tc.modified)
                        continue;
                tableFile << IOCommand[NODEID]<<delimiter<< id << endl;

                for (size_t order = minMCOrder; order< maxMCOrder ; order++){
                        if (!tc.orderTable[order-minMCOrder].modified)
                                continue;
                        tableFile << IOCommand[ORDER] << delimiter<<order << endl;
                        for (size_t i =0; i <tc.orderTable[order-minMCOrder].rows.size(); i++){
                                ChainRow& row = tc.orderTable[order-minMCOrder].rows[i];
                                tableFile <<IOCommand[NEWROW] <<delimiter  <<row.getTotalSampleSize() <<delimiter<< row.getLength()
                                                              << delimiter <<row.decisionVector.size()<< delimiter<< row.getMultiplicity()<< endl;
                                for ( size_t i=0; i<row.decisionVector.size(); i++)
                                {
                                        Decision nextNodes = row.decisionVector[i];

                                        tableFile <<IOCommand[NEWDECISION]<<delimiter;
                                        for (size_t j=0; j<row.nodeChain.size(); j++)
                                                tableFile <<row.nodeChain[j] <<delimiter;


                                        tableFile <<nextNodes.getNodeID() <<delimiter<<nextNodes.getFrequency() << delimiter <<nextNodes.getDecision()
                                                                          <<delimiter<<nextNodes.getOddRatioOfAcceptance() <<endl;

                                }
                                numOfSolvedRepeat++;
                        }

                }
        }
        tableFile << IOCommand[NUMOFRESOLVEDREPEAT]<< delimiter <<numOfSolvedRepeat<<endl;
        tableFile << IOCommand[EXPECTEDCOVERAGE] << delimiter;
        for (size_t i = 0; i < expectedCoverageByLen.size()-1; i++)
                tableFile << expectedCoverageByLen[i] << delimiter;
        tableFile << expectedCoverageByLen[expectedCoverageByLen.size()-1] << endl;
        
        tableFile.close();
        cout <<  numOfSolvedRepeat << " number of resolved repeats were written in file." <<endl;
}


void MarkovChainHandler::workerThread(size_t myID, LibraryContainer& libraries){

        MarkovChain markovChain(dbg, settings, maxMCOrder);
        // local storage of reads
        vector<ReadRecord> myReadBuf;

        bool result = true;
        while (result) {
                size_t blockID, recordID;
        
                result = libraries.getRecordChunk(myReadBuf, blockID, recordID);
                markovChain.extractTransitionMatrixInChunk(myReadBuf);
        }
        aggregateResult( markovChain);
}
void MarkovChainHandler::extractTransitionMatrix(LibraryContainer& libraries){

        Util::startChrono();
        cout << "Creating kmer lookup table... "; cout.flush();
        dbg.buildKmerNPPTable();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        Util::startChrono();

        const unsigned int& numThreads = settings.getNumThreads()> 4? 4:settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        libraries.startIOThreads(settings.getThreadWorkSize(),
                                 10 * settings.getThreadWorkSize() * settings.getNumThreads());

        // start worker threads

        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&MarkovChainHandler::workerThread,
                                          this, i, ref(libraries));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
        libraries.joinIOThreads();
        dbg.destroyKmerNPPTable();
        updateStatistic();
        //printTable();
        RepeatResolver repeatResolver (markovTable, numNodes, estiCoverage,avg_readLength, maxReadSize,settings);
        repeatResolver.fitModelIntoTable(true);
      
        for (size_t i = minMCOrder + Kmer::getK(); i <maxReadSize; i++){
               expectedCoverageByLen.push_back( repeatResolver.getExpectedCovByLen(i));
        }
}

void MarkovChainHandler::estimateInsertSize(const std::vector<size_t> inserSizeVe){
        // std::sort(inserSizeVe.begin(), inserSizeVe.end());
         int outlier = floor( inserSizeVe.size()/4);
         size_t num = 0;
         for (size_t i = outlier; i< inserSizeVe.size()- outlier; i++){
                        meanInsertSize = meanInsertSize + inserSizeVe[i];
                        num ++;
         }
         meanInsertSize = round( meanInsertSize /num );
  
         for (size_t i = outlier; i< inserSizeVe.size()- outlier; i++){
                        stdInsertSize = stdInsertSize + pow ((inserSizeVe[i] - meanInsertSize), 2);
                        
         }
         stdInsertSize = round (sqrt( stdInsertSize /num));
         cout << "Mena of insert size:" <<meanInsertSize<< endl;
         cout << "STD of insert size :" << stdInsertSize << endl;
         
}
void MarkovChainHandler::extractTransitionMatrixFromFile(LibraryContainer& libraries){
        vector < vector< NodeID >> chainNodes;
        MarkovChain markovChain(dbg, settings, maxMCOrder);
        uploadStatisticFromFile(settings.getTempDirectory()+"sta.txt");
        for (size_t i = 0; i < libraries.getSize(); i++) {
                const ReadLibrary &input = libraries.getInput(i);
                string fileName = input.getNodeChainFilename();
                readNodeChainFile ( chainNodes,  fileName);
        
        }
        markovChain.fillTransitionMatrixFromList(chainNodes);
        aggregateResult( markovChain );
       
        RepeatResolver repeatResolver (markovTable, numNodes, estiCoverage, avg_readLength, maxReadSize, settings);
        repeatResolver.fitModelIntoTable(false);
          
        for (size_t i = minMCOrder + Kmer::getK(); i <maxReadSize; i++){
               expectedCoverageByLen.push_back( repeatResolver.getExpectedCovByLen(i));
        }
        for (size_t i = minMCOrder + Kmer::getK(); i <maxReadSize; i++){
               expectedCoverageByLen.push_back( repeatResolver.getExpectedCovByLen(i));
        }

}
void MarkovChainHandler::readNodeChainFile ( vector < vector< NodeID >>& chainChunk, string fileName ){
        std::ifstream ifs(fileName.c_str());
        std::set<string> ids;
        vector<NodeID> firstNodeChain, secondNodeChain ;
        while(!ifs.eof())
        {
                string line = "";
                std::getline(ifs, line);
                if (ifs.eof()){
                        ifs.close();
                        break;
                }
                size_t first = line.find("C(") + 2;
                size_t last = line.find(")C");
                string chain = line.substr (first,last-first);
                vector<string> items = splitString(chain,' ');
             
                
                if (items.size()>= minMCOrder){
                        vector <NodeID> nodeIDs;
                        for (string it:items){
                                int id = std::atoi (it.c_str());
                                nodeIDs.push_back (id );
                        }
                        chainChunk.push_back(nodeIDs);
                }
            
        }
        ifs.close();
        
}
void MarkovChainHandler::testSearchFunction(){
        vector<NodeID> nodeList;
        ChainRow row;
        size_t numberOfTest  = 0;
        //look up  for the  randomly selected items of existing elements in the table
        while (numberOfTest <1000){
                size_t nodeid = rand() % (numNodes *2);
                if (nodeid ==0)
                        continue;
                TableContainer tc = markovTable[ nodeid ];
                if (!tc.modified)
                        continue;
                size_t order = rand() % ((maxMCOrder -minMCOrder)+1);
                TransitionTable tt = tc.orderTable[order-minMCOrder];
                if (tt.getRows().size() ==0)
                        continue;
                size_t rowNum = rand()% (tt.getRows().size());
                ChainRow row = tt.getRows()[rowNum];
                ChainRow result;

                if (!searchInMCTable(row.nodeChain, result)){
                        cout <<"Unfortunatly we couldn't find it in the existing items :(" <<endl;
                        tc.printTable();
                        cout << "Searching for this item :" <<endl;
                        row.printRow();
                        searchInMCTable(row.nodeChain, result);
                }

                numberOfTest ++;
        }
       numberOfTest = 0;
       // search for a non existence item that we make it randomly
       while (numberOfTest <1000){
                size_t nodeid = rand() % (numNodes *2);
                if (nodeid ==0)
                        continue;
                TableContainer tc = markovTable[ nodeid ];
                if (!tc.modified)
                        continue;
                size_t order = rand() % ((maxMCOrder -minMCOrder)+1);
                TransitionTable tt = tc.orderTable[order-minMCOrder];
                if (tt.getRows().size() ==0)
                        continue;
                size_t rowNum = rand()% (tt.getRows().size());
                ChainRow row = tt.getRows()[rowNum];
                ChainRow result;
                if (row.nodeChain.size() <2)
                        continue;
                for (size_t i = 0; i<row.nodeChain.size()-1; i++)
                {
                        int newId = rand() % (numNodes);
                        if (rand()%2 ==0)
                                newId = -newId ;
                        row.nodeChain[i] = newId;
                }
                if (searchInMCTable(row.nodeChain, result)){
                        cout <<"Unfortunatly we  find it in the existing items :(" <<endl;
                        tc.printTable();
                        cout << "Searching for this item :" <<endl;
                        row.printRow();
                        searchInMCTable(row.nodeChain, result);
                }

                numberOfTest ++;
        }
}
void MarkovChainHandler::updateStatistic(){


        estiCoverage = ( (double)numOfReads / (double) estiGenomeSize  ) * avg_readLength ;
        //cout       << "\nThe illumina Error Rate assumed to be : "             << illuminaErrorRate * 100 <<"%" <<endl;
        cout       << "Number of reads handled : "     << numOfReads         << endl;
        cout       << "The average read length is : "  <<fixed<< avg_readLength     << endl;
        cout       << "The estimated genome size is : " << estiGenomeSize     << endl;
        cout       << "The estimated coverage is : "   <<std::setprecision(1)<< std::fixed<< estiCoverage<<endl;
        cout       << "The minimum odd ratio is (log scale) : "   << minOddratioLog << endl;
}
void MarkovChainHandler::uploadStatisticFromFile(string fileName){
        std::ifstream ifs(fileName.c_str());
        std::set<string> ids;
        while(ifs.is_open() && !ifs.eof())
        {
                string line = "";
                std::getline(ifs, line);
                vector<string> items = splitString(line,':');
                if(items[0]==  "NumOfReads"){
                        numOfReads = atoi (items[1].c_str());
                }
                if(items[0]==  "AvgReadLen"){
                        avg_readLength = atoi (items[1].c_str());
                }
                if(items[0]==  "MaxReadLen"){
                        maxReadSize = atoi (items[1].c_str());
                }
                
        }
        ifs.close();
        
        estiCoverage = ( (double)numOfReads / (double) estiGenomeSize  ) * avg_readLength ;
        //cout << "\nThe illumina Error Rate assumed to be :" << illuminaErrorRate * 100 <<"%" <<endl;
        cout << "\nTotal number of reads : " << numOfReads << endl;
        cout << "The average read length is : " << avg_readLength <<endl;
        cout << "The estimated genome size is : " << estiGenomeSize << endl;
        cout << "The estimated coverage is : "  <<std::setprecision(1)<<std::fixed<< estiCoverage<<endl;  
        cout << "The minimum odd ratio is (log scale) : "   << minOddratioLog << endl;
}
bool MarkovChainHandler::searchInMCTable(vector<NodeID> nodeChain, ChainRow &row){
        NodeID currNode = nodeChain.back();
        if (currNode > numNodes || currNode< -numNodes || currNode ==0)
                return false;
        size_t order = nodeChain.size();
        while (order >= minMCOrder){
                TransitionTable transition = markovTable [currNode +numNodes].orderTable[order-minMCOrder];
                if (transition.modified && transition.smartSearchRow( nodeChain , row))
                        return true;
                nodeChain.erase(nodeChain.begin());
                order = nodeChain.size();
        }

       return false;
}


set<NodeID>  MarkovChainHandler::getPotentialNextNodes(vector<NodeID> nodeChain){
        
        set<NodeID> results;
        SSNode node = dbg.getSSNode( nodeChain.back());
        // First suppose all next nodes are eligeble to extend
        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++){
                results.insert(it->getNodeID());
        }
        if (nodeChain.size() < minMCOrder || node.getNumRightArcs() ==1 ){
                return results;
        }

        while (nodeChain.size()>maxMCOrder){
                nodeChain.erase(nodeChain.begin());
        }
        ChainRow row;

        // the chain dosn't exist in the table, retun all outgoing arcs
        if (!searchInMCTable(nodeChain, row))
                return results;
        
        
        results.clear();
        
        // the chain exist in the table so return only allowed outgoing nodes.
        
        for (size_t i = 0 ; i < row.decisionVector.size (); i++){
                if (row.decisionVector[i].getDecision()){
                        results.insert(row.decisionVector[i].getNodeID());
                }
        } 
        return results;
}   


void MarkovChainHandler::verifyChains (string fileName){
        std::ifstream ifs(fileName.c_str());
        
        while(ifs.is_open() && !ifs.eof())
        {
                string line = "";
                std::getline(ifs, line);
                if (line =="")
                        continue;
                size_t first = line.find("C(") + 2;
                size_t last = line.find(")C");
                string chain = line.substr (first,last-first);
                vector<string> items = splitString(chain,' ');
                if (items.size()> minMCOrder){
                        vector <NodeID> chain;
                        for (string it:items){
                                int id = std::atoi (it.c_str());
                                chain.push_back (id );
                        }
                        size_t i = chain.size()-1;
                        while ( i >= minMCOrder ){
                                NodeID next = chain.back();
                                chain.pop_back();
                                set<NodeID >nextNodes = getPotentialNextNodes(chain);
                                if (nextNodes.find(next) == nextNodes.end()){
                                        cout << "A correct paht not found in your MM" <<endl;
                                        cout << "grep -e \"" ;
                                        for (auto it:chain)
                                                cout << it << " ";
                                        cout << next  <<"\" reads.ncf " <<endl;
                                        getPotentialNextNodes(chain);
                                       
                                }
                                i --;
                        }
                }
                
        }
        ifs.close();
        
}

float MarkovChainHandler::getExpectedCovByLen(size_t len){
        if (len < Kmer::getK()+ minMCOrder || len > maxReadSize)
                return 0;
        return (expectedCoverageByLen[len - Kmer::getK()-minMCOrder])> 0 ? (expectedCoverageByLen[len - Kmer::getK()-minMCOrder]) :0;

}

void MarkovChainHandler::printTable(){

        for (int  id =- numNodes; id <= numNodes; id++){
                if (id ==0)
                        continue;

                TableContainer tc = markovTable[ id + numNodes];
                if (tc.modified){
                        cout << "*************************************************************************************************" <<endl;
                        cout << "These tables show the transition probabilities to open the new nodes from "<<id <<" based on previous nodes" << endl;
                        tc.printTable();

                }
        }
}


void MarkovChainHandler::aggregateResult(MarkovChain& markovChain ){

       lock_guard<mutex> lock(mergeMutex);
       for (int id = -numNodes; id <= numNodes; id++) {
               if (id == 0)     // skip node 0, doesn't exist
                       continue;
               TableContainer subTables = markovChain.markovTable[id+numNodes];
               markovTable[id+numNodes].mergeTables(subTables);
       }
       avg_readLength = (double) ( avg_readLength * numOfReads + markovChain.numReads * markovChain.avg_readLength)/(double)(numOfReads+ markovChain.numReads);
       numOfReads = numOfReads + markovChain.numReads ;
       if (maxReadSize < markovChain.maxReadSize)
               maxReadSize = markovChain.maxReadSize;


}
size_t MarkovChainHandler ::getGraphSize(){

        size_t totMargLength = 0;       // total marginal length of all nodes

        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = dbg.getSSNode(id);
                if (!node.isValid())
                        continue;
                totMargLength += node.getMarginalLength();
        }
        return totMargLength;
}

void MarkovChain::extractTransitionMatrixInChunk(vector<ReadRecord>& readChunk)
{
        for (ReadRecord& it : readChunk){
                extractTransitionMatrixFromSeed(it);
                updateStatistic(it.read.length());

        }

}
void MarkovChain::fillTransitionMatrixFromList(vector < vector< NodeID >> &chainNodes)
{
        for (vector<NodeID> chain : chainNodes){
                vector<vector<vector< NodeID>>> subchainLists = extractSubchains(chain);
                extractTransactionFromSubchainLists(subchainLists);
        }     
        
}
void MarkovChain ::updateStatistic(size_t readLength){

        avg_readLength = (double)(avg_readLength * numReads + readLength)/(double) (numReads + 1) ;
        numReads ++;
        if (readLength > maxReadSize)
                maxReadSize = readLength;

}
void MarkovChain::extractTransitionMatrixFromSeed(const ReadRecord& record){

        string read = record.read;

        if (read.length() < Kmer::getK() + minMCOrder)
                return;
        vector<KmerSeed> seeds;
        findSeedKmer(read, seeds);
        for (KmerSeed seed :seeds){
                vector<vector<vector< NodeID>>> subchainLists = extractSubchains(seed.nodeID);
                extractTransactionFromSubchainLists(subchainLists);
        }
}
void MarkovChain::extractTransactionFromSubchainLists (vector<vector<vector< NodeID>>> &subchainLists ){
        for (size_t i =0 ; i< subchainLists.size(); i++){
                for (size_t j =0; j< subchainLists[i].size(); j++){
                        vector<NodeID> nodeChain = subchainLists[i][j];
                        vector<NodeID> reverseSubchainNode;
                        size_t len = getChainLength(nodeChain);
                        NodeID currentNodeID = nodeChain[nodeChain.size()-2];
                        SSNode currentNode = dbg.getSSNode(currentNodeID);
                        if (currentNode.getNumRightArcs() >1)
                                markovTable[currentNodeID + numNodes].addNewChain(nodeChain, len );
                        
                        //now reverse complement of the chain
                        reverseSubchainNode = nodeChain;
                        reverse(reverseSubchainNode.begin(),reverseSubchainNode.end());
                        for (size_t i =0; i< reverseSubchainNode.size(); i++)
                                reverseSubchainNode[i] = -reverseSubchainNode[i];
                        
                        currentNodeID = reverseSubchainNode[reverseSubchainNode.size()-2];
                        currentNode = dbg.getSSNode(currentNodeID);

                        if (currentNode.getNumRightArcs() >1)
                                markovTable[currentNodeID + numNodes].addNewChain(reverseSubchainNode, len);
                        
                }
        }
}
size_t MarkovChain::getChainLength(vector<NodeID> nodeChain){
        size_t len = 0;
        if (nodeChain.size() <= minMCOrder)
                return len;

        len = Kmer::getK() + 1; // for the first and the last item in the chain

        for (size_t k =1; k<nodeChain.size()-1; k++){
                SSNode node = dbg.getSSNode(nodeChain[k]);
                len = len + node.getMarginalLength();
        }

        return len;
}
vector< vector< vector< NodeID >> > MarkovChain::extractSubchains(const vector< NodeID> nodeID ){
        
        vector<vector< vector< NodeID> > > result;
        if (nodeID.size() < minMCOrder+1 )
                return result;
        //extract the subChain if the seed has at least two nodes.
        size_t maxSubchainSize = min (maxMCOrder+1, nodeID.size());

        for (size_t subchainSize = minMCOrder+1; subchainSize <= maxSubchainSize; subchainSize++)
        {
                vector<vector< NodeID >> subChainLists;
                for (size_t startIndex =0; startIndex <= nodeID.size() - subchainSize; startIndex++ )
                {
                        vector<NodeID> subChainNode;

                        for (size_t index = startIndex; index< startIndex + subchainSize; index++)
                                subChainNode.push_back(nodeID[index]);
                        subChainLists.push_back(subChainNode);
                }
                result.push_back(subChainLists);
        }
        return result;

}
void MarkovChain::findSeedKmer(const std::string& read,
                                     vector<KmerSeed>& mergedSeeds){

        vector<NodePosPair> nppv(read.length() + 1 - Kmer::getK());

        // find the node position pairs using the kmer lookup table
        findNPPFast(read, nppv);

        // transform consistent npps to seeds
        vector<KmerSeed> seeds;
        extractSeeds(nppv, seeds);

        // sort seeds according to nodeID
        sort(seeds.begin(), seeds.end());

        // merge seeds
        KmerSeed::mergeSeeds(seeds, mergedSeeds);


}
void KmerSeed::mergeSeeds(const vector<KmerSeed>& seeds,
                      vector<KmerSeed>& mergedSeeds){

        if (seeds.empty())
                return;

        mergedSeeds.reserve(seeds.size());
        mergedSeeds.push_back(seeds.front());

        for (size_t i = 1; i < seeds.size(); i++) {
                const KmerSeed& left = seeds[i-1];
                const KmerSeed& right = seeds[i];
                
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


void MarkovChain::findNPPFast(const string& read, vector<NodePosPair>& nppv)
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
void MarkovChain::extractSeeds(const vector<NodePosPair>& nppv,
                                     vector<KmerSeed>& seeds){

        size_t prev = nppv.size();
        for (size_t i = 0; i < nppv.size(); i++) {
                if (!nppv[i].isValid())
                        continue;

                // is it the first time we encounter a valid npp?
                if (prev == nppv.size()) {
                        prev = i;
                        seeds.push_back(KmerSeed(nppv[i].getNodeID(), nppv[i].getPosition(), i, i + 1));
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
                        seeds.push_back(KmerSeed(nppv[i].getNodeID(), nppv[i].getPosition(), i, i + 1));
                } else {
                        seeds.back().readEnd = i + 1;
                }
        }
}

/**
 * ChainRow class routines description
 **/


ChainRow::ChainRow(vector<NodeID>nodeChain_, size_t len, size_t freq): totalSampleSize(0), order(nodeChain_.size()-1), multiplicity(0), singleOutgoint(false){

        if (nodeChain_.size()< minMCOrder + 1 )
                return;
        NodeID nextNode = nodeChain_.back();
        nodeChain_.pop_back();
        increaseFreqOfNode(nextNode, len, freq);
        nodeChain = nodeChain_;
}


ChainRow::ChainRow(vector<NodeID> nodeChain_,size_t length_,
                   size_t sampleSize ,size_t multiplicity, vector< Decision> decisionVec_): order(nodeChain_.size()), multiplicity(0), singleOutgoint(false)
{
        if (nodeChain_.size()< minMCOrder )
                return;
        decisionVector = decisionVec_;
        nodeChain = nodeChain_;
        setTotalSampleSize(sampleSize);
        setLength(length_);
        setMultiplicity(multiplicity);
}


void ChainRow::updateAVGlength(float len, size_t freq){

        if (totalSampleSize == 0)
                length_avg = len;
        else
                length_avg = float ((double)(length_avg * totalSampleSize + len * freq )/ (double)(totalSampleSize + freq) ) ;
}



bool ChainRow::makeDecision( float expectedCov){
        std::sort(decisionVector.begin(),decisionVector.end(), sortDecisions());
        //there is not enough sample to decide
        if ( expectedCov < minimumExpectedCov)
                return false;
        
        set<NodeID> invalidNodeSet;
        size_t numOfValidOutgoing = decisionVector.size (); 
        for (size_t i = 0 ; i < decisionVector.size (); i++){
                double probMul_1 = Util::poissonPDF((double) decisionVector[i].getFrequency(), expectedCov);
                double probMul_0 = Util::poissonPDF((double) decisionVector[i].getFrequency(), errorMean);

                long double oddratio = probMul_0 /probMul_1 ;
                
                if (log10(oddratio) > minOddratioLog && decisionVector[i].getFrequency() < expectedCov){
                        /*#ifdef DEBUG
                        printRow();
                        cout << "probMul_1    : " << probMul_1 <<endl;
                        cout << "probMul_0    : "<<  probMul_0 <<endl;
                        cout << "oddratio     : "<<  oddratio  <<endl;
                        cout << "frequency    : " << decisionVector[i].getFrequency() <<endl;
                        cout << "expectedCov  : " <<expectedCov <<endl;
                        #endif*/
                        decisionVector[i].setDecision(false);
                        decisionVector[i].setOddRatioOfAcceptance(oddratio);
                        invalidNodeSet.insert(decisionVector[i].getNodeID());
                        numOfValidOutgoing --;
                }

        }
        if (numOfValidOutgoing == 1)
                singleOutgoint = true;
        //the repeat is completely resolved, so announce it to the higher order of MM  tables 
        if (invalidNodeSet.size() >0 ){
                std::set<NodeID>::iterator it;
                for (size_t i = 0 ; i < decisionVector.size (); i++){
                        it = invalidNodeSet.find( decisionVector[i].getNodeID() );
                        if (it == invalidNodeSet.end()){
                                decisionVector[i].setDecision(true);                 
                        }else{
                                decisionVector[i].setDecision(false);                 
                        }
                }
                return true;
        }
        
        return false;
        
}
void ChainRow::increaseFreqOfNode(NodeID nodeID,size_t len, size_t freq  ){

        bool found = false;
        updateAVGlength(len, freq);
        for (size_t i=0; i < decisionVector.size(); i++){
                if (decisionVector[i].getNodeID() == nodeID){
                        decisionVector[i].setFrequency ( decisionVector[i].getFrequency() + freq);
                        found = true;
                }
        }
        if (!found){
                Decision newDecision(nodeID, freq);
                decisionVector.push_back(newDecision);
        }
        totalSampleSize = totalSampleSize + freq;
}

bool ChainRow::equal (vector<NodeID> nodeList){

        if (nodeChain.size()!= nodeList.size())
                return false;

        for (size_t i =0; i < nodeList.size(); i++){
                if (nodeList[i] != nodeChain[i]){
                        return false;
                }
        }
        return true;
}
bool ChainRow::equalEnd(vector<NodeID> nodeList){
        if (nodeList.size() > nodeChain.size())
                return false;
        for (size_t i = 0; i < nodeList.size(); i++){
                if (nodeList[i] != nodeChain[nodeChain.size()-nodeList.size() + i]){
                        return false;
                }
        }
        return true;
}

void ChainRow::printRow (){

        std::sort(decisionVector.begin(),decisionVector.end(), sortDecisions());
        for ( size_t i=0; i<decisionVector.size();i ++)
        {
                Decision nextNodes = decisionVector[i];
                for (size_t j=0; j<nodeChain.size(); j++){
                        cout <<nodeChain[j] << " , " ;
                }
                cout <<"==> " <<nextNodes.getNodeID() << " , freq ( " <<nextNodes.getFrequency() <<std::setprecision(2)
                << " ) , prob( " << float ((float)nextNodes.getFrequency() /(float) totalSampleSize) << " )";
                if (nextNodes.getDecision())
                        cout << " , this decision has been chosen. The odd ratio of correctness ( "<<nextNodes.getOddRatioOfAcceptance() <<" )"<<endl;
                else
                        cout <<endl;
        }

        cout << "Total sample size is : " <<totalSampleSize << "  , the min length of the seeds which cover this chain is: " <<int( length_avg)<<" The guessed multiplicity is: "<<getMultiplicity() <<endl ;
     
}

/**
 * TransitionTable class routines description
 **/


bool TransitionTable:: smartSearchRow(vector<NodeID> chainList, ChainRow& row){

        if (!modified)
                return false;
        if (rows.size() ==0)
                return false;
        size_t i = 0;
        int  mid = 0, low =0 , high = rows.size()-1;
        int  prv_low =0 , prv_high = rows.size()-1;
        while (i < chainList.size() ){
                NodeID value = chainList [i];
                while (low <= high){
                        mid = (low + high) / 2;
                        if (rows[mid].nodeChain[i] == value){
                                low = mid;
                                high = mid;
                                while (low-1 >= prv_low && rows[low-1].nodeChain[i] == value)
                                        low --;
                                while (high+1 <= prv_high && rows[high+1].nodeChain[i] == value)
                                        high ++;
                                break;
                        }else if (mid > prv_low && rows[mid].nodeChain[i] > value)
                                high = mid - 1;
                         else if (mid < prv_high && rows[mid].nodeChain[i] < value)
                                low = mid + 1 ;
                         else
                                 return false;
                }
                prv_high = high;
                prv_low = low;
                i++;

        }
        if ( rows[mid].equal(chainList)){
                row = rows[mid];
                return true;
        }
        return false;
}


void TransitionTable::addRowSlow(vector<NodeID> chainList, size_t len){

        if (!modified)
                modified = true;
        order = chainList.size() - 1;
        bool found = false;
        NodeID nextNode = chainList.back();
        chainList.pop_back();
        for (size_t i = 0; i <rows.size() ;i ++){
                if ( rows[i].equal(chainList)){
                        rows[i].increaseFreqOfNode(nextNode, len, 1);
                        found = true;
                        break ;
                }
        }
        if (!found){
                chainList.push_back(nextNode);
                ChainRow newRow(chainList, len, 1);
                rows.push_back(newRow);
        }

}

void TransitionTable::addNewRow(ChainRow row){

        if (!modified)
                modified = true;
        order = row.nodeChain.size();
        rows.push_back(row);
}


void TransitionTable::printRows(){
        if (rows.empty())
                return;

        sort();

        cout << "The order of this table is : " <<order <<endl ;
        for (ChainRow row : rows){
                row.printRow();
                cout << "               *********                       " <<endl;
        }

}

void TransitionTable::sort(){
         std::sort(rows.begin(),rows.end(), sortRows());
}

void TransitionTable::mergeTables(TransitionTable &subTable){

        if (subTable.modified && !modified)
                modified = true;

        for (ChainRow row : subTable.rows){

                for ( size_t i=0; i<row.decisionVector.size();i ++)
                {
                        Decision decision = row.decisionVector[i];
                        bool found = false;
                        vector<NodeID> chainList = row.nodeChain;

                        this->order = chainList.size();
                        NodeID nextNode = decision.getNodeID();
                        size_t freq = decision.getFrequency();
                        size_t avgLen = row.getLength();

                        for (size_t i = 0; i <rows.size() ;i ++){
                                if ( rows[i].equal(chainList)){
                                        rows[i].increaseFreqOfNode(nextNode,avgLen ,freq);
                                        found = true;
                                        break ;
                                }
                        }
                        if (!found){
                                chainList.push_back(nextNode);
                                ChainRow newRow(chainList,avgLen, freq);
                                rows.push_back(newRow);
                        }
                }

        }

}
void TransitionTable::updateDecisions (ChainRow row_){
        set <NodeID> invalidNodeSet ;
        for (size_t i = 0 ; i < row_.decisionVector.size (); i++){
                if (!row_.decisionVector[i].getDecision())
                        invalidNodeSet.insert(row_.decisionVector[i].getNodeID());
        }
        
        for (auto  row :rows){
                if( row.equalEnd(row_.nodeChain)){
                        for (size_t i = 0 ; i < row.decisionVector.size(); i++){
                                auto search = invalidNodeSet.find(row.decisionVector[i].getNodeID());
                                if (search != invalidNodeSet.end()){
                                        row.decisionVector[i].setDecision(false);  
                                }else{
                                        row.decisionVector[i].setDecision(true);  
                                }
                        }                 
                }
                
        }
                
}
void TransitionTable::forceRemove(vector<NodeID> nodeChain){

        for (vector<ChainRow>::iterator it= rows.begin(); it != rows.end(); )
        {

                if(it->equalEnd(nodeChain))
                        it = rows.erase(it);
                else
                        ++it;
        }
        if (rows.size() ==0 )
                modified = false;

}
/**
 * TableContainer class routines description
 **/
void TableContainer:: addNewChain(vector<NodeID> chainList, size_t len){

        if (modified == false)
                modified = true;
        size_t order = chainList.size() -1;
        if (order< minMCOrder || order > maxMCOrder)
                return ;
        orderTable [order -minMCOrder].addRowSlow(chainList, len);
}
void TableContainer::addNewRow(ChainRow &row, size_t order){
        if (modified == false)
                modified = true;
        orderTable [order - minMCOrder].addNewRow(row);

}
void TableContainer::printTable(){

        if (minMCOrder < 1 )
                return ;
        for (size_t order = minMCOrder; order< maxMCOrder ; order++){
                if (orderTable[order-minMCOrder].modified){
                        orderTable[order-minMCOrder].printRows();
                        cout << "       ******************" <<endl;
                }
        }

}
void TableContainer::mergeTables(TableContainer &subTables){

        if (minMCOrder < 1 )
                return ;
        if (subTables.modified && ! modified)
                modified = true;
        for (size_t order = minMCOrder; order< maxMCOrder ; order++){
                orderTable[order-minMCOrder].mergeTables(subTables.orderTable[order-minMCOrder]);
        }

}
void TableContainer::propagateDecision( ChainRow row){
        size_t currentOrder = row.getOrder();
        if (row.getSingleValidOut()){
                // if there is only one valid outgoing arc, so remove the higher ordres in the table
                for (size_t order = currentOrder + 1; order< maxMCOrder ; order++){
                        orderTable[order-minMCOrder].forceRemove(row.nodeChain);
                }
        }
        else{
                // if there is more than one valid outgoing arc, only uncheck the false ones
                for (size_t order = currentOrder + 1; order< maxMCOrder ; order++){
                        orderTable[order-minMCOrder].updateDecisions(row);
                }
        }
        
        
}
void TableContainer::sort()
{
        for (size_t order = minMCOrder ; order< maxMCOrder ; order++)
                orderTable[order-minMCOrder].sort();
}





void RepeatResolver::updateExpectedCovByLen(bool errorFreeAssumtion){

        expectedCoverageByLen [maxReadSize - Kmer::getK()- minMCOrder];
        
        for (size_t i = Kmer::getK() + minMCOrder; i<= maxReadSize; i++ ){
                double errorEffect = errorFreeAssumtion ? pow(1-illuminaErrorRate, i):1;
                
                double expectedCov = (float)((float)estiCoverage/(float)avg_readLength)* (float)(avg_readLength - i + 1)* errorEffect;
                expectedCoverageByLen.push_back(expectedCov);
                //cout <<"len: " <<i << " : expected Cov: " <<std::setprecision(1)<<std::fixed<< expectedCov <<endl;

        }
}
float RepeatResolver::getExpectedCovByLen(size_t len){
        if (len < Kmer::getK()+ minMCOrder || len > maxReadSize)
                return 0;
        return (expectedCoverageByLen[len - Kmer::getK()-minMCOrder])> 0 ? (expectedCoverageByLen[len - Kmer::getK()-minMCOrder]) :0;
}

void RepeatResolver ::fitModelIntoTable(TableContainer &tables){
        if(! tables.modified)
                return ;
        makeDecision(tables);

}

void RepeatResolver::makeDecision(TableContainer &tables){
        for (size_t order = minMCOrder; order< maxMCOrder ; order++){
                if (tables.orderTable[order-minMCOrder].modified){
                        for (size_t i = 0; i <tables.orderTable[order-minMCOrder].rows.size(); i++){
                                ChainRow& row = tables.orderTable[order-minMCOrder].rows[i];
                                size_t seedLength = row.getLength();
                                float expectedCov = getExpectedCovByLen(seedLength);
                                if (row.makeDecision(expectedCov)) {
                                        // if a decision is made, then inform the higher order
                                        tables.propagateDecision(row);
                                }
                                
                                
                        }
                }
        }
}



void RepeatResolver::fitModelIntoTableThread(size_t threadID, ParGraph& wlb){


        while (true) {
                size_t firstNode, chunkSize;
                wlb.getNodeChunk(firstNode, chunkSize);

                if (chunkSize == 0)
                        break;

                //cout << "Work from " << threadID << " " << firstNode << " to " << firstNode + chunkSize << endl;
                for (size_t id = firstNode; id < firstNode + chunkSize; id++) {
                  
                        if (id ==0)
                                continue;
                        fitModelIntoTable(markovTable[-id + numNodes]);
                        fitModelIntoTable(markovTable[ id + numNodes]);

                }
        }

}

void RepeatResolver::fitModelIntoTable(bool errorFreeAssumtion){
        
        
        updateExpectedCovByLen(errorFreeAssumtion);
        const unsigned int& numThreads =  settings.getNumThreads();
        ParGraph wlb(numNodes, settings.getThreadWorkSize());
        
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++){
                workerThreads[i] = thread(&RepeatResolver::fitModelIntoTableThread, this,i, ref(wlb));
        }
        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
        cout << "\tProcessing table (100%) " << endl;
}


void ParGraphMCH::getNodeChunk(size_t& chunkOffset, size_t& chunkSize)
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
        cout << std::fixed << std::setprecision(1) << "\tProcessing table (" << perc << "%)\r";
        cout.flush();

        lock.unlock();
}
