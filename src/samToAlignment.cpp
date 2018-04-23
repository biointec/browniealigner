
#include "samToAlignment.h"
#include "readfile/fastafile.h"
#include <thread>
#include <iostream>
#include <fstream>
#include <set>

 
 SamToAlignmentHandler::SamToAlignmentHandler ( string genomeFileName, string cigarChars): unmappedReads(0){
        loadGenomeFile(genomeFileName);
        for (size_t i = 0 ; i < cigarChars.length(); i++)
                supportedCigarChars.push_back(cigarChars[i]);
        
 }
 
 void SamToAlignmentHandler::loadGenomeFile(string genomeFileName){
         std::ifstream input(genomeFileName);
         
         if (!input.good()) {
                 std::cerr << "Error opening: " << genomeFileName << std::endl;
                 exit(0);
         }
         std::string line, id ="", preid = "", DNA_sequence ="" ;
         while (std::getline(input, line).good()) {
                 if (line[0] == '>') {
                         preid = id;
                         id = line.substr(1);
                         if (DNA_sequence != "")
                                refSeq_v.push_back( make_pair(preid, DNA_sequence));
                         DNA_sequence.clear();
                 }
                 else if (line[0] != '>'){
                         DNA_sequence += line;
                         
                 }
         }
         refSeq_v.push_back( make_pair(id, DNA_sequence));
         cout << refSeq_v.size() <<  " number of chromosome were added !!" <<endl;
 }

void SamToAlignmentHandler::workerThread(size_t myID, LibraryContainer &libContSam, string samFileName, bool checkAlternative){

        SamToAlignment sam2alignment(refSeq_v, supportedCigarChars,checkAlternative);
        // local storage of reads
        vector<ReadRecord> myReadBuf;
        vector<SamRecord>mysamBuf;
        size_t samFilePos = 0;
        bool result = true;
        while (result) {
                size_t blockID, recordID;
                result = libContSam.getRecordChunk(myReadBuf, blockID, recordID);
                //ReadRecord r = myReadBuf[myReadBuf.size()];
                sam2alignment.readSamFile(samFilePos, mysamBuf,myReadBuf.size(), samFileName);
                sam2alignment.convert(myReadBuf, mysamBuf);
                mysamBuf.clear();
                if (result)
                        libContSam.commitRecordChunk(myReadBuf, blockID, recordID);
                else
                        break;
        }
        updateStatistic(sam2alignment);
}
void SamToAlignmentHandler::updateStatistic(SamToAlignment &sam2alignment){
        lock_guard<mutex> lock(mergeMutex);
        unmappedReads  = unmappedReads+ sam2alignment.unmappedReads;
}
void SamToAlignmentHandler::extractMapped(string inputFileName, string outputFileName){
        std::ifstream input(inputFileName);
        if (!input.good()) {
                std::cerr << "Error opening: " << inputFileName << std::endl;
                exit(0);
        }
         std::string line, id ="", preid = "", DNA_sequence ="" ;
         while (std::getline(input, line).good()) {
                 if (line[0] == '>') {
                         preid = id;
                         id = line.substr(1);
                         //cout << id <<endl;
                         if (DNA_sequence != "")
                                refSeq_v.push_back( make_pair(preid, DNA_sequence));
                         DNA_sequence.clear();
                 }
                 else if (line[0] != '>'){
                         DNA_sequence += line;
                         
                 }
         }
        
    
}
void SamToAlignmentHandler::convert(string samFileName,string inputReadFileName,bool checkAlternative)
{
        size_t numOfThreads = 1;
        LibraryContainer libContRead;
        libContRead.insert(ReadLibrary(inputReadFileName,"","."));
         size_t threadSize = 100000;
        libContRead.startIOThreads(threadSize, 10 * threadSize * numOfThreads, true);    
         // start worker threads
        vector<thread> workerThreads(numOfThreads);
        for (size_t i = 0; i < workerThreads.size(); i++){
                workerThreads[i] = thread(&SamToAlignmentHandler::workerThread,this, i, ref(libContRead), samFileName, checkAlternative);
        }
        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
        libContRead.joinIOThreads();
        cout << "number unmapped reads :" << unmappedReads <<endl;
        
}
SamToAlignment::SamToAlignment ( vector< pair<string, string> >  &refSeq_in, vector<char>&  supportedCigarChars_in, bool checkAlternative_):refSeq_v(refSeq_in),supportedCigarChars(supportedCigarChars_in) ,  checkAlternative(checkAlternative_), unmappedReads(0){
        
}
string SamToAlignment::getGenomeStr(string refId, size_t pos ,size_t length){
        string str = "";
        for (size_t i = 0; i < refSeq_v.size(); i++){
                     if (splitString( refSeq_v[i].first, ' ')[0] == refId){
                                int minLen = (int)refSeq_v[i].second.length() <  (int)length + (int)pos ? (int)refSeq_v[i].second.length()- (int)pos :length  ;
                                if (minLen > 0)
                                        str = refSeq_v[i].second.substr ( pos , minLen);
                                else 
                                        return "";
                     }
        }
        return str; 
}
void SamToAlignment::convert(vector<ReadRecord> &readChunk, vector<SamRecord> mySamBuf)
{
        size_t i = 0;
        for (ReadRecord& it : readChunk){
                
                if (!convert(it,mySamBuf[i])){
                        mySamBuf[i].mapped = false;
                        unmappedReads ++;
                        it.preRead.insert(it.preRead.length()-1,"\tA(0)A");        
                }else{
                        it.preRead.insert(it.preRead.length()-1,"\tA(1)A");        
                }
                

                i++;
        }
}
string SamToAlignment::reverseComplement(string str)
{
        std::reverse(str.begin(), str.end());
        string complement = "";
        for (size_t i = 0; i < str.size(); i++) {
                if (str[i] == 'A' or  str[i] == 'a'){
                      complement = complement + 'T';  
                }
                else if (str[i] == 'T' or  str[i] == 't'){
                      complement = complement + 'A';  
                }
                else if (str[i] == 'C' or  str[i] == 'c'){
                      complement = complement + 'G';  
                }
                else if (str[i] == 'G' or  str[i] == 'g'){
                      complement = complement + 'C';  
                }
                else if (str[i] == 'N' or  str[i] == 'n'){
                      complement = complement + 'N';  
                }
                else {
                        complement = complement +str[i];
                }
                
        }
        return complement;
}

bool SamToAlignment::convert(ReadRecord &r , SamRecord s){
        //cout << splitString( splitString( splitString(r.preRead,'\n' )[0], '@')[1], '/')[0] <<endl;
 
        if (!s.mapped )
                return false;

        if ( s.QNAME !=splitString( splitString( splitString( splitString(r.preRead,'\n' )[0], '@')[1], '/')[0],' ')[0] ){
                cout << splitString( splitString( splitString( splitString(r.preRead,'\n' )[0], '@')[1], '/')[0],' ')[0] <<endl;;
                cout << "Different id  "<< s.QNAME <<endl;
                return false;
        }
        string read = r.read, corrected = "";
        if (s.read_rc)
                read = reverseComplement(r.read);
        string ref = getGenomeStr(s.RNAME,s.POS-1, r.read.length()*2);
        if (ref == "")
                return false;

        size_t refPos = 0;
        size_t readPos = 0;
        if (s.cigar_v[0].second == 'S'){
                readPos = s.cigar_v[0].first; 
                for (size_t i = 0; i< s.cigar_v[0].first; i++){
                        corrected = corrected + "N";
                }
        }

        for (size_t i  =0; i< s.cigar_v.size(); i++){
                if (s.cigar_v[i].second == 'M'){
                        if (ref.length()<refPos + s.cigar_v[i].first)
                                return false;
                        corrected = corrected + ref.substr(refPos, s.cigar_v[i].first) ;
                        readPos = readPos + s.cigar_v[i].first ;
                        refPos = refPos + s.cigar_v[i].first ;
                }
                if (s.cigar_v[i].second == 'S' && i >0){
                        readPos = readPos + s.cigar_v[i].first ;
                        refPos = refPos + s.cigar_v[i].first ;
                        for (size_t p = 0; p< s.cigar_v[i].first; p++){
                                corrected = corrected + "N";
                        }
                }
                if (s.cigar_v[i].second == 'I')
                        readPos = readPos + s.cigar_v[i].first ;
                
                if (s.cigar_v[i].second == 'D'){
                        if (refPos > ref.length()){
                                return false;
                        }
                        corrected = corrected + ref.substr(refPos, s.cigar_v[i].first) ;
                        refPos = refPos + s.cigar_v[i].first ;                                
                }
        }  
        if (s.read_rc)
                corrected = reverseComplement(corrected);
   
        // there is deletion in the cigar, so 
        if (corrected.length() > r.read.length()){
                if (!s.read_rc)
                        corrected = corrected.substr(0, r.read.length());
                else
                        corrected = corrected.substr(corrected.length() - r.read.length(), r.read.length());
        }
        // there is insertion in the cigar
        if (corrected.length() < r.read.length()){
                if (ref.length() < refPos + r.read.length() - corrected.length())
                        return false;
                if (s.read_rc)
                        corrected = corrected + ref.substr(refPos,r.read.length() - corrected.length()) ;
                else{
                        size_t dif = r.read.length() - corrected.length();
                        int newPos =(int)s.POS-1 -(int)dif;
                        newPos = max(newPos, 0);
                        string str = getGenomeStr(s.RNAME, newPos, dif);
                        corrected = str + corrected;
                }


        }

        r.read = corrected;
        
        return true;
}



void SamToAlignment::readSamFile ( size_t& pos, std::vector< SamRecord >& samChunk, size_t chunkSize, string samFileName )
{
        size_t numOfItems = 0;
        std::ifstream ifs(samFileName.c_str());
        ifs.seekg( pos );
        std::set<string> ids;
        while(ifs.is_open())
        {
                char c = ifs.peek();
                string line = "";

                while(c=='@'){
                        std::getline(ifs, line);
                        c = ifs.peek();
                }
                std::getline(ifs, line);
   
                if (line ==""){
                        pos  = ifs.tellg();
                        ifs.close();
                        return;
                }
                SamRecord record(supportedCigarChars,checkAlternative);

                if (!record.parseRecord(line))
                        record.mapped = false;
                
                if (!record.secondaryAlignment && !record.supplementaryAlignment){
                        samChunk.push_back(record);
                        numOfItems ++;
                }
                if (numOfItems >= chunkSize ){
                        pos  = ifs.tellg();
                        ifs.close();
                        return;
                }
        }
        ifs.close();

} 
vector<string> SamRecord::splitString(const string str,const char delimiter)
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
vector<string> SamToAlignment::splitString(const string str,const char delimiter)
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
vector<pair< string, char>> SamRecord::splitString(const string str,vector <char> delimiters)
{
        vector<pair< string, char>> parsed;
        size_t pos = 0;
        string currentStr = "";
        while (pos < str.length()) {
                char c = str [pos];
                std::vector<char>::iterator it;
                it = find (delimiters.begin(), delimiters.end(),c);
                if (it != delimiters.end()) {
                        parsed.push_back( make_pair(currentStr, c));
                        currentStr = "";
                        
                }else{
                        currentStr = currentStr + c;
                }
              pos++;
        }
        return parsed;
}

bool SamRecord::parseRecord(const string & line)
{
        vector<string> items  = splitString(line, '\t');   
        if (items.size() < 11){
                //cout << "wrong format" ;
                return false;
        }
        QNAME = items[0];             
        FLAG  = stoi(items[1]);               
        RNAME = items[2];              
        POS   = stoi(items[3]);                
        MAPQ  = stoi(items[4]);               
        CIGAR = items[5];              
        RNEXT = items[6];              
        PNEXT = stoi(items[7]);              
        TLEN  = stoi(items[8]);               
        SEQ   = items[9];                
        QUAL  = items[10]; 
        if (items.size() >11){
                string temp = items[11];
                if (temp[0] =='N' && temp[1]=='M'){
                        vector<string> elements = splitString(temp, ':');
                        if (elements.size ()<3)
                                return false;
                        NM = atoi(elements[2].c_str());
                        
                }
                if (temp[0] =='X' && temp[1]=='A'){
                        XA = temp;
                }
        }
        if (items.size() > 12){
                string temp = items[12];
                if (temp[0] =='N' && temp[1]=='M'){
                        vector<string> elements = splitString(temp, ':');
                        if (elements.size () < 3)
                                return false;
                        NM = atoi(elements[2].c_str());
                }
                if (temp[0] =='X' && temp[1]=='A'){
                        XA = temp;
                }
        }
        interpretFlag(); 
        if (XA != "" && checkAlternative){
                parseXAtag();
        }
        if (CIGAR == "*" )
                return false;
        
        vector<pair<string , char>> cigarItems = splitString(CIGAR, cigarChars);
        for (size_t i = 0; i< cigarItems.size(); i++){
                if (cigarItems[i].second == 'M'){
                         similarityScore = similarityScore +  atoi (cigarItems[i].first.c_str());
                }
        }
        similarityScore = similarityScore - NM;
        for (auto al : alternatives){
                if (al.similarityScore > similarityScore){
                         RNAME = al.RNAME;
                         POS = al.POS;
                         similarityScore = al.similarityScore;
                         read_rc = al.read_rc;
                         CIGAR = al.CIGAR;
                         NM = al.NM;
                         alternativeIsValid = true;
                    
                }
                
        }
        cigarItems = splitString(CIGAR, cigarChars);
        for (size_t i = 0; i< cigarItems.size(); i++){
                 cigar_v.push_back(make_pair(atoi( cigarItems[i].first.c_str()) , cigarItems[i].second) );
                 std::vector<char>::iterator it;
                 it = find (supportedCigarChars.begin(), supportedCigarChars.end(),cigarItems[i].second);
                 if (it ==supportedCigarChars.end()){
                         return false;
                 }
        }
        return true;
        
}
void SamRecord::parseXAtag()
{

        vector<string>  XAs= splitString(XA, ';');
        if (XAs.size() >5){
                int stop = 0;
                stop ++;
        }
        for (auto xa : XAs){
                vector <string> fields = splitString(xa, ',');
                if (fields.size() == 4){
                        string refName = fields[0];
                        if (refName.size() >4){
                                if( refName [0] == 'X' && refName [1]=='A' && refName [2]==':' && refName [3]=='Z'&&refName[4] ==':'){
                                        refName = refName.substr(5);
                                }
                        }
                        string pos_s = fields[1];
                        string cigar = fields[2];
                        string NM_s = fields[3];
                        size_t pos = 0;
                        size_t NM = 0, simScore = 0;
                    
                        bool reverseComplement = false;
                        if (cigar !="*" && pos_s !="" && refName !="" && NM_s !=""){
                                if (pos_s[0]=='-'){
                                        reverseComplement = true;    
                                }
                                else{
                                        reverseComplement = false;  
                                }
                                pos = atoi (pos_s.substr(1, pos_s.length()).c_str());
                                NM = atoi (NM_s.c_str());
                                bool supported = true;
                                vector<pair<string , char>> cigarItems = splitString(cigar, cigarChars);
                                for (size_t i = 0; i < cigarItems.size(); i++){
                                       
                                        std::vector<char>::iterator it;
                                        it = find (supportedCigarChars.begin(), supportedCigarChars.end(), cigarItems[i].second);
                                        if (it == supportedCigarChars.end()){
                                                supported = false;
                                                break;
                                        }
                                        if (cigarItems[i].second == 'M')
                                                simScore = simScore + atoi (cigarItems[i].first.c_str());
                                }
                                if (supported){
                                        simScore = simScore - NM;
                                        SamRecordAlternative alternative(refName, cigar, pos, NM, reverseComplement ,simScore);
                                        alternatives.push_back(alternative);
                                }
                        }
                }
        }
        
}

string SamRecord::getBinaryFlagRev(int flag){
    std::string r;
    while(flag != 0 ) {
            r = (flag % 2==0 ?"0":"1") + r; 
            flag/=2;
    }
    std::reverse(r.begin(), r.end());
    while(r.length() < 12)
            r = r + "0";
    return r;  
}

void SamRecord::interpretFlag(){
        
        FlagBitwiseRev = getBinaryFlagRev(FLAG);
        if (FlagBitwiseRev[0] == '1')
                read_paired = true;
        if (FlagBitwiseRev[1] == '1')
                read_mapped = true;
        if (FlagBitwiseRev[2] == '1')
                read_unmapped = true;
        if (FlagBitwiseRev[3] == '1')
                mate_unmapped = true;
        if (FlagBitwiseRev[4] == '1')
                read_rc = true;
        if (FlagBitwiseRev[5] == '1')
                mate_rc = true;
        if (FlagBitwiseRev[6] == '1')
                first_inPair = true;
        if (FlagBitwiseRev[7] == '1')
                second_inPair = true;
        if (FlagBitwiseRev[8] == '1')
                secondaryAlignment = true; 
        if (FlagBitwiseRev[11] == '1')
                supplementaryAlignment = true;
}
int main(int argc, char ** args)
{
        cout << "Welcome to the sam2alignment software\n" << endl;
        cout << "Today is " << Util::getDateTime() << endl;
        cout << "----------------------------------------------------\n" << endl;
        if (argc < 5) {
                cout << "ERROR: Wrong number of arguments!\n";
                cout << "Please specify the allowed cigar chars in the last argument."<<endl;
                cout << "Usage: ./sam2alignment samfile.sam genome.fasta inputRead.fastq MIDS <optional>" << endl << flush;
                cout << "optional:" <<endl;
                cout << "-A\t consider alternative alignment if it has a lower edit distance the primary one" <<endl;
                cout << "-d\t maximum length difference between the original and corrected read." <<endl;
                return 0;
        }
        string samFileName = args[1];
        string genomeFileName  = args[2];
        string inputReadFileName  = args[3];
        string cigarChars  = args[4];
        
        bool checkAlternative = false;
        size_t maxLenDif = 0;
        for (int i = 5 ; i< argc; i++){
                string optional = args[i];
                if (optional == "-A"){
                        checkAlternative = true;
                }
                if (optional == "-d"){
                        maxLenDif = stoi (args[i+1]);
                        cout << maxLenDif <<endl;
                        i++;
                }
        } 
        Util::startChrono();
        
        SamToAlignmentHandler s2a(genomeFileName, cigarChars);
        s2a.convert(samFileName, inputReadFileName , checkAlternative);
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << "Exiting... bye!" << endl << endl;
        return 0;
}
