
#include "library.h"
#include <thread>
#include "alignment.h"

//using namespace std;
class SamRecordAlternative{
public :
        string RNAME;              //Reference sequence NAME               \*|[!-()+-<>-~][!-~]* 
        string CIGAR;              //CIGAR string                          \*|([0-9]+[MIDNSHPX=])+
      
        size_t NM;
        size_t POS;                //1-based leftmost mapping POSition     [0,231-1]
        bool read_rc;

        int similarityScore;
        SamRecordAlternative(string RNAME_,string  CIGAR_, size_t POS_, size_t NM_, bool read_rc_ , size_t similarityScore_):RNAME(RNAME_),CIGAR(CIGAR_), NM(NM_), POS(POS_), read_rc(read_rc_), similarityScore(similarityScore_){};       
};
class SamRecord{
public:
     string QNAME;              //Query template NAME                   [!-?A-~]{1,254} 
     size_t FLAG;               //bitwise FLAG                          [0,216-1]
     string RNAME;              //Reference sequence NAME               \*|[!-()+-<>-~][!-~]*
     size_t POS;                //1-based leftmost mapping POSition     [0,231-1]
     size_t MAPQ;               //MAPping Quality                       [0,2^8-1]
     string CIGAR;              //CIGAR string                          \*|([0-9]+[MIDNSHPX=])+
     string RNEXT;              //Ref. name of the mate/next read       \*|=|[!-()+-<>-~][!-~]*
     size_t PNEXT;              //Position of the mate/next read        [0,2^31-1] 
     size_t TLEN;               //observed Template LENgth              
     string SEQ;                //segment SEQuence                      
     string QUAL;               //ASCII of Phred-scaled base QUALity+33  [!-~]+
     size_t NM;                 //edit distance 
     string XA;                 // alternatives hits
     int similarityScore;        // number of matches - NM
     bool alternativeIsValid;
     bool checkAlternative;
     vector <SamRecordAlternative> alternatives;
     SamRecord (vector<char> supportedCigar, bool cheackAlternative_):QNAME(""),FLAG(0), RNAME(""), POS(0),MAPQ (0),CIGAR("*"), RNEXT(""), TLEN(0),SEQ(""), QUAL(""),NM(0),XA(""),similarityScore(0) ,alternativeIsValid(false), checkAlternative(cheackAlternative_),read_paired(false), read_mapped(false),read_unmapped(false), mate_unmapped(false), read_rc(false), mate_rc(false),first_inPair(false), second_inPair(false), secondaryAlignment(false), supplementaryAlignment(false) , mapped(true){
             cigarChars.push_back('M');
             cigarChars.push_back('I');
             cigarChars.push_back('D');
             cigarChars.push_back('N');
             cigarChars.push_back('S');
             cigarChars.push_back('H');
             cigarChars.push_back('P');
             cigarChars.push_back('X');
             for (size_t i = 0; i<supportedCigar.size(); i++){
                std::vector<char>::iterator it;
                it = find (cigarChars.begin(), cigarChars.end(),toupper(supportedCigar[i]));
                if (it != cigarChars.end()){
                        supportedCigarChars.push_back(toupper(supportedCigar[i]));
                }
                     
        }
     };
     
     bool parseRecord (const string & line);
     void parseXAtag ();
     void interpretFlag();
     string getBinaryFlagRev(int flag);
    
     vector<string> splitString(const string str,const char delimiter);
     vector<pair< string, char>> splitString(const string str,vector <char> delimiters);
     vector<pair<size_t,char>> cigar_v;
     vector <char> cigarChars;
     vector<char> supportedCigarChars;
    
     string FlagBitwiseRev; 
     bool read_paired;            // 1 0x1 template having multiple segments in sequencing
     bool read_mapped;            // 2 0x2 each segment properly aligned according to the aligner
     bool read_unmapped;          // 4 0x4 segment unmapped
     bool mate_unmapped;          // 8 0x8 next segment in the template unmapped
     bool read_rc;                // 16 0x10 SEQ being reverse complemented
     bool mate_rc;                //32 0x20 SEQ of the next segment in the template being reverse complemented
     bool first_inPair;           //64 0x40 the first segment in the template
     bool second_inPair;          //128 0x80 the last segment in the template
     bool secondaryAlignment;     //256 0x100 secondary alignment
     bool supplementaryAlignment;  // 2048 supplementary alignment
     size_t cigarLen;
     bool mapped;

};


class SamToAlignment{
private:
         vector< pair<string, string> >  &refSeq_v;      // reference sequences
         vector<char> &supportedCigarChars;
         bool checkAlternative;
         
public:
        size_t unmappedReads ;
        SamToAlignment( vector< pair<string, string> >  &refSeq_in, vector<char>  &supportedCigarChars_in, bool checkAlternative_);
        void convert (  vector<ReadRecord> &myReadBuf, vector<SamRecord>mysamBuf);
        bool convert(ReadRecord &r , SamRecord s);
        void readSamFile(size_t& pos, std::vector< SamRecord >& samChunk, size_t chunkSize, string samFileName);
        string getGenomeStr(string refId, size_t pos, size_t length);
        vector<string> splitString(const string str, const char delimiter);
        string reverseComplement(string input);
};
class SamToAlignmentHandler{
private:
        vector< pair<string, string> >  refSeq_v;      // reference sequences
        size_t unmappedReads ;
        std::mutex mergeMutex;               // mutex
        vector<char> supportedCigarChars;
public:
        
        SamToAlignmentHandler(string genomeFileName, string cigarChars);
        void convert (string samFileName, string inputReadFileName, bool checkAlternative);
        void loadGenomeFile(string genomeFileName);
        void workerThread(size_t myID, LibraryContainer &libContSam, string samFileName, bool checkAlternative);
        void updateStatistic(SamToAlignment &sam2alignment);
        void extractMapped(string inputFileName, string outputFileName);
};


