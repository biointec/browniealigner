#include "library.h"
#include <thread>
#include <string>
#include "alignment.h"
using namespace std;

class ComparisonMetrics
{
private:
        size_t numReads;                // number of reads handled
        size_t numOfcomparedReads;
        std::mutex metricMutex;         // mutex for merging metrics
        size_t TP;
        size_t FP;
        size_t FN;
        size_t TN;
        size_t TPFullRec;
        size_t FPFullRec;
        size_t FNFullRec;
        size_t TNFullRec;
        size_t maximumNumOfError;
        size_t numNewIntroducedN;
        size_t numInitialN;
public:
        double gain;
        double gainFR;

        /**
         * Default constructor
         */
        ComparisonMetrics() : numReads(0),numOfcomparedReads(0),TP(0),FP(0),FN(0),TN(0)
        ,TPFullRec(0),FPFullRec(0),FNFullRec(0), TNFullRec(0),maximumNumOfError(0), numNewIntroducedN(0), numInitialN(0), gain(0),gainFR(0){}

        /**
         * compare two reads base per base and
         * @aligned extra information some tools might have , if the read is aligned by tool
         *
         */
        void compareRead(string  correctedRead , string  erroneousRead, string   perfectRead, bool aligned);

        /**
         * Add other metrics (thread-safe)
         * @param metrics Metrics to add
         */
        void addMetrics(ComparisonMetrics& tMetrics);

        /**
         * Get the number of reads
         * @return The number of reads
         */
        size_t getNumReads() const {
                return numReads;
        }

        size_t getNumOfCompReads(){
               return numOfcomparedReads;
        }
        /**
         * prints out all calculated metrics after comparison
         *
         */
        void printStatistics();
        float getGain(){
                updataGain();
                return gain;
        }
        float getFRGain()
        {
                updataFullRecGain();
                return gainFR;
        }
        void updataGain ()
        {
                if ( (TP + FN) != 0)
                        gain = 100 * (double)((int)TP- (int)FP)/(double)(TP + FN);
                else
                        gain = 0;
        }
        void updataFullRecGain()
        {
                if ((TPFullRec + FNFullRec) != 0)
                        gainFR = 100 * (double)((int)TPFullRec - (int)FPFullRec)/(double)(TPFullRec + FNFullRec);
                else
                        gainFR = 0;
        }
        size_t getFP()
        {
                return FP;
        }
        size_t getTN()
        {
                return TN;
        }
        size_t getTP()
        {
                return TP;
        }
        size_t getFN()
        {
                return FN;
        }
        size_t getFP_FullR()
        {
                return FPFullRec;
        }
        size_t getTN_FullR()
        {
                return TNFullRec;
        }
        size_t getTP_FullR()
        {
                return TPFullRec;
        }
        size_t getFN_FullR()
        {
                return FNFullRec ;
        }
        size_t getNumOfAllErrors ()
        {
                return TP+ FN;
        }
        size_t getMaxNumOfErrorsInOneRead()
        {
                return maximumNumOfError;
        }
        size_t getNumInitialN(){
                return numInitialN;
        }
        size_t getNumNewIntroducedN(){
                return numNewIntroducedN;
        }
        void settFP(size_t value)
        {
                FP = value;
        }
        void setTN(size_t value)
        {
                TN = value;
        }
        void setTP(size_t value)
        {
                TP = value;
        }
        void setFN(size_t value)
        {
                FN = value;
        }
        void setFP_FullR(size_t value)
        {
                FPFullRec = value;
        }
        void setTN_FullR(size_t value)
        {
                TNFullRec = value;
        }
        void setTP_FullR(size_t value)
        {
                TPFullRec = value;
        }
        void setFN_FullR(size_t value)
        {
                FNFullRec = value;
        }
        void setNumOfReads(size_t value)
        {
                numReads = value;
        }
        void setNumOfComReads(size_t value)
        {
                numOfcomparedReads = value;
        }
        void setNumOfInitialN(size_t value){
                numInitialN = value;
        }
        void setNumOfNewIntroducedN(size_t value){
                numNewIntroducedN = value;
        }


};

class Comparison{
private :
        string filePreName;
        size_t maxIndel;
        NWAligner alignment;
        std::mutex &writeInFileMutexLocal;
        size_t maxAlignmentRound;
        vector<string> splitString(const string str,const char delimiter);
        bool checkHeading(ReadRecord erroneousRecord,ReadRecord perfectRecord,ReadRecord correctedRecord);
        bool compareRead(ReadRecord erroneousRecord,ReadRecord perfectRecord,ReadRecord correctedRecord);

        bool checkIndelSize(const string s1, const string s2);
        /**
         * (thread-safe) method
         * write the reads which contains FP
         *
         */
        void writeWorse(ReadRecord erroneousRecord,ReadRecord perfectRecord,
                            ReadRecord correctedRecord, string erroneousRead, string perfectRead, string correctedRead, size_t FP, size_t TP, size_t TN, size_t FN);
        /**
         * shifts the gaps sign '-' to the end of reads if it dosn't change the alignment
         *
         */
        bool moveGapToEnd(string & erroneousRead,string &correctedRead, string &perfectRead );


public:


        ComparisonMetrics metrics;
        void compareChunk(vector<ReadRecord>& readChunkErroneous, vector<ReadRecord>& readChunkPerfect, vector<ReadRecord>& readChunkCorrected);
        /**
         * Default constructor
         * @param writeInFileMutexGlobal global mutex between threads to write in  file
         */
        Comparison (std::mutex &writeInFileMutexGlobal,string filePre):filePreName(filePre), maxIndel(3), alignment(maxIndel, 1, -1, -3),writeInFileMutexLocal(writeInFileMutexGlobal),maxAlignmentRound(10) {};
};

class ComparisonHandler {
        size_t numOfThreads;
        std::mutex readMutex;          // mutex for reading from files
        std::mutex writeMutex;         // mutex for writing in file common between all threads
        string filePreName;

        void workerThread(size_t myID, LibraryContainer &libContErroneous ,
                          LibraryContainer&libContPerfect, LibraryContainer &libContCorrected);


        void synchroniseReadNumbers(vector<ReadRecord> &myReadBufErroneous, vector<ReadRecord> &myReadBufPerfect, vector<ReadRecord> &myReadBufCorrected,
                                               vector<ReadRecord> & myReadBufErroneous_pre,vector<ReadRecord> & myReadBufPerfect_pre ,vector<ReadRecord>& myReadBufCorrected_pre);

        /**
         * corrected reads which contains FP, are written in these files
         *
         */
        void createEmptyFiles(string filePre ){
                ofstream mapped , perfec, corrected, alignment;
                mapped.open(filePre+"_mapped.fastq");
                mapped.close();
                perfec.open(filePre+"_perfect.fasta");
                perfec.close();
                corrected.open(filePre+"_corrected.fastq");
                corrected.close();
                alignment.open(filePre+"_alignment.dat");
                alignment.close();
        }
public:
        ComparisonMetrics metrics;

        /**
         * Default constructor
         */
        ComparisonHandler(string filePre,  size_t _numOfThreads =1): numOfThreads(_numOfThreads), filePreName(filePre){
                createEmptyFiles(filePre);
        };
        void doComparison(LibraryContainer &libContErroneous ,
                          LibraryContainer&libContPerfect, LibraryContainer &libContCorrected, string filePre);


};
