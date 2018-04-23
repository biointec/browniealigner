#include "comparison.h"


#include <string>
#include <thread>
#include "util.h"
#include <algorithm>
#include <iomanip>



void NWAligner::applyAlignmentToArgumentsBanded(string& s1, string& s2){
                const NWAligner& F = *this;   // shorthand notation

        string al1, al2;

        int i = s1.size();
        int j = s2.size();
        while (i > 0 || j > 0) {
                if ((i > 0) && (j > 0) && (F(i, j) == F(i-1, j-1) + S(s1[i-1], s2[j-1]))) {
                        al1.push_back(s1[i-1]);
                        al2.push_back(s2[j-1]);
                        i--;
                        j--;
                } else if ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G)) {
                        al1.push_back(s1[i-1]);
                        al2.push_back('-');
                        i--;
                } else {
                        al1.push_back('-');
                        al2.push_back(s2[j-1]);
                        j--;
                }
        }

        reverse(al1.begin(), al1.end());
        reverse(al2.begin(), al2.end());
        s1 = al1;
        s2 = al2;

}
void ComparisonMetrics::addMetrics( ComparisonMetrics& tMetrics)
{
        lock_guard<mutex> lock(metricMutex);
        numReads = numReads + tMetrics.numReads ;
        setNumOfComReads( getNumOfCompReads()+ tMetrics.getNumOfCompReads());
        setTP (getTP()+ tMetrics.getTP());
        setTN (getTN()+ tMetrics.getTN());
        settFP(getFP()+ tMetrics.getFP());
        setFN (getFN()+ tMetrics.getFN());
        setTP_FullR( getTP_FullR()+ tMetrics.getTP_FullR());
        setFP_FullR( getFP_FullR()+ tMetrics.getFP_FullR());
        setFN_FullR( getFN_FullR()+ tMetrics.getFN_FullR());
        setTN_FullR( getTN_FullR()+ tMetrics.getTN_FullR());


        setNumOfInitialN(getNumInitialN() + tMetrics.getNumInitialN());
        setNumOfNewIntroducedN(getNumNewIntroducedN() + tMetrics.getNumNewIntroducedN());
        if (getMaxNumOfErrorsInOneRead() < tMetrics.getMaxNumOfErrorsInOneRead())
                maximumNumOfError = tMetrics.getMaxNumOfErrorsInOneRead();
}

void ComparisonMetrics::printStatistics()
{
        std::cout <<std::setprecision(2)  << std::fixed;
        cout << endl << "<<<Report for reads>>>" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << "Number of reads for comparison: " <<getNumOfCompReads()<< " , "<<getNumReads() -getNumOfCompReads() << " reads were skipped because they don't meet the requirement."  << endl;
        cout << "Maximum number Of Errors in one read is: " << getMaxNumOfErrorsInOneRead();
        cout << "\n<<<Alignment based report>>>\n";
        cout << "----------------------------------------------------" << endl;
        cout << "The total number of errors in the uncorrected data: "<<getNumOfAllErrors()<<", (" <<(double)getNumOfAllErrors()/(double)getNumOfCompReads() <<" per read)"<<endl;
        cout << "The number of initial N in the uncorrected data: " << getNumInitialN() <<", (" <<(double)getNumInitialN() /(double)getNumOfCompReads() <<" per read)"<<endl;
        cout << "The number of newly introduced N in the corrected data " << getNumNewIntroducedN() <<", (" <<(double)getNumNewIntroducedN() /(double)getNumOfCompReads() <<" per read)"<<endl;
        cout << "\n<<<The evaluation report based on base pairs>>>\n";
        cout << "----------------------------------------------------" << endl;
        cout << "Among " << getNumOfAllErrors() <<" of Errors "<< getTP() << " number of them are corrected "<<endl;
        cout << "     TP:"<< getTP() <<"      TN:"<< getTN() <<"      FP:"<< getFP() <<"    FN:"<<getFN() <<endl;
        cout << "    The Gain value is: ("<<std::setprecision(2)<<getGain()<<"%)" <<endl;

        cout << "\n<<<The evaluation report based full read recovery>>>\n";
        cout << "----------------------------------------------------" << endl;
        cout << "Among "<<getNumOfCompReads()<<" of reads "<<(getTP_FullR()+ getTN_FullR())<<", number of them are fully recovered " <<endl;

        cout << "     TP: "<< getTP_FullR() <<"      TN: "<< getTN_FullR() <<"      FP: "<< getFP_FullR()<<"     FN: "<<getFN_FullR()<<endl;
        cout << "    The Gain value for full recovery of reads is: ("<< getFRGain() <<"%)" <<endl;
        cout<< endl<< "<<<Quality based reports>>>"<<endl;

}
void ComparisonMetrics::compareRead(string  correctedRead , string  erroneousRead, string   perfectRead, bool aligned)
{

        size_t initialErrors = 0, remainingErrors = 0;
        char correctChar , modifiedChar, erroneousChar;

        if (correctedRead.size() != erroneousRead.size())
                return;
        if (correctedRead.size() != perfectRead.size())
                return;
        setNumOfComReads(getNumOfCompReads() + 1);
       
       
        size_t len = perfectRead.size();
        size_t index  = 0 ;

        //indels at the begining and at the ends are ignored

        while ((index < len) && (correctedRead[index] == '-' || erroneousRead[index]  == '-'  ||  perfectRead[index]  == '-' ))
                index ++;
        while ((len > index) && (correctedRead[len-1] == '-' || erroneousRead[len -1]   == '-'  || perfectRead[len -1]  == '-' ))
                len --;
        while (index < len){
                modifiedChar = correctedRead[index];
                erroneousChar = erroneousRead[index];
                correctChar = perfectRead[index];
                if (modifiedChar == 'N' || erroneousChar=='N'|| correctChar=='N'){
                        index ++;
                        if (erroneousChar == 'N')
                                setNumOfInitialN(getNumInitialN()+1);
                        else{
                                if (modifiedChar == 'N' && correctChar!='N')
                                        setNumOfNewIntroducedN(getNumNewIntroducedN() +1);
                        }
                        continue;

                }
                if (erroneousChar != correctChar) {
                        initialErrors++;
                        if (modifiedChar == correctChar) {
                                TP ++;
                        }
                        else {
                                FN ++;
                                remainingErrors++;
                        }
                } else {
                        if(modifiedChar == correctChar) {
                                TN ++;
                        }
                        else {
                                FP ++;
                                remainingErrors++;
                        }
                }
                index ++;
        }

        if(initialErrors > maximumNumOfError)
                maximumNumOfError = initialErrors;
        
        // if the read is not aligned, then it is FN
        if (!aligned){
                FNFullRec++;
                return;
        }
        
        //full recovery updataSta
        if (remainingErrors != 0){

                if (initialErrors == 0)

                        FPFullRec ++;
                else

                        FNFullRec ++;
        }
        else {

                if (initialErrors > 0)
                        TPFullRec ++;
                else
                        
                        TNFullRec ++;
                                
        }
}

void Comparison::writeWorse(ReadRecord erroneousRecord,ReadRecord perfectRecord,
                            ReadRecord correctedRecord, string erroneousRead, string perfectRead, string correctedRead, size_t FP, size_t TP, size_t TN, size_t FN)
{
         // lock the mutex
        std::unique_lock<std::mutex> lock(writeInFileMutexLocal);
        ofstream mapped(filePreName+"_mapped.fastq",std::ofstream::app);
        ofstream perfec(filePreName+"_perfect.fasta",std::ofstream::app);
        ofstream corrected(filePreName+"_corrected.fastq",std::ofstream::app);
        ofstream alignment(filePreName+"_alignment.dat",std::ofstream::app);
        mapped<<erroneousRecord.preRead;
        mapped<<erroneousRecord.read;
        mapped<<erroneousRecord.postRead;
        mapped.close();
        perfec<<perfectRecord.preRead;
        perfec<<perfectRecord.read;
        perfec<<perfectRecord.postRead;
        perfec.close();
        corrected<<correctedRecord.preRead;
        corrected<<correctedRecord.read;
        corrected<<correctedRecord.postRead;
        corrected.close();
        alignment<<endl<< erroneousRecord.preRead;
        alignment<< "\tTP: " <<TP << "\t FP: " <<FP  <<"\tTN: "<<TN << "\tFN: " <<FN <<endl;
        alignment<< "E:" <<erroneousRead <<endl;
        alignment<< "P:" <<perfectRead<<endl;
        alignment<< "C:" <<correctedRead <<endl;
        string alignmentStr = "A:";
        for (size_t i = 0; i < perfectRead.size(); i++){
                if (erroneousRead[i] == perfectRead [i] &&  erroneousRead[i]==correctedRead[i])
                        alignmentStr = alignmentStr + "|" ;

                else
                        alignmentStr = alignmentStr + "*";
        }
        alignment << alignmentStr <<endl;
        alignment.close();
        lock.unlock();

}
bool Comparison::checkIndelSize(const string s1, const string s2)

{
        size_t maxSizeDiff1 = s1.size() > s2.size() ?
        s1.size()- s2.size():  s2.size() - s1.size();

        if (maxSizeDiff1 > maxIndel )
                return false;
        size_t g_1 = std::count(s1.begin(), s1.end(), '-');
        size_t g_2 = std::count(s2.begin(), s2.end(), '-');
        if (g_1 > maxIndel || g_2 >maxIndel)
                return false;
        return true;
}
vector< string > Comparison::splitString(const string str, const char delimiter)
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
bool Comparison::checkHeading(ReadRecord erroneousRecord,ReadRecord perfectRecord,ReadRecord correctedRecord){
       vector<string> erroneousHeading = splitString(erroneousRecord.preRead, '|');
       vector<string> perfectHeading = splitString(perfectRecord.preRead, '|');
       vector<string> correctedHeading = splitString(correctedRecord.preRead, '|');
       string erroneousHeader = splitString( splitString(erroneousHeading[4],'\n')[0],'-')[1];
       string perfectHeader   = splitString( splitString(perfectHeading[8],'\t')[0], '-')[1];
       string correctedHeader = splitString (splitString( correctedHeading[4],'-')[1], '\n')[0];

       if (perfectHeader != erroneousHeader || perfectHeader!=correctedHeader || correctedHeader!= erroneousHeader){
               cout << erroneousHeader <<endl;
               cout << perfectHeader <<endl;
               cout << correctedHeader <<endl;
               return false;
       }
       return true;
}

bool Comparison::compareRead(ReadRecord erroneousRecord, ReadRecord perfectRecord, ReadRecord correctedRecord)
{
        string perfectRead = perfectRecord.read;
        string correctedRead = correctedRecord.read;
        string erroneousRead = erroneousRecord.read;
        
        size_t first = correctedRecord.preRead.find("A(") + 2;
        size_t last = correctedRecord.preRead.find(")A");
        bool aligned = true;
        if (first!=std::string::npos && last!=std::string::npos)
                aligned = correctedRecord.preRead.substr (first,last-first) == "0" ? false : true;
        
        //checkHeading(erroneousRecord, perfectRecord, correctedRecord);
        perfectRead.erase(std::remove(perfectRead.begin(), perfectRead.end(), '-'), perfectRead.end());

        string pre_erroneous = erroneousRead;
        string pre_corrected = correctedRead;
        string pre_perfect = perfectRead;

        if (!checkIndelSize(perfectRead, correctedRead))
                return false;
        if (perfectRead != correctedRead){

                alignment.alignBanded(perfectRead,correctedRead);
                alignment.applyAlignmentToArgumentsBanded(perfectRead,correctedRead);

        }
        if (!checkIndelSize(perfectRead, erroneousRead))
                return false;
        if (perfectRead != erroneousRead){

                alignment.alignBanded(perfectRead, erroneousRead);
                alignment.applyAlignmentToArgumentsBanded(perfectRead, erroneousRead);
        }
        bool hasGap = false;
        if (perfectRead.find('-')!=string::npos || correctedRead.find('-')!=string::npos ||erroneousRead.find('-')!=string::npos)
                hasGap = true;
        size_t round = 0;
        while (hasGap && (pre_perfect != perfectRead || erroneousRead != pre_erroneous || pre_corrected !=correctedRead) ){
                pre_perfect = perfectRead;
                pre_corrected = correctedRead;
                pre_erroneous = erroneousRead;
                if (!checkIndelSize(perfectRead, correctedRead))
                        return false;
                alignment.alignBanded(perfectRead,correctedRead);
                alignment.applyAlignmentToArgumentsBanded(perfectRead,correctedRead);
                if (!checkIndelSize(perfectRead, erroneousRead))
                        return false;
                alignment.alignBanded(perfectRead, erroneousRead);
                alignment.applyAlignmentToArgumentsBanded(perfectRead, erroneousRead);
                round ++;
                if (round > maxAlignmentRound)
                        return false;
        }
        if (hasGap)
                if (!moveGapToEnd(erroneousRead,correctedRead,perfectRead))
                        return false;

        size_t fp = metrics.getFP();
        size_t tp = metrics.getTP();
        size_t fn = metrics.getFN();
        size_t tn = metrics.getTN();
        metrics.compareRead(correctedRead ,erroneousRead,perfectRead, aligned);
        if (metrics.getFP() > fp )
                writeWorse(erroneousRecord,perfectRecord,correctedRecord,erroneousRead,perfectRead,correctedRead,metrics.getFP() - fp,metrics.getTP() - tp, metrics.getTN()-tn, metrics.getFN() -fn );
        return true;
}
bool Comparison::moveGapToEnd(string & erroneousRead,string &correctedRead, string &perfectRead )
{
        size_t maxSize = max(erroneousRead.size(),correctedRead.size() );
        maxSize = max(maxSize, perfectRead.size());
        while(erroneousRead.size() < maxSize)
                erroneousRead = erroneousRead + "-";
        while(correctedRead.size() < maxSize)
                correctedRead = correctedRead + "-";
        while (perfectRead.size() < perfectRead.size())
                perfectRead = perfectRead + "-";
        size_t i = perfectRead.length() / 2;
        size_t resetNum =0;
        size_t maxGapCorrection = 100;
        while (i > 0){
                if ( perfectRead [i] == '-'  ) {
                        if (perfectRead [i-1] == correctedRead [i] && perfectRead [i-1] == erroneousRead[i]){
                                perfectRead [i] = perfectRead [i-1];
                                perfectRead [i-1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }

                }
                if ( erroneousRead[i] == '-'){
                        if (erroneousRead [i-1] == correctedRead [i] && erroneousRead [i-1] == perfectRead[i]){
                                erroneousRead [i] = erroneousRead [i-1];
                                erroneousRead [i-1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }

                }
                if (correctedRead [i] =='-'){
                        if (correctedRead [i-1] == erroneousRead [i] && correctedRead [i-1] == perfectRead[i]){
                                correctedRead [i] = correctedRead [i-1];
                                correctedRead [i-1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }

                }

                if (perfectRead [i] == '-' && erroneousRead[i] == '-') {
                        if (perfectRead[i-1] == correctedRead[i] ){//perfectRead[i-1] == erroneousRead[i-1] 
                                perfectRead [i] = perfectRead [i-1];
                                perfectRead [i-1] = '-';
                                erroneousRead[i] = erroneousRead[i-1];
                                erroneousRead[i-1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }
                }

                if (correctedRead [i] == '-' && erroneousRead[i] == '-') {
                        if (correctedRead[i-1] == erroneousRead[i-1] && correctedRead[i-1] == perfectRead[i] ){
                                correctedRead [i] = correctedRead [i-1];
                                correctedRead [i-1] = '-';
                                erroneousRead[i] = erroneousRead[i-1];
                                erroneousRead[i-1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }
                }

                if (correctedRead [i] == '-' && perfectRead[i] == '-') {
                        if (correctedRead[i-1] == perfectRead[i-1] && correctedRead[i-1] == erroneousRead[i] ){
                                correctedRead [i] = correctedRead [i-1];
                                correctedRead [i-1] = '-';
                                perfectRead[i] = perfectRead[i-1];
                                perfectRead[i-1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }
                }
                if (i>1 && correctedRead[i]=='-' && correctedRead[i-1] =='-'){
                    if (correctedRead[i-2] == perfectRead[i] ){ //&& correctedRead[i-2]== erroneousRead[i]
                        correctedRead[i] = correctedRead[i-2];
                        correctedRead[i-2] = '-';
                        i = perfectRead.length() / 2;
                        resetNum ++;
                    }
                }
                if (resetNum > maxGapCorrection)
                        return false;

                i --;
        }
        resetNum = 0;
        i = perfectRead.length() / 2;
        while (i < perfectRead.length()){
                if ( perfectRead [i] == '-'  ) {
                        if (perfectRead [i+1] == correctedRead [i] && perfectRead [i+1] == erroneousRead[i]){
                                perfectRead [i] = perfectRead [i+1];
                                perfectRead [i+1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }

                }
                if ( erroneousRead[i] == '-'){
                        if (erroneousRead [i+1] == correctedRead [i] && erroneousRead [i+1] == perfectRead[i]){
                                erroneousRead [i] = erroneousRead [i+1];
                                erroneousRead [i+1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }

                }
                if (correctedRead [i] =='-'){
                        if (correctedRead [i+1] == erroneousRead [i] && correctedRead [i+1] == perfectRead[i]){
                                correctedRead [i] = correctedRead [i+1];
                                correctedRead [i+1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }

                }

                if (perfectRead [i] == '-' && erroneousRead[i] == '-') {
                        if (perfectRead[i+1] == correctedRead[i] ){ //perfectRead[i+1] == erroneousRead[i+1] && 
                                perfectRead [i] = perfectRead [i+1];
                                perfectRead [i+1] = '-';
                                erroneousRead[i] = erroneousRead[i+1];
                                erroneousRead[i+1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }
                }

                if (correctedRead [i] == '-' && erroneousRead[i] == '-') {
                        if (correctedRead[i+1] == erroneousRead[i+1] && correctedRead[i+1] == perfectRead[i] ){
                                correctedRead [i] = correctedRead [i+1];
                                correctedRead [i+1] = '-';
                                erroneousRead[i] = erroneousRead[i+1];
                                erroneousRead[i+1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }
                }

                if (correctedRead [i] == '-' && perfectRead[i] == '-') {
                        if (correctedRead[i+1] == perfectRead[i+1] && correctedRead[i+1] == erroneousRead[i] ){
                                correctedRead [i] = correctedRead [i+1];
                                correctedRead [i+1] = '-';
                                perfectRead[i] = perfectRead[i+1];
                                perfectRead[i+1] = '-';
                                i = perfectRead.length() / 2;
                                resetNum ++;
                        }
                }
                if (i <correctedRead.length() -2 && correctedRead[i] == '-' && correctedRead [i+1] =='-'){
                    if (correctedRead[i+2] ==perfectRead[i] ){ //&&correctedRead[i+2]==erroneousRead[i]
                        correctedRead [i] = correctedRead[i+2];
                        correctedRead [i+2] = '-';
                        i = perfectRead.length() / 2;
                        resetNum ++;
                        
                    }
                }
                
                
                
                if (resetNum > maxGapCorrection)
                        return false;

                i ++;
        }
        return true;
}


void Comparison::compareChunk(vector<ReadRecord>& readChunkErroneous, vector<ReadRecord>& readChunkPerfect,
                              vector<ReadRecord>& readChunkCorrected)
{


        size_t minSize = min(readChunkErroneous.size(), readChunkPerfect.size());
        minSize = min (minSize,readChunkCorrected.size());
        for (size_t i =0; i< readChunkErroneous.size(); i++){
               //checkHeading(readChunkErroneous[i], readChunkPerfect[i], readChunkCorrected[i]);
               compareRead(readChunkErroneous[i],readChunkPerfect[i], readChunkCorrected[i]);
        }
        cout << std::fixed << std::setprecision(2) << "\t\t\t\t\t\t Gain (" << metrics.getGain() << "%),\t\t PF: "<<metrics.getFP()<< "\r";
        cout.flush();

        metrics.setNumOfReads(metrics.getNumReads() + readChunkErroneous.size());

}
void ComparisonHandler::synchroniseReadNumbers(vector<ReadRecord> &myReadBufErroneous, vector<ReadRecord> &myReadBufPerfect, vector<ReadRecord> &myReadBufCorrected,
                                               vector<ReadRecord> & myReadBufErroneous_pre,vector<ReadRecord> & myReadBufPerfect_pre ,vector<ReadRecord>& myReadBufCorrected_pre){


        myReadBufErroneous.insert(myReadBufErroneous.begin(),myReadBufErroneous_pre.begin(),myReadBufErroneous_pre.end());
        myReadBufPerfect.insert(myReadBufPerfect.begin(),myReadBufPerfect_pre.begin(),myReadBufPerfect_pre.end());
        myReadBufCorrected.insert(myReadBufCorrected.begin(),myReadBufCorrected_pre.begin(),myReadBufCorrected_pre.end());

        size_t min_size = min(myReadBufErroneous.size(),myReadBufPerfect.size());
        min_size = min (min_size, myReadBufCorrected.size());


        vector<ReadRecord> myReadBufErroneous_original = myReadBufErroneous;
        vector<ReadRecord> myReadBufCorrected_original = myReadBufCorrected;
        vector<ReadRecord> myReadBufPerfect_original   = myReadBufPerfect;

        myReadBufCorrected.clear();
        myReadBufPerfect.clear();
        myReadBufErroneous.clear();
        myReadBufCorrected_pre.clear();
        myReadBufErroneous_pre.clear();
        myReadBufPerfect_pre.clear();

        myReadBufErroneous.reserve(min_size);
        myReadBufCorrected.reserve(min_size);
        myReadBufPerfect.reserve(min_size);
        myReadBufErroneous_pre.reserve(myReadBufErroneous_original.size()-min_size);
        myReadBufCorrected_pre.reserve(myReadBufCorrected_original.size()-min_size);
        myReadBufPerfect_pre.reserve(myReadBufPerfect_original.size()-min_size);

        for (size_t i =0; i <myReadBufErroneous_original.size(); i++ ){
                if (i <min_size)
                        myReadBufErroneous.push_back(myReadBufErroneous_original[i]);

                else
                        myReadBufErroneous_pre.push_back(myReadBufErroneous_original[i]);
        }

        for (size_t i =0; i < myReadBufCorrected_original.size(); i++ ){
                if (i <min_size)
                        myReadBufCorrected.push_back(myReadBufCorrected_original[i]);

                else
                        myReadBufCorrected_pre.push_back(myReadBufCorrected_original[i]);
        }

        for (size_t i =0; i < myReadBufPerfect_original.size(); i++ ){
                if (i <min_size)
                        myReadBufPerfect.push_back(myReadBufPerfect_original[i]);

                else
                        myReadBufPerfect_pre.push_back(myReadBufPerfect_original[i]);
        }


}

void ComparisonHandler::workerThread(size_t myID, LibraryContainer &libContErroneous ,
                          LibraryContainer&libContPerfect, LibraryContainer &libContCorrected)
{
        Comparison copmare(writeMutex, filePreName);

        // local storage of reads
        vector<ReadRecord> myReadBufErroneous, myReadBufPerfect, myReadBufCorrected;
        vector<ReadRecord> myReadBufErroneous_pre, myReadBufPerfect_pre, myReadBufCorrected_pre;
        // performance counters per thread
        bool resulte = true, resultc = true, resultp = true;
        size_t maxBufferSize = 1000;
        bool sync = false;
        while (true ) {
                // lock the mutex
                std::unique_lock<std::mutex> lock(readMutex);
                size_t blockID, recordID;
                if (myReadBufErroneous_pre.size() < maxBufferSize)
                        resulte = libContErroneous.getRecordChunk(myReadBufErroneous, blockID, recordID);
                else{
                        myReadBufErroneous = myReadBufErroneous_pre;
                        myReadBufErroneous_pre.clear();
                }
                if (myReadBufPerfect_pre.size() < maxBufferSize)
                        resultp = libContPerfect.getRecordChunk(myReadBufPerfect, blockID, recordID);
                else
                {
                        myReadBufPerfect = myReadBufPerfect_pre;
                        myReadBufPerfect_pre.clear();
                }
                if (myReadBufCorrected_pre.size() < maxBufferSize)
                        resultc = libContCorrected.getRecordChunk(myReadBufCorrected, blockID, recordID);
                else
                {
                        myReadBufCorrected = myReadBufCorrected_pre;
                        myReadBufCorrected_pre.clear();
                }

                maxBufferSize = max(myReadBufCorrected.size(), myReadBufPerfect.size());
                maxBufferSize = max (maxBufferSize,myReadBufErroneous.size());
                if (sync || myReadBufCorrected.size() != myReadBufErroneous.size() ||
                        myReadBufCorrected.size() != myReadBufPerfect.size() ||
                        myReadBufErroneous.size() != myReadBufPerfect.size()){
                        synchroniseReadNumbers(myReadBufErroneous,myReadBufPerfect,myReadBufCorrected, myReadBufErroneous_pre,myReadBufPerfect_pre ,myReadBufCorrected_pre);
                                sync = true;
                        }
                lock.unlock();
                if (resultc && resulte && resultp)
                        copmare.compareChunk(myReadBufErroneous, myReadBufPerfect,myReadBufCorrected);
                else{
                        myReadBufCorrected.insert(myReadBufCorrected.begin(), myReadBufCorrected_pre.begin(), myReadBufCorrected_pre.end());
                        myReadBufErroneous.insert(myReadBufErroneous.begin(), myReadBufErroneous_pre.begin(), myReadBufErroneous_pre.end());
                        myReadBufPerfect.insert(myReadBufPerfect.begin(), myReadBufPerfect_pre.begin(), myReadBufPerfect_pre.end());
                        copmare.compareChunk(myReadBufErroneous, myReadBufPerfect,myReadBufCorrected);
                        break;
                }

        }

        // update the global metrics with the thread info (thread-safe)
        metrics.addMetrics(copmare.metrics);

}


void ComparisonHandler::doComparison(LibraryContainer &libContErroneous ,
                          LibraryContainer&libContPerfect, LibraryContainer &libContCorrected, string filePre)
{

        // cout << "Number of threads: " << numOfThreads << endl;
        size_t threadSize = 100000;
        libContErroneous.startIOThreads(threadSize,
                                 10 * threadSize * numOfThreads,
                                 false);
        libContPerfect.startIOThreads(threadSize,
                                 10 * threadSize * numOfThreads,
                                 false);
        libContCorrected.startIOThreads(threadSize,
                                 10 * threadSize * numOfThreads,
                                 false);
        // start worker threads
        vector<thread> workerThreads(numOfThreads);
        for (size_t i = 0; i < workerThreads.size(); i++){
                workerThreads[i] = thread(&ComparisonHandler::workerThread,
                                          this, i, ref(libContErroneous),ref(libContPerfect),ref(libContCorrected));
                
        }
        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
                
        metrics.printStatistics();

        libContErroneous.joinIOThreads();
        libContCorrected.joinIOThreads();
        libContPerfect.joinIOThreads();

}




int main(int argc, char ** args)
{
        cout << "Welcome to the comparison software\n" << endl;
        cout << "Today is " << Util::getDateTime() << endl;
        cout << "----------------------------------------------------\n" << endl;
        if (argc < 4) {
                cout << "ERROR: Wrong number of arguments!\n";
                cout << "Usage: ./comparison erroneousReads.fq perfectReads.fa correctedReads.fq" << endl << flush;
                cout << "\n\t\t******ATTENTION******"<<endl;
                cout << "Reads in these files should be in the same order" <<endl;
                cout << "Files should have a same number of reads" <<endl;
                cout << "Maximum number of allowed indels is fixed to 3\n" <<endl;
                return 0;
        }
        string erroneousReadsFileName = args[1];
        string perfectReadsFileName = args[2];
        string correctedReadsFileName = args[3];
        string filePre = "";
        if (argc == 5)
                filePre = args[4];
        size_t numberThreads = 1;

        
        vector<pair<string, string> > libraries_erroneous;

        LibraryContainer libContErroneous;
        libContErroneous.insert(ReadLibrary(erroneousReadsFileName,"","."));
        LibraryContainer libContPerfect;
        libContPerfect.insert(ReadLibrary(perfectReadsFileName,"","."));
        LibraryContainer libContCorrected;
        libContCorrected.insert(ReadLibrary(correctedReadsFileName,"","."));


        Util::startChrono();
        ComparisonHandler cmh(filePre, numberThreads);
        cmh.doComparison(libContErroneous , libContPerfect, libContCorrected,filePre);

        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        cout << "Exiting... bye!" << endl << endl;
        return 0;

}
