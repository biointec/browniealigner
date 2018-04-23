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

#include "library.h"
#include "tkmer.h"

#include "readfile/fastafile.h"
#include "readfile/fastqfile.h"
#include "readfile/rawfile.h"
#include "readfile/sequencefile.h"
#include "readfile/samfile.h"

#include <iostream>
#include <algorithm>
#include <functional>

using namespace std;

// ============================================================================
// FILETYPE ENUM
// ============================================================================

std::ostream &operator<<(std::ostream &out, const FileType &fileType) {
        switch (fileType) {
                case FASTA :
                        out << "fasta";
                        break;
                case FASTA_GZ :
                        out << "fasta.gz";
                        break;
                case FASTQ :
                        out << "fastq";
                        break;
                case FASTQ_GZ :
                        out << "fastq.gz";
                        break;
                case SAM:
                        out << "sam";
                        break;
                case SAM_GZ:
                        out << "sam.gz";
                        break;
                case RAW:
                        out << "raw";
                        break;
                case RAW_GZ:
                        out << "raw.gz";
                        break;
                case UNKNOWN_FT:
                        out << "unknown";
                        break;
        }

        return out;
}

void ReadLibrary::readMetadata(const string& path)
{
        ifstream ifs(getMetadataFilename().c_str());
        if (!ifs.good())
                return;
        string line, tmp;
        getline(ifs, line);
        stringstream(line) >> tmp >> tmp >> tmp >> numReads;
        getline(ifs, line);
        stringstream(line) >> tmp >> tmp >> tmp >> avgReadLength;
        ifs.close();
}

void ReadLibrary::writeMetadata(const string& path) const
{
        ofstream ofs(getMetadataFilename().c_str());
        ofs << "Number of reads: " << numReads << "\n"
            << "Average read length: " << avgReadLength << endl;
        ofs.close();
}

// ============================================================================
// READLIBRARY CLASS
// ============================================================================

ReadLibrary::ReadLibrary(const std::string& inputFilename_,
                         const std::string& outputFilename_,
                         const std::string& tempDir_) :
        inputFilename(inputFilename_), outputFilename(outputFilename_),
        tempDir(tempDir_), fileType(UNKNOWN_FT),
        numReads(0), avgReadLength(0.0)
{
        // try to figure out the file format based on the extension
        string extension;

        if (inputFilename.length() >= 4)
                extension = inputFilename.substr(inputFilename.length() - 4);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".SAM") {
                fileType = SAM;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 4);
        } else if (extension == ".RAW") {
                fileType = RAW;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 4);
        }

        if (inputFilename.length() >= 6)
                extension = inputFilename.substr(inputFilename.length() - 6);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".FASTA") {
                fileType = FASTA;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 6);
        } else if (extension == ".FASTQ") {
                fileType = FASTQ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 6);
        }

        if (inputFilename.length() >= 7)
                extension = inputFilename.substr(inputFilename.length() - 7);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".RAW.GZ") {
                fileType = RAW_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 7);
        } else if (extension == ".SAM.GZ") {
                fileType = SAM_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 7);
        }

        if (inputFilename.length() >= 9)
                extension = inputFilename.substr(inputFilename.length() - 9);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".FASTA.GZ") {
                fileType = FASTA_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 9);
        } else if (extension == ".FASTQ.GZ") {
                fileType = FASTQ_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 9);
        }

        if (fileType == UNKNOWN_FT) {
                cerr << "Brownie: don't know how to open file: '" << inputFilename << "'\n";
                cerr << "Expected one of the following extensions: .fastq, .fasta, .sam, .raw (or .gz variants thereof)\n";
                exit(EXIT_FAILURE);
        }

        if (!Util::fileExists(inputFilename)) {
                cerr << "Brownie: cannot open input read file: '" << inputFilename << "'\n";
                cerr << "Exiting... " << endl;
                exit(EXIT_FAILURE);
        }

        // strip path from base filename
        size_t last_slash_idx = baseFilename.find_last_of("\\/");
        if (last_slash_idx != std::string::npos)
                baseFilename.erase(0, last_slash_idx + 1);

        // set the outputFilename only if not specified by the user
        if (outputFilename.empty()) {
                ostringstream oss;
                oss << baseFilename << ".corr." << fileType;
                outputFilename = oss.str();
        }
}

ReadFile* ReadLibrary::allocateReadFile() const
{
        switch (getFileType()) {
                case FASTA :
                        return new FastAFile(false);
                case FASTA_GZ :
                        return new FastAFile(true);
                case FASTQ :
                        return new FastQFile(false);
                case FASTQ_GZ :
                        return new FastQFile(true);
                case SAM :
                        return new SamFile(false);
                case SAM_GZ :
                        return new SamFile(true);
                case RAW :
                        return new RawFile(false);
                case RAW_GZ :
                        return new RawFile(true);
                case UNKNOWN_FT:
                        assert(false);
                        return NULL;
        }

        assert(false);
        return NULL;
}

// ============================================================================
// RECORDBLOCK CLASS
// ============================================================================

void RecordBlock::getRecordChunk(vector< ReadRecord >& buffer,
                                 size_t& chunkOffset, size_t targetChunkSize)
{
        assert(numChunksRead < numChunks);

        chunkOffset = nextChunkOffset;

        // find out how many records to copy
        size_t thisChunkSize = 0;
        for (size_t i = nextChunkOffset; i < recordBuffer.size(); i++) {
                nextChunkOffset++;
                buffer.push_back(recordBuffer[i]);
                thisChunkSize += recordBuffer[i].getReadLength() + 1 - Kmer::getK();
                if (thisChunkSize >= targetChunkSize)
                        break;
        }

        numChunksRead++;
}

void RecordBlock::getReadChunk(vector< string >& buffer,
                               size_t& chunkOffset, size_t targetChunkSize)
{
        chunkOffset = nextChunkOffset;

        // find out how many records to copy
        size_t thisChunkSize = 0;
        for (size_t i = nextChunkOffset; i < recordBuffer.size(); i++) {
                nextChunkOffset++;
                buffer.push_back(recordBuffer[i].getRead());
                thisChunkSize += recordBuffer[i].getReadLength() + 1 - Kmer::getK();
                if (thisChunkSize >= targetChunkSize)
                        break;
        }

        numChunksRead++;
}

void RecordBlock::writeRecordChunk(const vector< ReadRecord >& buffer,
                                   size_t chunkOffset)
{
        copy(buffer.begin(), buffer.end(), recordBuffer.begin() + chunkOffset);
        numChunksWritten++;
}

// ============================================================================
// LIBRARY CONTAINER CLASS
// ============================================================================

double LibraryContainer::getAvgReadLength() const
{
        size_t sumReadLength = 0;
        size_t sumReads = 0;
        for (const ReadLibrary& lib : container) {
                sumReadLength += lib.getAvgReadLength() * lib.getNumReads();
                sumReads += lib.getNumReads();
        }

        return (double)sumReadLength/(double)sumReads;
}

bool LibraryContainer::getRecordChunk(vector<ReadRecord>& buffer,
                                      size_t& blockID, size_t& recordOffset)
{
        // clear the buffer
        buffer.clear();

        // A) wait until work becomes available
        std::unique_lock<std::mutex> workLock(workMutex);
        workReady.wait(workLock, [this]{return workBlocks.find(currWorkBlockID) !=
                                               workBlocks.end(); });

        // B) get the current working block
        RecordBlock *block = workBlocks[currWorkBlockID];

        // check if it is a termination message (NULL)
        if (block == NULL) {
                workLock.unlock();

                // send a termination message to the output thread
                std::unique_lock<std::mutex> outputLock(outputMutex);
                outputBlocks[currWorkBlockID] = NULL;
                outputLock.unlock();

                // and get out
                return false;
        }

        // if not, get a chunk of work
        blockID = block->getBlockID();
        block->getRecordChunk(buffer, recordOffset, targetChunkSize);

        // C) if you took the last chunk, increase currWorkBlockID
        bool moveToNextBlock = !block->chunkAvailable();
        if (moveToNextBlock)
                currWorkBlockID++;

        // D.a) if no output thread is active, remove block from work queue
        if (moveToNextBlock && !outputThreadActive)
                workBlocks.erase(block->getBlockID());

        workLock.unlock();

        // D.b) if no output thread is active, return block to input queue
        if (moveToNextBlock && !outputThreadActive) {
                block->reset();
                std::unique_lock<std::mutex> inputLock(inputMutex);
                inputBlocks.push_back(block);
                inputReady.notify_one();
                inputLock.unlock();
        }

        assert(!buffer.empty());
        return true;
}

bool LibraryContainer::getReadChunk(vector<string>& buffer,
                                    size_t& blockID, size_t& recordOffset)
{
        // clear the buffer
        buffer.clear();

        // A) wait until work becomes available
        std::unique_lock<std::mutex> workLock(workMutex);
        workReady.wait(workLock, [this]{return workBlocks.find(currWorkBlockID) !=
                                               workBlocks.end(); });

        // B) get the current working block
        RecordBlock *block = workBlocks[currWorkBlockID];

        // check if it is a termination message (NULL)
        if (block == NULL) {
                workLock.unlock();

                // send a termination message to the output thread
                std::unique_lock<std::mutex> outputLock(outputMutex);
                outputBlocks[currWorkBlockID] = NULL;
                outputLock.unlock();

                // and get out
                return false;
        }

        // if not, get a chunk of work
        blockID = block->getBlockID();
        block->getReadChunk(buffer, recordOffset, targetChunkSize);

        // C) if you took the last chunk, increase currWorkBlockID
        bool moveToNextBlock = !block->chunkAvailable();
        if (moveToNextBlock)
                currWorkBlockID++;

        // D.a) if no output thread is active, remove block from work queue
        if (moveToNextBlock && !outputThreadActive)
                workBlocks.erase(block->getBlockID());

        workLock.unlock();

        // D.b) if no output thread is active, return block to input queue
        if (moveToNextBlock && !outputThreadActive) {
                block->reset();
                std::unique_lock<std::mutex> inputLock(inputMutex);
                inputBlocks.push_back(block);
                inputReady.notify_one();
                inputLock.unlock();
        }

        return true;
}


void LibraryContainer::inputThreadLibrary(ReadLibrary& input)
{
        // read counters
        size_t totNumReads = 0, totReadLength = 0;

        // open the read file
        ReadFile *readFile = input.allocateReadFile();
        readFile->open(input.getInputFilename());

        // aux variables
        ReadRecord record;

        while (readFile->good()) {
                // A) wait until a record block is available
                std::unique_lock<std::mutex> inputLock(inputMutex);
                inputReady.wait(inputLock, [this]{return !inputBlocks.empty();});

                // get a record block
                RecordBlock *block = inputBlocks.back();
                inputBlocks.pop_back();

                inputLock.unlock();

                // B) fill up the block (only this thread has access)
                block->reset();

                vector<ReadRecord> &recordBuffer = block->getRecordBuffer();
                size_t thisBlockSize = 0, thisBlockNumRecords = 0;
                size_t thisBlockNumChunks = 0, thisChunkSize = 0;
                while ((thisBlockSize < targetBlockSize) && (readFile->getNextRecord(record))) {
                        thisBlockNumRecords++;
                        totNumReads++;
                        totReadLength += record.getReadLength();

                        recordBuffer.push_back(record);
                        thisBlockSize += record.getReadLength() + 1 - Kmer::getK();

                        // count the number of chunks in this block
                        thisChunkSize += record.getReadLength() + 1 - Kmer::getK();
                        if (thisChunkSize >= targetChunkSize) {
                                thisBlockNumChunks++;
                                thisChunkSize = 0;
                        }
                }

                // also count the final chunk
                if (thisChunkSize > 0)
                        thisBlockNumChunks++;

                block->setBlockInfo(currInputFileID, currInputBlockID++, thisBlockNumChunks);

                // C) push the record block onto the worker queue
                std::unique_lock<std::mutex> workLock(workMutex);
                assert(workBlocks.find(block->getBlockID()) == workBlocks.end());
                workBlocks[block->getBlockID()] = block;

                // notify workers that more work is available
                workReady.notify_all();
                workLock.unlock();

                // D) update stdout information
                cout << "Number of reads processed: " << totNumReads << "\r"; cout.flush();
        }

        cout << "Number of reads processed: " << totNumReads << endl;

        double avgReadLength = (totNumReads == 0) ? 0 :
                (double)totReadLength/(double)totNumReads;
        input.setAvgReadLength(avgReadLength);
        input.setNumReads(totNumReads);

        readFile->close();

        // free temporary memory
        delete readFile;
}

void LibraryContainer::commitRecordChunk(const vector<ReadRecord>& buffer,
                                         size_t blockID, size_t recordOffset)
{
        // A) obtain the work mutex
        std::unique_lock<std::mutex> workLock(workMutex);

        // B) get the corresponding block and write the contents
        assert(workBlocks.find(blockID) != workBlocks.end());
        RecordBlock *block = workBlocks[blockID];
        block->writeRecordChunk(buffer, recordOffset);

        // C.a) if the block is ready, pass it to the output thread
        bool blockReady = block->isProcessed();
        if (blockReady)
                workBlocks.erase(block->getBlockID());

        workLock.unlock();

        // C.b) if the block is ready, pass it to the output thread
        if (blockReady) {
                std::unique_lock<std::mutex> outputLock(outputMutex);
                outputBlocks[block->getBlockID()] = block;

                outputReady.notify_one();
                outputLock.unlock();
        }
}

void LibraryContainer::outputThreadLibrary(ReadLibrary& lib)
{
        // open the read file
        ReadFile *readFile = lib.allocateReadFile();
        readFile->open(lib.getOutputFileName(), WRITE);

        ofstream nodeChainFile(lib.getNodeChainFilename());

        while (true) {
                // A) wait until an output block is ready
                std::unique_lock<std::mutex> outputLock(outputMutex);
                outputReady.wait(outputLock, [this]{return (outputBlocks.find(currOutputBlockID) !=
                                                            outputBlocks.end()); });

                // B) get the next output block
                RecordBlock* block = outputBlocks[currOutputBlockID];
                assert(block == outputBlocks.begin()->second);

                // check if it is a termination message (NULL)
                if (block == NULL) {
                        outputLock.unlock();
                        break;
                }

                // C) if we need to start a new file: get out
                if (block->getFileID() != currOutputFileID) {
                        outputLock.unlock();
                        break;
                }

                // D) remove the block from the output queue
                outputBlocks.erase(currOutputBlockID);
                currOutputBlockID++;

                outputLock.unlock();

                // E) write the actual data (only this thread has access)
                const vector<ReadRecord>& recordBuffer = block->getRecordBuffer();
                for (size_t i = 0; i < recordBuffer.size(); i++) {
                        readFile->writeRecord(recordBuffer[i]);

                        /*const vector<int>& nc = recordBuffer[i].getNodeChain();
                        //if (nc.size() < 3)
                         //       continue;

                        if (nc.empty()) {
                                nodeChainFile << "0\n";
                                continue;
                        }

                        for (size_t j = 0; j < nc.size() - 1; j++)
                                nodeChainFile << nc[j] << " ";
                        nodeChainFile << nc.back() << "\n";*/
                        nodeChainFile << recordBuffer[i].alignmentInfo <<endl;
                        
                }

                if (!readFile->good())
                        throw ios_base::failure("Cannot write to " + lib.getOutputFileName());

                // F) return the block to the input queue
                block->reset();
                std::unique_lock<std::mutex> inputLock(inputMutex);
                inputBlocks.push_back(block);
                inputReady.notify_one();
                inputLock.unlock();
        }

        nodeChainFile.close();
        readFile->close();
        delete readFile;
}

void LibraryContainer::inputThreadEntry()
{
        // read all input data
        for (size_t i = 0; i < getSize(); i++) {
                ReadLibrary &input = getInput(i);

                cout << "Processing file " << i+1 << "/" << getSize() << ": "
                     << input.getInputFilename() << ", type: "
                     << input.getFileType() << endl;

                inputThreadLibrary(input);
                currInputFileID++;
        }

        // send a termination message to the workers
        unique_lock<std::mutex> workLock(workMutex);
        workBlocks[currInputBlockID] = NULL;
        workReady.notify_all();
        workLock.unlock();
}

void LibraryContainer::outputThreadEntry()
{
        // write all output data
        for (size_t i = 0; i < getSize(); i++) {
                ReadLibrary &input = getInput(i);

                outputThreadLibrary(input);
                currOutputFileID++;
        }
}

void LibraryContainer::startIOThreads(size_t targetChunkSize_,
                                      size_t targetBlockSize_,
                                      bool outputThreadActive_)
{
        targetChunkSize = targetChunkSize_;
        targetBlockSize = targetBlockSize_;
        outputThreadActive = outputThreadActive_;

        // initialize input variables
        currInputFileID = currInputBlockID = 0;

        // initialize work thread variables
        currWorkBlockID = 0;

        // initialize output variables
        currOutputFileID = currOutputBlockID = 0;

        // initialize blocks and push them on the input stack
        for (int i = 0; i < NUM_RECORD_BLOCKS; i++)
                inputBlocks.push_back(new RecordBlock());

        // start IO threads
        iThread = thread(&LibraryContainer::inputThreadEntry, this);
        if (outputThreadActive)
                oThread = thread(&LibraryContainer::outputThreadEntry, this);
}

void LibraryContainer::joinIOThreads()
{
        iThread.join();

        if (oThread.joinable())         // output thread might not be active
                oThread.join();

        assert(inputBlocks.size() == NUM_RECORD_BLOCKS);
        // delete blocks and clear the input stack
        for (size_t i = 0; i < inputBlocks.size(); i++)
                delete inputBlocks[i];
        inputBlocks.clear();            // contains all blocks
        workBlocks.clear();             // contains only a termination msg
        outputBlocks.clear();           // contains only a termination msg
}

void LibraryContainer::readMetadata(const string& path)
{
        for (auto& it : container)
                it.readMetadata(path);
}

void LibraryContainer::writeMetadata(const string& path) const
{
        for (const auto& it : container)
                it.writeMetadata(path);
}
