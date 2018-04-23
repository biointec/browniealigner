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

#ifndef READLIBRARY_H
#define READLIBRARY_H

#include "global.h"
#include "readfile/readfile.h"

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include <mutex>
#include <condition_variable>

// ============================================================================
// FILETYPE ENUM
// ============================================================================

typedef enum { FASTQ, FASTA, FASTQ_GZ, FASTA_GZ,
               SAM, SAM_GZ, RAW, RAW_GZ, UNKNOWN_FT} FileType;

std::ostream &operator<<(std::ostream &out, const FileType &fileType);

// ============================================================================
// READ LIBRARY CLASS
// ============================================================================

class ReadLibrary
{
private:
        std::string inputFilename;      // name of the input file
        std::string outputFilename;     // name of the output file
        std::string baseFilename;       // base file name (no extension, no path)
        std::string tempDir;            // path to store temporary files

        FileType fileType;              // file type (FASTQ, FASTA, etc.)


        size_t numReads;                // number of reads in library
        double avgReadLength;           // average length of the reads

public:
        /**
         * Default constructor
         * @param inputFilename Input filename
         * @param outputFilename Output filename
         * @param tempDir Path to store temporary files
         */
        ReadLibrary(const std::string& inputFilename,
                    const std::string& outputFilename,
                    const std::string& tempDir);

        /**
         * Allocate and return the correct read file for a certain input
         * @return An allocated readfile
         */
        ReadFile* allocateReadFile() const;

        /**
         * Get the input filename
         * @return The input filename
         */
        std::string getInputFilename() const {
                return inputFilename;
        }
        /**
         * Get the base filename without extension or directory
         * @return The baseFilename
         */
        std::string getBaseFilename() const {
                return baseFilename;
        }
        /**
         * Get the output filename
         * @return The output filename
         */
        std::string getOutputFileName() const {
                return outputFilename;
        }

        /**
         * Get the node chain filename
         * @return The node chain filename
         */
        std::string getNodeChainFilename() const {
                return tempDir + baseFilename + ".ncf";
        }

        /**
         * Get the meta data filename
         * @return The meta data filename
         */
        std::string getMetadataFilename() const {
                return tempDir + baseFilename + ".met";
        }

        /**
         * Get the file type
         * @return The filetype
         */
        FileType getFileType() const {
                return fileType;
        }

        /**
         * Set the number of reads in this library
         * @return The number of reads in this library
         */
        void setNumReads(size_t target) {
                numReads = target;
        }

        /**
         * Get the number of reads in this library
         * @return The number of reads in this library
         */
        size_t getNumReads() const {
                return numReads;
        }

        /**
         * Set the read length
         * @param target The read length
         */
        void setAvgReadLength(double target) {
                avgReadLength = target;
        }

        /**
         * Get the read length
         * @return The read length
         */
        double getAvgReadLength() const {
                return avgReadLength;
        }

        /**
         * Read the metadata to disk
         * @param path Path to dir where the metadata should be read
         */
        void readMetadata(const std::string& path);

        /**
         * Write the metadata to disk
         * @param path Path to dir where the metadata should be written
         */
        void writeMetadata(const std::string& path) const;
};

// ============================================================================
// RECORDBLOCK CLASS
// ============================================================================

class RecordBlock
{
private:
        size_t fileID;                  // file identifier
        size_t blockID;                 // block identifier
        size_t numChunks;               // number of chunks in this block
        size_t numChunksRead;           // number of chunks read
        size_t numChunksWritten;        // number of chunks written
        size_t nextChunkOffset;         // offset of the next chunk to be read

        std::vector<ReadRecord> recordBuffer;   // actual records

public:
        /**
         * Default constructor
         */
        RecordBlock() : fileID(0), blockID(0), numChunks(0),
                numChunksRead(0), numChunksWritten(0), nextChunkOffset(0) {}

        /**
         * Clear the contents of this block
         */
        void reset() {
                fileID = blockID = numChunks = numChunksRead = 0;
                numChunksWritten = nextChunkOffset = 0;
                recordBuffer.clear();
        }

        /**
         * Return the blockID
         * @return The blockID
         */
        size_t getBlockID() const {
                return blockID;
        }

        /**
         * Return the fileID
         * @return The fileID
         */
        size_t getFileID() const {
                return fileID;
        }

        /**
         * Return true of at least one chunk if available for reading
         * @return true or false
         */
        bool chunkAvailable() const {
                return numChunksRead < numChunks;
        }

        /**
         * Returns true if the block is fully processed
         * @return true of false
         */
        bool isProcessed() const {
                return numChunksWritten == numChunks;
        }

        /**
         * Accessor to the record buffer
         * @return A reference to the record buffer
         */
        std::vector<ReadRecord>& getRecordBuffer() {
                return recordBuffer;
        }

        /**
         * Set the block meta information
         * @param fileID_ File identifier
         * @param blockID_ Block identifier
         * @param numChunks_ Number of chunks in this block
         */
        void setBlockInfo(size_t fileID_, size_t blockID_, size_t numChunks_) {
                fileID = fileID_;
                blockID = blockID_;
                numChunks = numChunks_;
        }

        /**
         * Get next available input record chunk
         * @param buffer Record buffer to write to (contents will be appended)
         * @param chunkOffset Chunk offset (output)
         * @param targetChunkSize Target chunk size (input)
         */
        void getRecordChunk(std::vector<ReadRecord>& buffer,
                            size_t& chunkOffset, size_t targetChunkSize);

        /**
         * Get next available input read chunk
         * @param buffer Record buffer to write to (contents will be appended)
         * @param chunkOffset Chunk offset (output)
         * @param targetChunkSize Target chunk size (input)
         */
        void getReadChunk(std::vector<std::string>& buffer,
                          size_t& chunkOffset, size_t targetChunkSize);

        /**
         * Replace a chunk with the contents provided by the buffer
         * @param buffer Record buffer with new contents
         * @param chunkOffset Record offset of the chunk
         */
        void writeRecordChunk(const std::vector<ReadRecord>& buffer,
                              size_t chunkOffset);
};

// ============================================================================
// LIBRARYCONTAINER CLASS
// ============================================================================

class LibraryContainer
{
private:
        std::vector<ReadLibrary> container;             // all library files

        size_t targetBlockSize;                         // a block is a number of (overlapping) k-mers read in a single run
        size_t targetChunkSize;                         // a chunk is a piece of work attributed to one thread

        std::thread iThread;                            // input thread
        std::thread oThread;                            // output thread
        bool outputThreadActive;                        // is the output thread active?

        // input thread variables (protect by inputMutex)
        std::mutex inputMutex;                          // input thread mutex
        std::condition_variable inputReady;             // input blocks ready condition
        std::vector<RecordBlock*> inputBlocks;          // input blocks that are empty
        size_t currInputFileID;                         // identifier of the current input file
        size_t currInputBlockID;                        // identifier of the current input block

        // worker thread variables (protect by workMutex)
        std::mutex workMutex;                           // worker thread mutex
        std::condition_variable workReady;              // work ready condition
        std::map<size_t, RecordBlock*> workBlocks;      // blocks being processed
        size_t currWorkBlockID;                         // identifier of the block currently being processed

        // output thread variables (protect by outputMutex)
        std::mutex outputMutex;                         // output thread mutex
        std::condition_variable outputReady;            // write buffer full condition
        std::map<size_t, RecordBlock*> outputBlocks;    // processed blocks
        size_t currOutputFileID;                        // identifier of the current output file
        size_t currOutputBlockID;                       // identifier of the current output block

        /**
         * Write the inputs for a specific library
         * @param library Reference to a specific library
         */
        void inputThreadLibrary(ReadLibrary& library);

        /**
         * Write the outputs for a specific library
         * @param library Reference to a specific library
         */
        void outputThreadLibrary(ReadLibrary& library);

        /**
         * Entry routine for the input threads
         */
        void inputThreadEntry();

        /**
         * Entry routine for the output threads
         */
        void outputThreadEntry();

public:
        /**
         * Add a read library to the container
         * @arg library Library to add
         */
        void insert(const ReadLibrary& library) {
                container.push_back(library);
        }

        /**
         * Get the number of inputs
         * @return The number of inputs
         */
        size_t getSize() const {
                return container.size();
        }

        /**
         * Get a specified input
         * @param index Input identifier
         * @return A reference to the input
         */
        ReadLibrary& getInput(size_t index) {
                assert(index < container.size());
                return container[index];
        }

        /**
         * Get a specified input
         * @param index Input identifier
         * @return A reference to the input
         */
        const ReadLibrary& getInput(size_t index) const {
                assert(index < container.size());
                return container[index];
        }

        /**
         * The (weighted) average read length over all libraries
         * @return The (weighted) average read length over all libraries
         */
        double getAvgReadLength() const;

        /**
         * Get next record chunk from the input
         * @param buffer Buffer in which to store the records (output)
         * @param blockID Block identifier (output)
         * @param recordOffset Record offset within the block (output)
         * @return False if no more reads are available, true otherwise
         */
        bool getRecordChunk(std::vector<ReadRecord>& buffer,
                            size_t& blockID, size_t& recordOffset);

        /**
         * Get next read chunk from the input
         * @param buffer Buffer in which to store the reads (output)
         * @param blockID Block identifier (output)
         * @param recordOffset Record offset within the block (output)
         * @return False if no more reads are available, true otherwise
         */
        bool getReadChunk(std::vector<std::string>& buffer,
                          size_t& blockID, size_t& recordOffset);

        /**
         * Commit a ReadRecord chunk to the write buffer
         * @param blockID Block identfier obtained from get...Chunk()
         * @param recordOffset Record offset obtained from get...Chunk()
         */
        void commitRecordChunk(const std::vector<ReadRecord>& buffer,
                               size_t blockID,
                               size_t recordOffset);

        /**
         * Initialize read threading
         * @param targetChunkSize Target size for a single chunk
         * @param targetBlockSize Target size for a single block
         * @param writeReads True if the input reads are again to be written
         */
        void startIOThreads(size_t targetChunkSize,
                            size_t targetBlockSize,
                            bool writeReads = false);

        /**
         * Join input (and optionally) also the output thread
         * Make sure that all worker threads are joined before calling this
         * routine.  Calling getReadChunk() or commitRecordChunk() after this
         * routine is called will result in undefined behavior
         */
        void joinIOThreads();

        /**
         * Read the metadata to disk
         * @param path Path to dir where the metadata should be read
         */
        void readMetadata(const std::string& path);

        /**
         * Write the metadata to disk
         * @param path Path to dir where the metadata should be written
         */
        void writeMetadata(const std::string& path) const;
};

#endif
