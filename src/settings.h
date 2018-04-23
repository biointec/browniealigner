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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class LibraryContainer;

// ============================================================================
// COMMAND CLASS (ENUM)
// ============================================================================

enum class Command { none, assemble };

// ============================================================================
// SETTINGS CLASS
// ============================================================================

class Settings
{
private:
        /**
         * Print Brownie program version
         */
        void printProgramVersion() const;

        /**
         * Print Brownie usage instructions
         */
        void printUsage() const;

        /**
         * Print Brownie usage instructions for the assemble module
         */
        void printUsageAssemble() const;

        /**
         * Print Brownie usage information for alignment module
         **/
         
        void printUsageAlign()const;
        /**
         * Print Brownie usage instructions for the visualize module
         */
        void printUsageVisualize() const;

        /**
         * Print Brownie usage instructions for the assemble module
         */
        void printUsageCompare() const;

       
        unsigned int kmerSize;          // user specified kmer size
      
        
        std::string pathtotemp;         // directory specified by user
        Command command;                // type of command
        size_t numThreads;              // number of threads
        bool doubleStranded;            // double stranded sequences
        int runSpecificStage;           // takes a non-zero value to run a specific stage
        std::string referenceFilename;  // reference filename
        
        int essaMEMSparsenessFactor;    // sparseness factor for essaMEM
        int bubbleDFSNodeLimit;         // maximum number of visited nodes during bubble detection
        size_t readCorrDFSNodeLimit;    // maximal number of visited nodes during read mapping
        double covCutoff;               // coverage cutoff value to separate true and false nodes based on their node-kmer-coverage
        size_t abundanceMin;            // a kmer should occure this many time to be used in th graph construction. only 1 and 2 are possible
        bool branchAndBound;            // use branch and bound in dfs algorithm.
        bool MarkovFilter;              // align reads in two steps based on the Markov Model
public:
        /**
         * Default constructor
         */
        Settings();
        Settings(unsigned int kmerSize,std::string pathtotemp);

        /**
         * Parse the command line arguments
         * @param argc Command line argument count
         * @param args Command line arguments
         * @param libCont Library container (output)
         */
        void parseCommandLineArguments(int argc, char **args,
                                       LibraryContainer& libCont);

        /**
         * Parse the command line arguments for the assemble module
         * @param argc Command line argument count
         * @param args Command line arguments
         * @param libCont Library container (output)
         */
        void parseCommandLineArgAssemble(int argc, char **args,
                                         LibraryContainer& libCont);

        /**
         * Parse the command line arguments for the vizualize module
         * @param argc Command line argument count
         * @param args Command line arguments
         * @param libCont Library container (output)
         */
        void parseCommandLineArgVisualize(int argc, char **args,
                                          LibraryContainer& libCont);

        /**
         * Parse the command line arguments for the compare module
         * @param argc Command line argument count
         * @param args Command line arguments
         * @param libCont Library container (output)
         */
        void parseCommandLineArgCompare(int argc, char **args,
                                        LibraryContainer& libCont);

        /**
         * Get the type of command
         * @return The type of command
         */
        Command getCommand() const {
                return command;
        }

        /**
         * Get the specific stage to run
         * @return The specific stage to run
         */
        int getRunSpecificStage() const {
                return runSpecificStage;
        }

        /**
         * Get the reference filename
         * @return The reference filename
         */
        std::string getReferenceFilename() const {
                return referenceFilename;
        }

        /**
         * Get the user-specified kmer hash length
         * @return The hash length
         */
        unsigned int getK() const {
                return kmerSize;
        }

        /**
         * Get the number of threads
         * @return The number of threads
         */
        size_t getNumThreads() const {
                return numThreads;
        }

        /**
         * Check of if the reads are double stranded
         * @return true of false
         */
        bool isDoubleStranded() const {
                return doubleStranded;
        }

        /**
         * Get the temporary working directory
         * @return The temporary working directory
         */
        std::string getTempDirectory() const {
                return pathtotemp;
        }

        /**
         * Prepend the temporary working directory to a file
         * @return The expanded directory
         */
        std::string addTempDirectory(const std::string filename) const {
                return pathtotemp + filename;
        }

        /**
         * Get the IO block size in number of kmers
         * @return The IO block size in number of kmers
         */
        size_t getThreadWorkSize() const {
                return 100000;
        }

        /**
         * Get the IO block size in number of kmers
         * @return The IO block size in number of kmers
         */
        size_t getThreadBubbleWorkSize() const {
                return 32768;
        }

        /**
         * Get the essaMEM sparseness factor
         * @return The essaMEM sparseness factor
         */
        int getEssaMEMSparsenessFactor() const {
                return essaMEMSparsenessFactor;
        }

        /**
         * Get the maximum number of nodes visited during a DFS during bubble detection
         * @return The maximum number of nodes visited during bubble detection
         */
        int getBubbleDFSNodeLimit() const {
                return bubbleDFSNodeLimit;
        }

        /**
         * Get the maximum number of nodes visited during a DFS during read correction
         * @return The maximum number of nodes visited during read correction
         */
        size_t getReadCorrDFSNodeLimit() const {
                return readCorrDFSNodeLimit;
        }

        /**
         * Get the coverage cutoff value for a node
         * @return The coverage cutoff value for a node
         */
        double getCutOffValue(){
                return covCutoff;
        }

         /**
         * Get the abundanceMin cutoff value for a node
         * @return The abundanceMinvalue for a node
         */
        double getAbundanceMinValue(){
                return abundanceMin;
        }
        
        /**
         * Get the branchAndBound filter value 
         * @return The branchAndBound value 
         **/
        bool getBranchAndBound() const{
                return branchAndBound;
        }
        /**
         * Get the MarkovFilter value 
         * @return MarkovFilter
         * 
         **/
        bool getMarkovFilter () const{
                return MarkovFilter;
        }
        
};

#endif
