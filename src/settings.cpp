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

#include "settings.h"
#include "library.h"

#include <iostream>
#include <fstream>

using namespace std;

// ============================================================================
// SETTINGS CLASS PRIVATE
// ============================================================================

void Settings::printProgramVersion() const
{
        cout << "brownie version " << BROWNIE_MAJOR_VERSION << "."
             << BROWNIE_MINOR_VERSION << "." << BROWNIE_PATCH_LEVEL << "\n";

        cout << "Copyright (C) 2016 Jan Fostier (jan.fostier@ugent.be)\n";
        cout << "This is free software; see the source for copying conditions. "
                "There is NO\nwarranty; not even for MERCHANTABILITY or "
                "FITNESS FOR A PARTICULAR PURPOSE.\n" << endl;

        cout << "Compilation settings:\n";
        cout << "  MAXKMERLENGTH = " << MAXKMERLENGTH << "\n";
#ifdef HAVE_ZLIB
        cout << "  ZLIB support = enabled";
#else
        cout << "  ZLIB support = disabled";
#endif
        cout << endl;
}

void Settings::printUsage() const
{
        cout << "Usage: brownie command [options]\n\n";

        cout << " command\n";
        //cout << "  assemble\t\tassemble short-read data\n";
        //cout << "  visualize\t\tproduce a Cytoscape-compatible de Bruijn graph visualization\n";
        //cout << "  compare\t\tcompare reference sequences to the de Brijn graph\n";
        cout << "  index\t\t build the de Bruijn graph\n"; 
        cout << "  align\t\t align short-read data to the graph\n\n";
        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help page\n";
        cout << "  -v\t--version\tdisplay version\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void Settings::printUsageAssemble() const
{
        cout << "Usage: brownie assemble [options] [file_options] file1 [[file_options] file2]...\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help page\n";
        //cout << "  -s\t--singlestranded\tenable single-stranded DNA [default = false]\n";  // This has never been tested!
        cout << "  -sX\t--stageX\trun only stage X (X = 1,2,...,6)\n\n";

        cout << " [options arg]\n";
        cout << "  -k\t--kmersize\tkmer size [default = 31]\n";
        cout << "  -t\t--threads\tnumber of threads [default = available cores]\n";
        cout << "  -v\t--visits\tsearch depth during bubble detection [default = 1000]\n";
        cout << "  -d\t--depth\t\tsearch depth during read correction [default = 1000]\n";
        cout << "  -e\t--essa\t\tsparseness factor of index structure [default = 1]\n";
        cout << "  -c\t--cutoff\tcoverage threshold during error correction [default = auto]\n";
        cout << "  -p\t--pathtotmp\tpath to directory to store temporary files [default = .]\n\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\tcorrected output read file name [default = inputfile.corr]\n\n";

        cout << " examples:\n";
        cout << "  ./brownie assemble inputA.fastq\n";
        cout << "  ./brownie assemble -k 29 -t 4 -o outputA.fasta inputA.fasta -o outputB.fasta inputB.fastq\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void Settings::printUsageAlign() const
{       
        cout << "Usage: brownie index [options] file1 [file2]...\n\n";
        cout << "Usage: brownie align [options] [file_options] file1 [[file_options] file2]...\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help page\n";
        cout << " [options arg]\n";
        cout << "  -k\t--kmersize\tkmer size [default = 31]\n";
        cout << "  -t\t--threads\tnumber of threads [default = available cores]\n";
        cout << "  -e\t--essa\t\tsparseness factor of index structure [default = 1]\n";
        cout << "  -nBB\t--noBranchAndBound\t do not use branch&Bound pruning\n" ;  
        cout << "  -nMM\t--noMarkovModel\t do not use Markov Model\n" ; 
                
        cout << " [file_options]\n";
        cout << "  -o\t--output\tcorrected output read file name [default = inputfile.corr]\n\n";

        cout << " examples:\n";
        cout << "  ./brownie index -k 31 -t 4 -p tempDir genome.fasta\n\n";
        cout << "  ./brownie align -k 31 -t 4 -p tempDir -o outputB.fastq inputB.fastq\n\n"; 
        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}
void Settings::printUsageVisualize() const
{
        cout << "Usage: brownie visualize [options] output_file\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help page\n";
        cout << "  -sX\t--stageX\tcompare with stage X (X = 3,4) [default = -s4]\n\n";

        cout << " [options arg]\n";
        cout << "  -s\t--seed\t\tseed node ID [default = first valid node]\n";
        cout << "  -d\t--depth\t\tmaximum depth from seed node [default = 10]\n\n";

        cout << " examples:\n";
        cout << "  ./brownie visualize output\n";
        cout << "  ./brownie visualize -s 10 -d 20 output\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void Settings::printUsageCompare() const
{
        cout << "Usage: brownie compare [options] reference.fasta\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help page\n";
        cout << "  -sX\t--stageX\tcompare with stage X (X = 3,4) [default = -s4]\n\n";

        cout << " [options arg]\n";
        cout << "  -k\t--kmersize\tkmer size [default = 31]\n";
        cout << "  -t\t--threads\tnumber of threads [default = available cores]\n";
        cout << "  -v\t--visits\tsearch depth during bubble detection [default = 1000]\n";
        cout << "  -d\t--depth\t\tsearch depth during read correction [default = 1000]\n";
        cout << "  -e\t--essa\t\tsparseness factor of index structure [default = 1]\n";
        cout << "  -c\t--cutoff\tcoverage threshold during error correction [default = auto]\n";
        cout << "  -p\t--pathtotmp\tpath to directory to store temporary files [default = .]\n\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\tcorrected output read file name [default = inputfile.corr]\n\n";

        cout << " examples:\n";
        //cout << "  ./brownie assemble inputA.fastq\n";
        //cout << "  ./brownie assemble -k 31 -t 4 -o outputA.fasta inputA.fasta -o outputB.fasta inputB.fastq\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

// ============================================================================
// SETTINGS CLASS PUBLIC
// ============================================================================

Settings::Settings() :  kmerSize(31),command(Command::none),
        numThreads(std::thread::hardware_concurrency()), doubleStranded(true),
        runSpecificStage(0),
        essaMEMSparsenessFactor(1), bubbleDFSNodeLimit(1000),
        readCorrDFSNodeLimit(5000), covCutoff(0), abundanceMin(2) , branchAndBound(true), MarkovFilter(true), essaMEM(true)  {}
 
Settings::Settings(unsigned int kmerSize,std::string pathtotemp) : kmerSize(kmerSize),pathtotemp(pathtotemp),command(Command::none),
        numThreads(std::thread::hardware_concurrency()), doubleStranded(true),
        runSpecificStage(0),
        essaMEMSparsenessFactor(1), bubbleDFSNodeLimit(1000),
        readCorrDFSNodeLimit(5000), covCutoff(0), abundanceMin(2) ,  branchAndBound(true), MarkovFilter(true), essaMEM(true)   {}
void Settings::parseCommandLineArguments(int argc, char** args,
                                         LibraryContainer& libCont)
{
        if (argc > 1) {
                string arg(args[1]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if ((arg == "-v") || (arg == "--version")) {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else if (arg == "assemble" || arg == "index" || arg == "align") {
                        command = Command::assemble;
                        parseCommandLineArgAssemble(argc, args, libCont);
                } 
        }
}

void Settings::parseCommandLineArgAssemble(int argc, char** args,
                                           LibraryContainer& libCont)
{
        
        // parse all input arguments
        string inputFilename, outputFilename;
        vector<pair<string, string> > libraries;
        for (int i = 2; i < argc; i++) {
                string arg(args[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsageAlign();
                        exit(EXIT_SUCCESS);
                } else if ((arg == "-s1") || (arg == "--stage1")) {
                        runSpecificStage = 1;
                } else if ((arg == "-s2") || (arg == "--stage2")) {
                        runSpecificStage = 2;
                } else if ((arg == "-s3") || (arg == "--stage3")) {
                        runSpecificStage = 3;
                } else if ((arg == "-s4") || (arg == "--stage4")) {
                        runSpecificStage = 4;
                } else if ((arg == "-s5") || (arg == "--stage5")) {
                        runSpecificStage = 5;
                } else if ((arg == "-s6") || (arg == "--stage6")) {
                        runSpecificStage = 6;
                } else if ((arg == "-nBB") || (arg == "--noBranchAndBound"))  {
                        branchAndBound = false;
                } else if ((arg == "-nMM") || (arg == "--noMarkovModel"))  {
                        MarkovFilter = false;
                } 
                else if ((arg == "-nMEM") || (arg == "--noEssaMEM"))  {
                        essaMEM = false;
                } 
                else if ((arg == "-k") || (arg == "--kmersize")) {
                        i++;
                        if (i < argc)
                                kmerSize = atoi(args[i]);
                } else if ((arg == "-t") || (arg == "--threads")) {
                        i++;
                        if (i < argc)
                                numThreads = atoi(args[i]);
                } else if ((arg == "-v") || (arg == "--visits")) {
                        i++;
                        if (i < argc)
                                bubbleDFSNodeLimit = atoi(args[i]);
                } else if ((arg == "-d") || (arg == "--depth")) {
                        i++;
                        if (i < argc)
                                readCorrDFSNodeLimit = atoi(args[i]);
                } else if ((arg == "-e") || (arg == "--essa")) {
                        i++;
                        if (i < argc)
                                essaMEMSparsenessFactor = atoi(args[i]);
                } else if ((arg == "-c") || (arg == "--cutoff")) {
                        i++;
                        if (i < argc)
                                covCutoff = atoi(args[i]);
                } 
                /*else if ((arg == "-s") || (arg == "--singlestranded")) {
                        doubleStranded = false;
                }*/ else if ((arg == "-p") || (arg == "--pathtotmp")) {
                        i++;
                        if (i < argc)
                                pathtotemp = args[i];
                } else if ((arg == "-o") || (arg == "--output")) {
                        i++;
                        if (i < argc)
                                outputFilename = args[i];
                } else {        // it must be an input file
                        inputFilename = args[i];
                        libraries.push_back(pair<string, string>(inputFilename, outputFilename));
                        outputFilename.clear();
                }
        }

        // perform sanity check on input parameters
        if (numThreads == 0)
                numThreads = 1;

        if (numThreads > std::thread::hardware_concurrency()) {
                cerr << "WARNING: number of threads is bigger than the available number of cores" << endl;
        }

        if (kmerSize <= KMERBYTEREDUCTION *4) {
                cerr << "The kmer size must be at least " << 4*KMERBYTEREDUCTION + 1 << endl;
                throw ("Invalid argument");
        }

        if (kmerSize % 2 == 0) {
                cerr << "The kmer size must be odd" << endl;
                throw ("Invalid argument");
        }

        if (kmerSize > MAXKMERLENGTH) {
                size_t maxKmerLength = MAXKMERLENGTH;
                if (maxKmerLength % 2 == 0)
                        maxKmerLength--;
                cerr << "The kmer size can be at most " << maxKmerLength << endl;
                cerr << "Recompile Brownie with a higher MAXKMERLENGTH if a higher kmer size is desired" << endl;
                throw ("Invalid argument");
        }

        if (!pathtotemp.empty()) {
                if ((pathtotemp.back() != '/') && (pathtotemp.back() != '\\'))
                        pathtotemp.push_back('/');
        }

        // final check: see if we can write to the temporary directory
        ofstream ofs(pathtotemp + "log.txt");
        if (!ofs.good()) {
                cerr << "brownie: cannot write to directory: " << pathtotemp << "\n";
                cerr << "Please make sure the path exists" << endl;
                exit(EXIT_FAILURE);
        }

        for (int i = 0; i < argc; i++)
                ofs << args[i] << " ";
        ofs << endl;

        ofs.close();

        // add the libaries to the library container
        for (auto it : libraries) {
                ReadLibrary lib = ReadLibrary(it.first, it.second, getTempDirectory());
                libCont.insert(lib);
        }

        if (libCont.getSize() == 0) {
                cerr << "brownie: missing input read file(s)\n";
                cerr << "Try 'brownie assemble --help' for more information" << endl;
                exit(EXIT_FAILURE);
        }

        // try to read the metadata for each library
        libCont.readMetadata(getTempDirectory());
        string arg(args[1]);
        if (arg =="index"){
                runSpecificStage = 2;
                abundanceMin = 1;
        }
        if (arg =="align")
                runSpecificStage = 3;
        
}

void Settings::parseCommandLineArgVisualize(int argc, char** args,
                                            LibraryContainer& libCont)
{
        // parse all input arguments
        string inputFilename, outputFilename;
        vector<pair<string, string> > libraries;
        for (int i = 2; i < argc; i++) {
                string arg(args[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsageVisualize();
                        exit(EXIT_SUCCESS);
                }
        }
}

void Settings::parseCommandLineArgCompare(int argc, char** args,
                                          LibraryContainer& libCont)
{
        // set default values
        runSpecificStage = 4;

        // parse all input arguments
        string inputFilename, outputFilename;
        vector<pair<string, string> > libraries;
        for (int i = 2; i < argc; i++) {
                string arg(args[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsageCompare();
                        exit(EXIT_SUCCESS);
                } else if ((arg == "-s3") || (arg == "--stage3")) {
                        runSpecificStage = 3;
                } else if ((arg == "-s4") || (arg == "--stage4")) {
                        runSpecificStage = 4;
                } else if ((arg == "-k") || (arg == "--kmersize")) {
                        i++;
                        if (i < argc)
                                kmerSize = atoi(args[i]);
                } else if ((arg == "-p") || (arg == "--pathtotmp")) {
                        i++;
                        if (i < argc)
                                pathtotemp = args[i];
                } else {        // it must be a reference filename
                        referenceFilename = string(args[i]);
                }
        }

        if (kmerSize <= KMERBYTEREDUCTION *4) {
                cerr << "The kmer size must be at least " << 4*KMERBYTEREDUCTION + 1 << endl;
                throw ("Invalid argument");
        }

        if (kmerSize % 2 == 0) {
                cerr << "The kmer size must be odd" << endl;
                throw ("Invalid argument");
        }

        if (kmerSize > MAXKMERLENGTH) {
                size_t maxKmerLength = MAXKMERLENGTH;
                if (maxKmerLength % 2 == 0)
                        maxKmerLength--;
                cerr << "The kmer size can be at most " << maxKmerLength << endl;
                cerr << "Recompile Brownie with a higher MAXKMERLENGTH if a higher kmer size is desired" << endl;
                throw ("Invalid argument");
        }

        if (!pathtotemp.empty()) {
                if ((pathtotemp.back() != '/') && (pathtotemp.back() != '\\'))
                        pathtotemp.push_back('/');
        }

        // final check: see if we can write to the temporary directory
        ofstream ofs(pathtotemp + "log.txt");
        if (!ofs.good()) {
                cerr << "brownie: cannot write to directory: " << pathtotemp << "\n";
                cerr << "Please make sure the path exists" << endl;
                exit(EXIT_FAILURE);
        }

        for (int i = 0; i < argc; i++)
                ofs << args[i] << " ";
        ofs << endl;

        ofs.close();
}
