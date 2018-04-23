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

#include "brownie.h"

#include "kmeroverlap.h"
#include "kmeroverlaptable.h"
#include "kmertable.h"
#include "readcorrection.h"


#include <cmath>

using namespace std;

Brownie::Brownie(int argc, char** args)
{
        settings.parseCommandLineArguments(argc, args, libraries);
        Kmer::setWordSize(settings.getK());
        RKmer::setWordSize(settings.getK() - KMERBYTEREDUCTION * 4);
}

void Brownie::stageOne()
{
        // ============================================================
        // STAGE 1 : PARSE THE READS
        // ============================================================

        cout << "\nEntering stage 1" << endl;
        cout << "================" << endl;

        KmerTable *readParser = new KmerTable(settings);
        cout << "Generating kmers with k = " << Kmer::getK()
             << " from input files..." << endl;
        Util::startChrono();
        readParser->parseInputFiles(libraries);

        size_t kmerGOne = readParser->getNumKmersCovGTOne();
        size_t allKmers = readParser->getNumKmers() ;
        cout << "Parsed input files (" << Util::stopChronoStr() << ")" << endl;
        cout << "Total number of unique kmers in table: "
             << allKmers << " (" << kmerGOne << " with coverage > 1)" << endl;

#ifdef DEBUG
        readParser->validateStage1();
#endif

        // write kmers file containing all kmers with cov > 1
        cout << "Writing kmer file...";
        cout.flush();
        Util::startChrono();
        if (settings.getAbundanceMinValue() ==1)
                readParser->writeAllKmers(getKmerFilename());
        else
        readParser->writeKmersWithCovGTOne(getKmerFilename());
        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        delete readParser;

        // write metadata for all libraries
        libraries.writeMetadata(settings.getTempDirectory());

        cout << "Stage 1 finished.\n" << endl;
}

void Brownie::stageTwo()
{
        // ============================================================
        // STAGE 2 : KMER OVERLAP TABLE
        // ============================================================

        cout << "\nEntering stage 2" << endl;
        cout << "================" << endl;

        // create a kmer table from the reads
        KmerOverlapTable overlapTable(settings);
        Util::startChrono();
        cout << "Building kmer overlap table...";
        overlapTable.loadKmersFromDisc(getKmerFilename());
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << "Number of kmers loaded: " << overlapTable.size() << endl;

        // find the overlap between kmers
        Util::startChrono();
        cout << "Finding overlaps between kmers..." << endl;
        overlapTable.parseInputFiles(libraries);
        cout << "Done building overlap table ("
             << Util::stopChronoStr() << ")" << endl;
        cout << "Overlap table contains " << overlapTable.size()
             << " nodes" << endl;

#ifdef DEBUG
        overlapTable.validateStage2();
#endif

        // extract nodes and arcs from the kmer table
        overlapTable.extractNodes(getNodeFilename(2),
                                  getArcFilename(2),
                                  getMetaDataFilename(2));

        overlapTable.clear();   // clear memory !
        cout << "Stage 2 finished.\n" << endl;
}

void Brownie::stageThree()
{
      
        cout << "\nEntering stage 3" << endl;
        cout << "================" << endl;


        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... "; cout.flush();
        graph.loadGraph(getNodeFilename(2),
                        getArcFilename(2),
                        getMetaDataFilename(2));

        cout << "Graph contains " << graph.getNumNodes() << " nodes and "
             << graph.getNumArcs() << " arcs" << endl;
        ReadCorrectionHandler rcHandler(graph, settings, settings.getMarkovFilter());
        rcHandler.doErrorCorrection(libraries);
        rcHandler.doReadRefinement(libraries);
        //MarkovChainHandler mch(graph,settings);
        //mch.loadFromFile();
        //mch.verifyChains(settings.getTempDirectory()+ "/perfectReads.ncf");
        cout << "Read alignment completed in " << Util::stopChronoStr() << endl;
        cout << "Stage 3 finished\n" << endl;
        graph.clear();
      
}



void Brownie::assembleModule()
{       
        if (stageOneNecessary())
                stageOne();
        else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 1...\n";
        if (stageTwoNecessary())
                stageTwo();
        else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 2...\n";
        if (stageThreeNecessary())
                stageThree();
        else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 3...\n";

}

void Brownie::run()
{
        switch (settings.getCommand()) {
                case Command::assemble:
                        assembleModule();
                        break;
                case Command::none:
                        cerr << "brownie: no command specified\n";
                        cerr << "Try 'brownie --help' for more information" << endl;
                        exit(EXIT_FAILURE);
                        break;
        }
}

int main(int argc, char** args)
{
        try {
                Brownie brownie(argc, args);

                cout << "Welcome to Brownie v." << BROWNIE_MAJOR_VERSION << "."
                     << BROWNIE_MINOR_VERSION << "." << BROWNIE_PATCH_LEVEL;
#ifdef DEBUG
                cout << " (debug mode)" << endl;
#else
                cout << " (release mode)" << endl;
#endif
                cout << "Today is " << Util::getDateTime() << endl;
                brownie.run();
        } catch (exception &e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl << endl;
        return EXIT_SUCCESS;
}
