/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
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

#include "global.h"
#include "kmeroverlaptable.h"
#include "kmeroverlap.h"
#include "settings.h"
#include "tstring.h"
#include "library.h"

#include "readfile/fastafile.h"
#include "readfile/fastqfile.h"
#include "readfile/rawfile.h"
#include "readfile/sequencefile.h"
#include "readfile/samfile.h"

#include <deque>
#include <set>
#include <map>
#include <fstream>
#include <iostream>

using namespace std;

// ============================================================================
// KMER OVERLAP TABLE
// ============================================================================

KmerOverlapRef KmerOverlapTable::insert(const Kmer &kmer)
{
        // chose a representative kmer
        Kmer representative = settings.isDoubleStranded() ?
                kmer.getRepresentative() : kmer;

        bool reverse = (kmer != representative);

        KmerOverlapPair val(representative, KmerOverlap());
        pair<KmerOverlapIt, bool> insResult = table.insert(val);
        KmerOverlapRef result(insResult.first, reverse);

        return result;
}

KmerOverlapRef KmerOverlapTable::find(const Kmer &kmer) const
{
        // chose a representative kmer
        Kmer representative = settings.isDoubleStranded() ?
                kmer.getRepresentative() : kmer;

        bool reverse = (kmer != representative);
        return KmerOverlapRef(table.find(representative), reverse);
}

bool KmerOverlapTable::getLeftUniqueKmer(const KmerOverlapRef& rKmerRef,
                                  KmerOverlapRef& lKmerRef) const
{
        // initialise the right kmer reference to point to nothing
        lKmerRef = KmerOverlapRef(table.end(), false);

        char nucleotide;
        if (!rKmerRef.hasLeftUniqueOverlap(nucleotide))
                return false;

        Kmer lKmer = rKmerRef.getKmer();
        lKmer.pushNucleotideLeft(nucleotide);
        lKmerRef = find(lKmer);

        assert(lKmerRef.first != table.end());

        if (!lKmerRef.hasRightUniqueOverlap(nucleotide))
                return false;

        return true;
}

bool KmerOverlapTable::getRightUniqueKmer(const KmerOverlapRef& lKmerRef,
                                   KmerOverlapRef& rKmerRef) const
{
        // initialise the right kmer reference to point to nothing
        rKmerRef = KmerOverlapRef(table.end(), false);

        char nucleotide;
        if (!lKmerRef.hasRightUniqueOverlap(nucleotide))
                return false;

        Kmer rKmer = lKmerRef.getKmer();
        rKmer.pushNucleotideRight(nucleotide);
        rKmerRef = find(rKmer);

        assert(rKmerRef.first != table.end());

        if (!rKmerRef.hasLeftUniqueOverlap(nucleotide))
                return false;

        return true;
}

void KmerOverlapTable::convertKmersToString(const deque<KmerOverlapRef> &kmerSeq,
                                            string &output)
{
        output = kmerSeq[0].getKmer().str();

        for (size_t i = 1; i < kmerSeq.size(); i++)
                output.push_back(kmerSeq[i].getKmer().peekNucleotideRight());
}

void KmerOverlapTable::loadKmersFromDisc(const std::string& filename)
{
        // load all the kmers from disc
        ifstream ifs(filename.c_str(), ios::in | ios::binary);
        size_t numKmers;
        ifs.read((char*)&numKmers, sizeof(size_t));

        table.resize(numKmers);

        for (size_t i = 0; i < numKmers; i++) {
                Kmer kmer(ifs);
                insert(kmer);
        }

        ifs.close();
}

void KmerOverlapTable::parseRead(string& read,
                                 vector<pair<Kmer, KmerOverlap> >& kmerBuffer)
{
        // get out early
        if (read.size() < Kmer::getK())
                return;

        // reserve space for kmers and flags
        vector<KmerOverlapRef> refs(read.size() + 1 - Kmer::getK());

        // find the kmers in the table
        for (KmerIt it(read); it.isValid(); it++)
                refs[it.getOffset()] = find(it.getKmer());
        // now mark the overlap implied by the read
        size_t lastIndex = 0;
        for (KmerIt it(read); it.isValid(); it++) {

                if (refs[it.getOffset()].first == table.end())
                        continue;
                if (it.hasRightOverlap())
                        if (refs[it.getOffset()+1].first != table.end())
                                refs[it.getOffset()].markRightOverlap(it.getRightOverlap());
                if (it.hasLeftOverlap())
                        if (refs[it.getOffset()-1].first != table.end())
                                refs[it.getOffset()].markLeftOverlap(it.getLeftOverlap());
                lastIndex = it.getOffset();
        }

        // get the first index of the read that should be kept
        size_t firstIndex = refs.size();
        for (KmerIt it(read); it.isValid(); it++) {
                if (refs[it.getOffset()].first == table.end())
                        continue;
                firstIndex = it.getOffset();
                break;
        }

}


/*
void KmerOverlapTable::parseRead(string& read,
                                 vector<pair<Kmer, KmerOverlap> >& kmerBuffer)
{


        // get out early
        if (read.size() < Kmer::getK())
                return;

        // reserve space for kmers and flags
        vector<KmerOverlapRef> refs(read.size() + 1 - Kmer::getK());

        // find the kmers in the table
        for (KmerIt it(read); it.isValid(); it++)
                refs[it.getOffset()] = find(it.getKmer());
        // now mark the overlap implied by the read
        size_t lastIndex = 0;
        for (KmerIt it(read); it.isValid(); it++) {

                if (refs[it.getOffset()].first == table.end())
                        continue;
                lastIndex = it.getOffset();
        }

        // get the first index of the read that should be kept
        size_t firstIndex = refs.size();
        for (KmerIt it(read); it.isValid(); it++) {
                if (refs[it.getOffset()].first == table.end())
                        continue;
                firstIndex = it.getOffset();
                break;
        }
        for (KmerIt it(read); it.isValid(); it++) {

                if (it.getOffset() < firstIndex)
                        continue;
                if (it.getOffset() > lastIndex)
                        continue;
                if (refs[it.getOffset()].first == table.end())
                        KmerOverlapRef ref = insert(it.getKmer());
        }

        for (KmerIt it(read); it.isValid(); it++)
                refs[it.getOffset()] = find(it.getKmer());

        for (KmerIt it(read); it.isValid(); it++) {
                if (refs[it.getOffset()].first == table.end())
                        continue;
                if (it.hasRightOverlap())
                        if (refs[it.getOffset()+1].first != table.end())
                                refs[it.getOffset()].markRightOverlap(it.getRightOverlap());
                if (it.hasLeftOverlap())
                        if (refs[it.getOffset()-1].first != table.end())
                                refs[it.getOffset()].markLeftOverlap(it.getLeftOverlap());
        }

}*/

void KmerOverlapTable::parseReads(size_t thisThread,
                                  vector<string>& readBuffer,
                                  vector<pair<Kmer, KmerOverlap> >& kmerBuffer)
{
        for (size_t i = 0; i < readBuffer.size(); i++)
                parseRead(readBuffer[i], kmerBuffer);
}

void KmerOverlapTable::workerThread(size_t thisThread, LibraryContainer* inputs)
{
        // aux variables
        vector<pair<Kmer, KmerOverlap> > kmerBuffer;
        // local storage of reads
        vector<string> myReadBuf;
        size_t blockID, recordOffset;
        while (inputs->getReadChunk(myReadBuf, blockID, recordOffset))
                parseReads(thisThread, myReadBuf, kmerBuffer);
}

void KmerOverlapTable::parseInputFiles(LibraryContainer &inputs)
{
        const unsigned int& numThreads =  settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        inputs.startIOThreads(settings.getThreadWorkSize(),
                              settings.getThreadWorkSize() * numThreads);

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&KmerOverlapTable::workerThread, this, i, &inputs);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        inputs.joinIOThreads();
}

void KmerOverlapTable::extractNodes(const string& nodeFilename,
                                    const string& arcFilename,
                                    const string& metaDataFilename)
{
        ofstream nodeFile(nodeFilename.c_str());
        ofstream arcFile(arcFilename.c_str());

        string descriptor;
        size_t numNodes = 0, numArcs = 0;
        size_t progressID = 1;
        for (KmerOverlapIt it = table.begin(); it != table.end(); it++) {

                KmerOverlapRef currKmer(it, false), nextKmer;

                if (progressID++ % OUTPUT_FREQUENCY == 0)
                        (cout << "Extracting node " << numNodes
                             << " from graph.\r").flush();

                // check if the node has been processed before
                if (currKmer.isProcessed()) continue;

                deque<KmerOverlapRef> kmerSeq;
                kmerSeq.push_back(currKmer);

                // extend node to the right
                while (getRightUniqueKmer(currKmer, nextKmer)) {

                        // check for a loop
                        if (nextKmer == kmerSeq.front())
                                break;
                        // check for a hairpin
                        if (nextKmer.first == kmerSeq.back().first)
                                break;

                        kmerSeq.push_back(nextKmer);
                        currKmer = nextKmer;
                }

                // extend node to the left
                currKmer = KmerOverlapRef(it, false);

                while (getLeftUniqueKmer(currKmer, nextKmer)) {

                        // check for a loop
                        if (nextKmer == kmerSeq.back())
                                break;
                        // check for a hairpin
                        if (nextKmer.first == kmerSeq.front().first)
                                break;

                        kmerSeq.push_front(nextKmer);
                        currKmer = nextKmer;
                }

                // mark all kmers as processed
                for (size_t i = 0; i < kmerSeq.size(); i++)
                        kmerSeq[i].setProcessed(true);

                progressID += kmerSeq.size();

                convertKmersToString(kmerSeq, descriptor);

                nodeFile << "NODE" << "\t" << numNodes++ << "\t"
                         << kmerSeq.size() << "\t" << "0" << "\t" << "0"
                         << "\n" << descriptor << "\n";

                arcFile << numNodes << "\t" << (int)kmerSeq[0].getLeftOverlap()
                        << "\t" << (int)kmerSeq.back().getRightOverlap();

                for (NucleotideID j = 0; j < 4; j++) {
                        char n = Nucleotide::nucleotideToChar(j);
                        if (kmerSeq[0].hasLeftOverlap(n))
                                arcFile << "\t" << 0;
                        if (kmerSeq.back().hasRightOverlap(n))
                                arcFile << "\t" << 0;
                }

                arcFile << "\n";

                numArcs += kmerSeq[0].getNumLeftOverlap() +
                           kmerSeq.back().getNumRightOverlap();
        }

        nodeFile.close();
        arcFile.close();

        ofstream mdFile(metaDataFilename.c_str());
        mdFile << numNodes << "\t" << numArcs << endl;
        mdFile.close();

        cout << "Extracted " << numNodes << " nodes and "
             << numArcs << " arcs." << endl;
}

#ifdef DEBUG
void KmerOverlapTable::validateStage2()
{
        FastAFile ass(false);
        ass.open("genome.fasta");

        size_t numKmersReal = 0, numKmersFound = 0, numOLFound = 0, numOLReal = 0;

        string read;
        while (ass.getNextRead (read)) {
                for (KmerIt it(read); it.isValid(); it++) {
                        Kmer kmer = it.getKmer();
                        numKmersReal++;

                        KmerOverlapRef result = find(kmer);

                        // not found, continue
                        if (result.first == table.end())
                                continue;

                        // we did find the kmer
                        numKmersFound++;

                        if (it.hasRightOverlap()) {
                                numOLReal++;
                                if (result.hasRightOverlap(it.getRightOverlap()))
                                        numOLFound++;
                        }
                }
        }

        ass.close();

        // validity check for kmer table
        for (auto& it : table) {
                Kmer kmer = it.first;
                KmerOverlap ol = it.second;

                if (ol.hasLeftOverlap('A')) {
                        Kmer left = kmer;
                        left.pushNucleotideLeft('A');
                        const auto& it2 = find(left);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has left OL with A but " << left << " not found in table" << endl;
                        if (!it2.hasRightOverlap(kmer.peekNucleotideRight()))
                                cerr << "Kmer " << kmer << " has left OL with A but " << left << " has no right OL with " << kmer.peekNucleotideRight() << endl;
                }

                if (ol.hasLeftOverlap('C')) {
                        Kmer left = kmer;
                        left.pushNucleotideLeft('C');
                        const auto& it2 = find(left);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has left OL with C but " << left << " not found in table" << endl;
                        if (!it2.hasRightOverlap(kmer.peekNucleotideRight()))
                                cerr << "Kmer " << kmer << " has left OL with C but " << left << " has no right OL with " << kmer.peekNucleotideRight() << endl;
                }

                if (ol.hasLeftOverlap('G')) {
                        Kmer left = kmer;
                        left.pushNucleotideLeft('G');
                        const auto& it2 = find(left);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has left OL with G but " << left << " not found in table" << endl;
                        if (!it2.hasRightOverlap(kmer.peekNucleotideRight()))
                                cerr << "Kmer " << kmer << " has left OL with G but " << left << " has no right OL with " << kmer.peekNucleotideRight() << endl;
                }

                if (ol.hasLeftOverlap('T')) {
                        Kmer left = kmer;
                        left.pushNucleotideLeft('T');
                        const auto& it2 = find(left);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has left OL with T but " << left << " not found in table" << endl;
                        if (!it2.hasRightOverlap(kmer.peekNucleotideRight()))
                                cerr << "Kmer " << kmer << " has left OL with T but " << left << " has no right OL with " << kmer.peekNucleotideRight() << endl;
                }

                if (ol.hasRightOverlap('A')) {
                        Kmer right = kmer;
                        right.pushNucleotideRight('A');
                        const auto& it2 = find(right);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has right OL with A but " << right << " not found in table" << endl;
                        if (!it2.hasLeftOverlap(kmer.peekNucleotideLeft()))
                                cerr << "Kmer " << kmer << " has right OL with A but " << right << " has no left OL with " << kmer.peekNucleotideLeft() << endl;
                }

                if (ol.hasRightOverlap('C')) {
                        Kmer right = kmer;
                        right.pushNucleotideRight('C');
                        const auto& it2 = find(right);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has right OL with C but " << right << " not found in table" << endl;
                        if (!it2.hasLeftOverlap(kmer.peekNucleotideLeft()))
                                cerr << "Kmer " << kmer << " has right OL with C but " << right << " has no left OL with " << kmer.peekNucleotideLeft() << endl;
                }

                if (ol.hasRightOverlap('G')) {
                        Kmer right = kmer;
                        right.pushNucleotideRight('G');
                        const auto& it2 = find(right);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has right OL with G but " << right << " not found in table" << endl;
                        if (!it2.hasLeftOverlap(kmer.peekNucleotideLeft()))
                                cerr << "Kmer " << kmer << " has right OL with G but " << right << " has no left OL with " << kmer.peekNucleotideLeft() << endl;
                }

                if (ol.hasRightOverlap('T')) {
                        Kmer right = kmer;
                        right.pushNucleotideRight('T');
                        const auto& it2 = find(right);
                        if (it2.first == table.end())
                                cerr << "Kmer " << kmer << " has right OL with T but " << right << " not found in table" << endl;
                        if (!it2.hasLeftOverlap(kmer.peekNucleotideLeft()))
                                cerr << "Kmer " << kmer << " has right OL with T but " << right << " has no left OL with " << kmer.peekNucleotideLeft() << endl;
                }
        }

        cout << "Validation report: " << endl;
        cout << "\tExist in table: " << numKmersFound << "/" << numKmersReal << " (" << 100.00*(double)numKmersFound/(double)numKmersReal << "%)" << endl;
        cout << "\tHave correct overlap: " << numOLFound << "/" << numOLReal << " (" << 100.00*(double)numOLFound/(double)numOLReal << "%)" << endl;
}
#endif
