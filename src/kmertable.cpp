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

#include "kmertable.h"
#include "tkmer.h"
#include "tstring.h"
#include "library.h"
#include "settings.h"
#include "readfile/sequencefile.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <random>

#define OUTPUT_FREQUENCY 32768
#define NUMTABLES (KmerLSB(1) << 8*KMERBYTEREDUCTION)

using namespace std;

// ============================================================================
// MIXING FUNCTION CLASS (PULBIC)
// ============================================================================

MixingLSB::MixingLSB()
{
        mixingF = vector<KmerLSB>(1 << KMERBYTEREDUCTION*8);
        imixingF = vector<KmerLSB>(1 << KMERBYTEREDUCTION*8);

        std::random_device rd;
        std::mt19937 gen(rd());

        // generate unshuffled input
        for (size_t i = 0; i < (1 << KMERBYTEREDUCTION*8); i++)
                mixingF[i] = i;

        // Fisher-Yates algorithm to shuffle input
        for (size_t i = 0; i < ((1 << KMERBYTEREDUCTION*8)-1); i++) {
                std::uniform_int_distribution<> dis(i, (1 << (KMERBYTEREDUCTION*8))-1);
                size_t j = dis(gen);
                swap(mixingF[i], mixingF[j]);
        }

        // compute the inverse table
        for (size_t i = 0; i < (1 << KMERBYTEREDUCTION*8); i++)
                imixingF[mixingF[i]] = i;
}

KmerLSB MixingLSB::mix(KmerLSB input) const
{
        return mixingF[input];
}

KmerLSB MixingLSB::invmix(KmerLSB input) const
{
        return imixingF[input];
}

// ============================================================================
// READ PARSER (PRIVATE)
// ============================================================================

size_t KmerTable::getThreadIDForKmer(const Kmer& kmer) const
{
        // split kmer = [lsb][reducedKmer]
        KmerLSB lsb;
        RKmer reducedKmer(kmer, lsb);

        // mix the lsb to obtain a more uniform distribution of the workload
        lsb = mixFunction.mix(lsb);

        // get the threadID based on the lsb
        size_t threadID = ((size_t)lsb * settings.getNumThreads()) / NUMTABLES;
        if ((((threadID+1)*NUMTABLES) / settings.getNumThreads()) <= (size_t)lsb)
                threadID++;

        return threadID;
}

size_t KmerTable::parseRead(string &read, vector<Kmer> *kmerBuffer)
{
        // read too short ?
        if (read.size() < Kmer::getK())
                return 0;

        // transform to uppercase
        transform(read.begin(), read.end(), read.begin(), ::toupper);

        size_t numKmers = 0;
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();

                // choose a representative kmer
                Kmer representative = settings.isDoubleStranded() ?
                        kmer.getRepresentative() : kmer;

                size_t threadID = getThreadIDForKmer(representative);

                kmerBuffer[threadID].push_back(representative);
                numKmers++;
        }

        return numKmers;
}

void KmerTable::parseReads(size_t thisThread,
                            vector<string>& readBuffer,
                            vector<Kmer>* tempKmerBuffer,
                            vector<Kmer>& myKmerBuf)
{
        for (size_t i = 0; i < readBuffer.size(); i++)
                parseRead(readBuffer[i], tempKmerBuffer);

        // push temporary k-mers onto the correct stacks
        for (size_t i = 0; i < settings.getNumThreads(); i++) {
                sharedKmerBufMutex[i].lock();
                if (i != thisThread) {  // copy temp kmers onto other thread's workstack
                        sharedKmerBuf[i].insert(sharedKmerBuf[i].end(),
                                                tempKmerBuffer[i].begin(),
                                                tempKmerBuffer[i].end());
                } else {                // copy other thread's workstack onto local stack
                        myKmerBuf.insert(myKmerBuf.end(),
                                         sharedKmerBuf[i].begin(),
                                         sharedKmerBuf[i].end());
                        sharedKmerBuf[i].clear();
                }
                sharedKmerBufMutex[i].unlock();
                if (i == thisThread)    // copy temp kmers onto local workstack
                        myKmerBuf.insert(myKmerBuf.end(),
                                         tempKmerBuffer[i].begin(),
                                         tempKmerBuffer[i].end());
                tempKmerBuffer[i].clear();
        }
}

void KmerTable::storeKmersInTable(size_t thisThread,
                                   const vector<Kmer>& myKmerBuf)
{
        size_t firstTable = (thisThread * NUMTABLES) / settings.getNumThreads();

        // store all kmers in the hash table
        for (size_t i = 0; i < myKmerBuf.size(); i++) {
                KmerLSB lsb;
                RKmer reducedKmer(myKmerBuf[i], lsb);
                lsb = mixFunction.mix(lsb);
                auto insResult = tableThread[thisThread][lsb-firstTable].insert(reducedKmer);

                // if the kmer was inserted for the first time, do nothing
                if (insResult.second)
                        continue;

                // else, toggle a bit to indicate multiple occurences
                RKmer &foundKmer = const_cast<RKmer&>(*insResult.first);
                foundKmer.setFlag1(true);
        }
}

void KmerTable::workerThread(size_t thisThread, LibraryContainer* inputs)
{
        const unsigned int& numThreads = settings.getNumThreads();

        // hash tables
        size_t firstTable = (thisThread * NUMTABLES) / numThreads;
        size_t lastTable = ((thisThread + 1) * NUMTABLES) / numThreads;
        size_t numTables = lastTable - firstTable;
        tableThread[thisThread] = new RKmerHashTable[numTables];

        for (size_t i = firstTable; i < lastTable; i++)
                tables[i] = &tableThread[thisThread][i-firstTable];

        // local storage of reads
        vector<string> myReadBuf;

        // temporary buffers
        vector<Kmer> myKmerBuf;
        vector<Kmer> *tempKmerBuf = new vector<Kmer>[numThreads];

        size_t numKmersInTable = 0;

        while (true) {
                size_t thisNumReads = 0;

                // get work from other threads
                sharedKmerBufMutex[thisThread].lock();
                myKmerBuf.insert(myKmerBuf.end(),
                                 sharedKmerBuf[thisThread].begin(),
                                 sharedKmerBuf[thisThread].end());
                sharedKmerBuf[thisThread].clear();
                sharedKmerBufMutex[thisThread].unlock();

                // if there are no kmers to store, produce local k-mers
                if (myKmerBuf.size() == 0) {
                        // get a number of reads (mutex lock)
                        size_t blockID, recordOffset;
                        inputs->getReadChunk(myReadBuf, blockID, recordOffset);

                        // process these input reads (lock-free)
                        parseReads(thisThread, myReadBuf, tempKmerBuf, myKmerBuf);

                        thisNumReads = myReadBuf.size();
                        myReadBuf.clear();

                        if (thisNumReads > 0) {
                                unique_lock<mutex> lock(terminateMutex);
                                terminateCV.notify_all();
                                lock.unlock();
                        }
                }

                // actually store the kmers in a table
                storeKmersInTable(thisThread, myKmerBuf);

                numKmersInTable += myKmerBuf.size();
                size_t thisNumKmers = myKmerBuf.size();
                myKmerBuf.clear();

                // if we were able to do something this iteration, continue
                if ((thisNumKmers != 0) || (thisNumReads != 0))
                        continue;

                // at this point, the worker thread could do nothing
                unique_lock<mutex> lock(terminateMutex);
                numThreadReady++;       // ready and waiting
                if (numThreadReady == settings.getNumThreads()) {
                        terminateCV.notify_all();
                        break;
                }

                terminateCV.wait(lock);
                if ((sharedKmerBuf[thisThread].empty()) &&
                    (numThreadReady == settings.getNumThreads()))
                        break;

                numThreadReady--;
                lock.unlock();
        }

        delete [] tempKmerBuf;

        //cout << "Thread " << thisThread << " exiting ... "
        //     << "(num kmers = " << numKmersInTable << ")" << endl;
}

// ============================================================================
// READ PARSER (PUBLIC)
// ============================================================================

KmerTable::~KmerTable()
{
        if (tableThread != NULL) {
                for (size_t i = 0; i < settings.getNumThreads(); i++) {
                        if (tableThread[i] != NULL)
                                delete [] tableThread[i];
                        tableThread[i] = NULL;
                }
        }

        delete [] tableThread; tableThread = NULL;
        delete [] tables; tables = NULL;
}

void KmerTable::parseInputFiles(LibraryContainer &inputs)
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        sharedKmerBuf = vector<vector<Kmer> >(numThreads);
        sharedKmerBufMutex = vector<mutex>(numThreads);

        tableThread = new RKmerHashTable*[numThreads];
        tables = new RKmerHashTable*[NUMTABLES];

        numThreadReady = 0;

        inputs.startIOThreads(settings.getThreadWorkSize(),
                              settings.getThreadWorkSize() * settings.getNumThreads());

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&KmerTable::workerThread, this, i, &inputs);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        inputs.joinIOThreads();
}

void KmerTable::clear()
{
        if (tables == NULL)
                return;

        for (size_t i = 0; i < NUMTABLES; i++)
                tables[i]->clear();
}

size_t KmerTable::getNumKmers() const
{
        if (tables == NULL)
                return 0;

        size_t numKmers = 0;
        for (size_t i = 0; i < NUMTABLES; i++)
                numKmers += tables[i]->size();

        return numKmers;
}

size_t KmerTable::getNumKmersCovGTOne() const
{
        if (tables == NULL)
                return 0;

        size_t numKmers = 0;
        for (KmerLSB lsb = 0; lsb < NUMTABLES; lsb++) {
                KmerLSB lsbinv = mixFunction.invmix(lsb);
                for (auto it : *tables[lsb]) {
                        Kmer kmer(it, lsbinv);
                        if (kmer.getFlag1())
                                numKmers++;
                }
        }

        return numKmers;
}

void KmerTable::writeAllKmers(const string& filename)
{
        // first, write the number of kmers to the file
        size_t size = getNumKmers();
        ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
        ofs.write((char*)(&size), sizeof(size_t));

        // write all the kmers
        for (KmerLSB lsb = 0; lsb < NUMTABLES; lsb++) {
                KmerLSB lsbinv = mixFunction.invmix(lsb);
                google::sparse_hash_set<RKmer, RKmerHash> &table = *tables[lsb];
                for (auto it : table) {
                        Kmer kmer(it, lsbinv);
                        kmer.writeNoFlags(ofs);
                }
        }

        ofs.close();
}

void KmerTable::writeKmersWithCovGTOne(const string& filename)
{
        // first, write the number of kmers to the file
        size_t size = getNumKmersCovGTOne();
        ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
        ofs.write((char*)(&size), sizeof(size_t));

        // write all the kmers
        for (KmerLSB lsb = 0; lsb < NUMTABLES; lsb++) {
                KmerLSB lsbinv = mixFunction.invmix(lsb);
                google::sparse_hash_set<RKmer, RKmerHash> &table = *tables[lsb];
                for (auto it : table) {
                        if (!it.getFlag1()) continue;  // uniqueness flag
                        Kmer kmer(it, lsbinv);
                        kmer.writeNoFlags(ofs);
                }
        }

        ofs.close();
}

std::pair< bool, bool > KmerTable::find(const Kmer& kmer) const
{
        // chose a representative kmer
        Kmer representative = settings.isDoubleStranded() ?
                kmer.getRepresentative() : kmer;

        // create and store the reduced kmer
        KmerLSB lsb;
        RKmer reducedKmer(representative, lsb);
        lsb = mixFunction.mix(lsb);

        auto it = tables[lsb]->find(reducedKmer);
        if (it == tables[lsb]->end())
                return pair<bool, bool>(false, false);

        RKmer &foundKmer = const_cast<RKmer&>(*it);
        return pair<bool, bool>(true, foundKmer.getFlag1());
}

#ifdef DEBUG
#include "readfile/fastafile.h"

void KmerTable::validateStage1()
{
        FastAFile ass(false);
        ass.open("genome.fasta");

        size_t numKmers = 0, numFound = 0, numCovGTOne = 0;

        string read;
        while (ass.getNextRead (read)) {
                for (KmerIt it(read); it.isValid(); it++) {
                        Kmer kmer = it.getKmer();

                        pair<bool, bool> result = find(kmer);
                        if ( result.first ) {
                                numFound++;
                                if (result.second)
                                        numCovGTOne++;
                        }
                        numKmers++;
                }
        }

        ass.close();

        double fracFound = 100.0*(double)numFound/(double) numKmers;
        double fracCovGTOne = 100.0*(double)numCovGTOne/(double)numKmers;

        cout.precision (4);
        cout << "Validation report: " << endl;
        cout << "\tk-mers in table: " << numFound << "/" << numKmers
        << " (" << fracFound << "%)" << endl;
        cout << "\tk-mers in table cov GT one: " << numCovGTOne << "/"
        << numKmers << " (" << fracCovGTOne << "%)" << endl;
}

#endif
