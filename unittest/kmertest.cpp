#include <gtest/gtest.h>
#include <cstdio>
#include "tkmer.h"
#include "tstring.h"

using namespace std;

const int maxKmer = 53;
const int numBytes = (maxKmer + 3) / 4;

typedef TKmer<numBytes> TestKmer;
typedef TKmer<numBytes -2 > RedKmer;
string source("AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA");

TEST(kmer, constrCppTest)
{
        // c++ string constructor
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);
                TestKmer kmer(kmerString);

                EXPECT_EQ(kmer.str() == kmerString, true);
        }

        // c++ string constructor with offset
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                for (size_t offset = 1; offset < source.size() - i; offset++) {
                        string kmerString = source.substr(0, i+offset);
                        TestKmer kmer(kmerString, offset);

                        EXPECT_EQ(kmer.str() == kmerString.substr(offset), true);
                }
        }
}

TEST(kmer, constrCTest)
{
        // c string constructor
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);
                TestKmer kmer(kmerString.c_str());

                EXPECT_EQ(kmer.str() == kmerString, true);
        }
}

TEST(kmer, constrTStringTest)
{
        // tight string contructor
        TString tStr(source);

        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                for (size_t offset = 0; offset < source.size() - i; offset++) {
                        // construct a tstring and check if construction succeeded
                        TestKmer kmer1(tStr, offset);
                        TestKmer kmer2(source, offset);

                        EXPECT_EQ(kmer1, kmer2);
                }
        }
}

TEST(kmer, KmerToReducedAndBack)
{
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                // construct a tstring and check if construction succeeded
                TestKmer kmer(source);

                KmerLSB lsb;
                RedKmer rkmer(kmer, lsb);
                TestKmer kmer2(rkmer, lsb);

                EXPECT_EQ(kmer, kmer2);
        }
}

TEST(kmer, streamTest)
{
        // construct a kmer and check if construction succeeded
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);
                TestKmer kmer(kmerString);

                stringstream ss;
                ss << kmer;

                EXPECT_EQ(ss.str() == kmerString, true);
        }
}

TEST(kmer, pushNucleotideLeft)
{
        // c++ string constructor
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);
                TestKmer kmer(kmerString);
                kmer.setFlag1(true);

                kmer.pushNucleotideLeft('A');
                kmer.pushNucleotideLeft('C');
                kmer.pushNucleotideLeft('G');
                kmer.pushNucleotideLeft('T');

                string target = "TGCA" + kmerString;

                EXPECT_EQ(kmer.str() == target.substr(0, i), true);
                EXPECT_EQ(kmer.getFlag1(), true);
                EXPECT_EQ(kmer.getFlag2(), false);
        }
}

TEST(kmer, pushNucleotideRight)
{
        // c++ string constructor
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);
                TestKmer kmer(kmerString);
                kmer.setFlag2(true);

                kmer.pushNucleotideRight('A');
                kmer.pushNucleotideRight('C');
                kmer.pushNucleotideRight('G');
                kmer.pushNucleotideRight('T');

                string target = kmerString + "ACGT";

                EXPECT_EQ(kmer.str() == target.substr(4, i), true);
                EXPECT_EQ(kmer.getFlag2(), true);
                EXPECT_EQ(kmer.getFlag1(), false);
        }
}

TEST(kmer, peekNucleotideLeft)
{
        // c++ string constructor
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);
                TestKmer kmer(kmerString);

                EXPECT_EQ(kmer.peekNucleotideLeft() == kmerString[0], true);
        }
}

TEST(kmer, peekNucleotideRight)
{
        // c++ string constructor
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);
                TestKmer kmer(kmerString);

                EXPECT_EQ(kmer.peekNucleotideRight() == kmerString[i-1], true);
        }
}

TEST(kmer, reverseTest)
{
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);

                TestKmer kmer(kmerString);
                kmer.setFlag2(true);
                kmer.reverse();

                reverse(kmerString.begin(), kmerString.end());

                EXPECT_EQ(kmer.str() == kmerString, true);
                EXPECT_EQ(kmer.getFlag2(), true);
                EXPECT_EQ(kmer.getFlag1(), false);
        }
}

TEST(kmer, complementTest)
{
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);

                TestKmer kmer(kmerString);
                kmer.setFlag1(true);
                kmer.complement();

                Nucleotide::complement(kmerString);
                EXPECT_EQ(kmer.str() == kmerString, true);
                EXPECT_EQ(kmer.getFlag1(), true);
                EXPECT_EQ(kmer.getFlag2(), false);
        }
}

TEST(kmer, revComplTest)
{
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);

                TestKmer kmer(kmerString);
                kmer.setFlag1(true);
                kmer.reverseComplement();

                Nucleotide::revCompl(kmerString);
                EXPECT_EQ(kmer.str() == kmerString, true);
                EXPECT_EQ(kmer.getFlag1(), true);
                EXPECT_EQ(kmer.getFlag2(), false);
        }
}

TEST(kmer, compareTest)
{
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);
                string kmerString = source.substr(0, i);

                string complement = kmerString;
                Nucleotide::revCompl(kmerString);

                TestKmer kmer(kmerString);
                kmer.setFlag1(true);
                TestKmer kmerRC = kmer;
                kmerRC.reverseComplement();
                kmer.setFlag2(true);
                TestKmer kmer2 = kmer;
                kmer.setFlag1(false);

                EXPECT_EQ(kmer < kmerRC, kmerString < complement);
                EXPECT_EQ(kmer == kmerRC, kmerString == complement);
                EXPECT_EQ(kmer == kmer2, true);
                EXPECT_EQ(kmer > kmerRC, kmerString > complement);
        }
}

TEST(kmer, hashTest)
{
        for (int i = 1; i <= maxKmer; i += 2) {
                TestKmer::setWordSize(i);

                TestKmer A(source);
                TestKmer B(source);
                B.setFlag1(true);
               // A.setFlag2(true);

                EXPECT_EQ(A.getHash() == B.getHash(), true);
        }
}
