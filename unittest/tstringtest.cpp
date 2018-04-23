#include <gtest/gtest.h>
#include <cstdio>
#include "tstring.h"

using namespace std;

TEST(TString, appendTest)
{
        string source1("ACGTACGTACGTGGATTCCCGA");
        string source2("CGGTAGGCTTAAAATTGCC");

        // construct a tstring and check if construction succeeded
        TString tStr1(source1);
        TString tStr2(source2);

        tStr1.append(tStr2);

        ostringstream iss;
        iss << tStr1;

        EXPECT_EQ(iss.str() == source1.append(source2), true);
}

TEST(TString, squareBracket)
{
        string source1("ACGTACGTACGTGGATTCCCGA");

        TString tStr1(source1);

        EXPECT_EQ(tStr1[0], 'A');
        EXPECT_EQ(tStr1[1], 'C');
        EXPECT_EQ(tStr1[2], 'G');
        EXPECT_EQ(tStr1[3], 'T');
        EXPECT_EQ(tStr1[4], 'A');
        EXPECT_EQ(tStr1[5], 'C');
        EXPECT_EQ(tStr1[6], 'G');
        EXPECT_EQ(tStr1[7], 'T');
        EXPECT_EQ(tStr1[8], 'A');
        EXPECT_EQ(tStr1[9], 'C');
        EXPECT_EQ(tStr1[10], 'G');
        EXPECT_EQ(tStr1[11], 'T');
        EXPECT_EQ(tStr1[12], 'G');
        EXPECT_EQ(tStr1[13], 'G');
        EXPECT_EQ(tStr1[14], 'A');
        EXPECT_EQ(tStr1[15], 'T');
        EXPECT_EQ(tStr1[16], 'T');
        EXPECT_EQ(tStr1[17], 'C');
        EXPECT_EQ(tStr1[18], 'C');
        EXPECT_EQ(tStr1[19], 'C');
        EXPECT_EQ(tStr1[20], 'G');
        EXPECT_EQ(tStr1[21], 'A');
}

TEST(TString, subString)
{
        string source1("ACGTACGTACGTGGATTCCCGA");

        TString tStr1(source1);

        EXPECT_EQ(tStr1.substr(0, 1), "A");
        EXPECT_EQ(tStr1.substr(10, 1), "G");
        EXPECT_EQ(tStr1.substr(10, 3), "GTG");
        EXPECT_EQ(tStr1.substr(10, 30), "GTGGATTCCCGA");
}

TEST(TString, toStringTest)
{
        string source1("ACGTACGTACGTGGATTCCCGA");

        TString tStr1(source1);

        EXPECT_EQ(tStr1.getSequence() == source1, true);
}

TEST(TString, peekNucleotideTest)
{
        string source1("GCGTACGTACGTGGATTCCCGT");

        TString tStr1(source1);

        EXPECT_EQ(tStr1.peekNucleotideLeft(), 'G');
        EXPECT_EQ(tStr1.peekNucleotideRight(), 'T');
}

TEST(TString, constructorTest)
{
        string source("ACGTACGTACGTGGATTCCCGA");

        // construct a tstring and check if construction succeeded
        TString tStr(source);

        ostringstream iss;
        iss << tStr;

        EXPECT_EQ(iss.str() == source, true);
}

TEST(TString, iteratorTest)
{
        Kmer::setWordSize(13);
        string source("ACGTACGTACGTGGATTCTTAGCCGTACGCCG");
        TString tStr(source);

        size_t cnt = 0;
        for (TString::iterator it = tStr.begin(); it != tStr.end(); it++)
                EXPECT_EQ(*it, Kmer(source, cnt++));

        EXPECT_EQ(cnt, source.size() - Kmer::getK() + 1);
}

TEST(TString, reverseIteratorTest)
{
        Kmer::setWordSize(13);
        string source("ACGTACGTACGTGGATTCTTAGCCGTACGCCG");
        TString tStr(source);

        int cnt = source.size() - Kmer::getK();
        for (TString::reverse_iterator it = tStr.rbegin(); it != tStr.rend(); it++)
                EXPECT_EQ(*it, Kmer(source, cnt--));

        EXPECT_EQ(cnt, -1);
}

TEST(TString, complementTest)
{
        string source("ACGTACGTACGTGGATTCTTAGCCGTACGCCGA");
        TString tstring(source);
        tstring.complement();

        EXPECT_EQ(Nucleotide::getComplement(source) == tstring.getSequence(), true);
}

TEST(TString, reverseTest)
{
        string source("ACGTACGTACGTGGATTCTTAGCCGTACGCCGA");
        TString tstring(source);
        tstring.reverse();

        EXPECT_EQ(Nucleotide::getReverse(source) == tstring.getSequence(), true);
}

TEST(TString, revComplTest)
{
        string source("ACGTACGTACGTGGATTCTTAGCCGTACGCCGA");
        TString tstring(source);
        tstring.reverseComplement();

        EXPECT_EQ(Nucleotide::getRevCompl(source) == tstring.getSequence(), true);
}
