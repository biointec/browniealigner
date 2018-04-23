#include <gtest/gtest.h>
#include <cstdio>
#include "kmeroverlap.h"
#include "tkmer.h"

TEST(kmerMD, kmerMDTest)
{
        KmerOverlap md;
        char nucleotide;

        EXPECT_EQ(md.hasLeftOverlap('A'), false);
        md.markLeftOverlap('A');
        EXPECT_EQ(md.hasUniqueLeftOverlap(nucleotide), true);
        EXPECT_EQ(nucleotide, 'A');
        EXPECT_EQ(md.hasUniqueRightOverlap(nucleotide), false);
        EXPECT_EQ(md.hasLeftOverlap('A'), true);
        md.markRightOverlap('C');
        EXPECT_EQ(md.hasUniqueRightOverlap(nucleotide), true);
        EXPECT_EQ(nucleotide, 'C');
        md.markRightOverlap('G');
        EXPECT_EQ(md.hasUniqueRightOverlap(nucleotide), false);
        EXPECT_EQ(md.hasRightOverlap('C'), true);
        EXPECT_EQ(md.hasLeftOverlap('G'), false);
}
