#include "readfile/fastafile.h"
#include "readfile/fastqfile.h"
#include "readfile/rawfile.h"
#include "readfile/samfile.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <string>

using namespace std;

// ============================================================================
// SAM FORMAT
// ============================================================================

TEST(readFile, SAMGZTest)
{
        ReadFile *file = new SamFile(true);
        file->open("test.sam.gz");

        string target("TAATCCCCGCCAAATTCGTGACCTGTCATTCGTCG");
        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "CGACAGGGATAGTGTAGCTGACCGTTGTGACTGGC";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 10){
                        EXPECT_EQ(read, target);
                        
                }
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, SAMTest)
{
        ReadFile *file = new SamFile(false);
        file->open("test.sam");

        string target("TAATCCCCGCCAAATTCGTGACCTGTCATTCGTCG");
        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "CGACAGGGATAGTGTAGCTGACCGTTGTGACTGGC";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 10){
                        EXPECT_EQ(read, target);
                        
                }
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, SAMCopyTest)
{
        ReadFile *in = new SamFile(false);
        in->open("test.sam");

        ReadFile *out = new SamFile(false);
        out->open("test.copy.sam", WRITE);

        ReadRecord record;
        while (in->getNextRecord(record)) {
                out->writeRecord(record);
        }

        in->close();
        delete in;

        out->close();
        delete out;
}

// ============================================================================
// FASTA FORMAT
// ============================================================================

TEST(readFile, FastAGZTest)
{
        ReadFile *file = new FastAFile(true);
        file->open("test.fa.gz");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");
        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 10){
                        EXPECT_EQ(read, target);
                        
                }
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;

        file = new FastAFile(true);
        file->open("test.fa.gz");

        target = "GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC";
        ReadRecord record;
        file->getNextRecord(record);

        EXPECT_EQ(record.getRead(), target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        numReads = 1;
        while (file->getNextRecord(record)) {
                numReads++;

                if (numReads == 10){
                        EXPECT_EQ(record.getRead(), target);
                        
                }
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, FastATest)
{
        ReadFile *file = new FastAFile(false);
        file->open("test.fa");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");
        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 10){
                        EXPECT_EQ(read, target);
                        
                }
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, FastACopyTest)
{
        ReadFile *in = new FastAFile(false);
        in->open("test.fa");

        ReadFile *out = new FastAFile(false);
        out->open("test.copy.fa", WRITE);

        ReadRecord record;
        while (in->getNextRecord(record)) {
                out->writeRecord(record);
        }

        in->close();
        delete in;

        out->close();
        delete out;
}

// ============================================================================
// FASTQ FORMAT
// ============================================================================

TEST(readFile, FastQTest)
{
        ReadFile *file = new FastQFile(true);
        file->open("test.fastq.gz");

        string target("ATATAGATGTACATAAATTAGTTGAAGTATATGAACG");
        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "TAGGAAAGCGAAGCCATTCAATACGAAGTATTGTATA";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 9){
                        EXPECT_EQ(read, target);
                        
                }
        }

        EXPECT_EQ(numReads, 9);

        file->close();
        delete file;
}


TEST(readFile, FastQGZTest)
{
        ReadFile *file = new FastQFile(false);
        file->open("test.fastq");

        string target("ATATAGATGTACATAAATTAGTTGAAGTATATGAACG");
        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "TAGGAAAGCGAAGCCATTCAATACGAAGTATTGTATA";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 9){
                        EXPECT_EQ(read, target);
                        
                }
        }

        EXPECT_EQ(numReads, 9);

        file->close();
        delete file;
}

TEST(readFile, FastQCopyTest)
{
        ReadFile *in = new FastQFile(false);
        in->open("test.fastq");

        ReadFile *out = new FastQFile(false);
        out->open("test.copy.fastq", WRITE);

        ReadRecord record;
        while (in->getNextRecord(record)) {
                out->writeRecord(record);
        }

        in->close();
        delete in;

        out->close();
        delete out;
}

// ============================================================================
// RAW FORMAT
// ============================================================================

TEST(readFile, RawGZTest)
{
        ReadFile *file = new RawFile(true);
        file->open("test.raw.gz");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");
        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 10){
                        EXPECT_EQ(read, target);
                }
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, RawTest)
{
        ReadFile *file = new RawFile(false);
        file->open("test.raw");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");

        string read;
        file->getNextRead(read);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read)) {
                numReads++;

                if (numReads == 10){
                        EXPECT_EQ(read, target);
                        
                }
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, FastRawTest)
{
        ReadFile *in = new RawFile(false);
        in->open("test.raw");

        ReadFile *out = new RawFile(false);
        out->open("test.copy.raw", WRITE);

        ReadRecord record;
        while (in->getNextRecord(record)) {
                out->writeRecord(record);
        }

        in->close();
        delete in;

        out->close();
        delete out;
}
