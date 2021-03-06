

# Short description:

browniealigner alings short reads to a reference genome which is an important task in many genome analysis pipelines. It is based on a branch and bound alignment algorithm that uses the seed-and-extend paradigm to accurately align short Illumina reads to a graph. Given a seed, the algorithm greedily explores all branches of the tree until the optimal alignment path is found. To reduce the search space it computes the upper bounds to the alignment score for each branch and discards the branch if it cannot improve the best solution found so far. Additionally, by using a two-pass alignment strategy and a higher-order Markov model, paths in the de Bruijn graph that do not represent a subsequence in the original reference genome are discarded from the search procedure.

#  Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

#  Prerequisites

This package requires a number of packages to be install on your system. Required: CMake; Google's Sparsehash; gcc (GCC 4.7 or a more recent version) Optional: ZLIB; Googletest Unit Testing

How to install these packages:

As a root, execute the following commands:

on Redhat / Fedora distributions

    yum install cmake
    yum install sparsehash-devel
    yum install zlib-devel (optional)

on Ubuntu / Debian distributions

    aptitude install cmake
    aptitude install libsparsehash-dev
    aptitude install libghc-zlib-dev (optional)


# Installing browniealigner:

The installation is now simple. First, clone browniealigner from the Github address

    git clone "https://github.com/biointec/browniealigner.git"

From this directory, run the following commands:

    mkdir build
    cd build
    cmake ..
    make install

By executing ./brownie you will see

    brownie: no command specified
    Try 'brownie --help' for more information

# Usage:
brownieAligner aligns short Illumina reads to the de Bruijn graphs. To do this you need first builds the graph based on an input data (could be a reference genome or other Illumina dataset for example)

     brownie index [options] file1 [file2]...
Then align your reads to the graph with this command:

     brownie align [options] [file_options] file1 [[file_options] file2]...
options:

    -h    --help    display help page
    options arg
    -k      --kmersize            kmer size [default = 31]
    -t      --threads             number of threads [default = available cores]
    -e      --essa                sparseness factor of index structure [default = 1]
    -nBB    --noBranchAndBound    do not use branch&Bound pruning
    -nMM    --noMarkovModel       do not use Markov Model
    -nMEM   --noEssaMEM           do not use EssaMEM
Using -nBB is not recommended at all. However, with -nMM browniealigner runs almost twice faster with a little drop in the accuracy. Using -nMEM also not recommended because it reduces the accuracy especially if your kmer size is large.  However, by this command browniealigner  runs faster with less memory demand.

     file_options:

    -o    output    aligned output read file name [default = inputfile.corr]
Examples :

     brownie index -k 31 -t 4 -p tempDir genome.fasta
     brownie align -k 31 -t 4 -p tempDir -o output.fastq input.fastq -nMM
 
The output file contains aligned reads which are the corresponding subsequences from the reference sequence extracted from the graph after the alignment. In the tempDir directory, you can find nodes of the graph in nodes.stage2 file. There is also an ncf file with the same name as the output file which contains the alignment path for each read. In this file for every read, there is a tab separated line like this:

    >SRR1206093.37/1  A(1)A   S(251)S B(0)B   P(42,33957)P     C(-200 -927 -137)C
    >SRR1206093.37/2  A(1)A   S(250)S B(0)B   P(33818,117)P    C(137)C
1.  The first part shows the read ID (here in this example we show the pair forward and reverse read).
2.  The second part A(*)A shows if the read aligned (A(1)A) to the graph or not (A(0)A).
3.  The third part S(*)S, is the similarity score after alignment.
4.  The fourth part B(*)B shows if the read is aligned with the branch and bound algorithm (B(1)B) or not (B(0)B).
5.  The fifth part P(*)P shows the offset position of the node that the read aligns to (for both the forward and the reverse complement).
6.  The last part C(*)C, is the chain of nodes that the read aligns to. Nodes with the negative sign are reverse complement of those stored in nodes.stage2.

# Some Implementation Details:
Building graph:

To reduce memory requirements, k-mers are encoded by 2k bits and stored in a memory-efficient hash map
with only 2 bits overhead per entry. Overlap between k-mers is encoded by 8 bits: 4 bits to indicate if the k-mers can be left-extended with A, C, T or G and similarly 4 bits to represent right overlap.

Indexing:

To store k-mers and retrieve them later in an efficient way, we use google sparsehash.  Per each k-mer, we store the node that it occurs and also the offset that it begins.  However, it is possible that a read does not have a  valid k-mer the can be found in this table (i.e., if the read has too many errors such that you cannot find a single substring of length k without any error ). In this case, we use essaMEM. We first concatenate the sequences in all the nodes into a single string. The content of each node is followed by its NodeID. It makes it possible to find the corresponding node that a given pattern occurs within. 


Alignment:

The reads are aligned back to the DBG using a seed-and-extend paradigm. In case a read contains at least one true k-mer, this k-mer is used as a seed that uniquely maps the read to a certain node in the DBG. A depth-first search (DFS) on the graph is performed to align both ends of the read beyond the seed(s). Pairwise alignments are used to find the optimal alignment path. Branch-and-bound conditions are used to limit the search space. For each branch, an upper bound is computed to the alignment score that could be obtained in that branch. The branch is discarded from the search procedure if it cannot improve the best solution found so far. In order to rapidly find candidate solutions with a high score, the DFS greedily prioritizes towards the node that appears best. 

Every time we need to compute the current similarity score of an aligning a read to some branches in the graph, we don't' compute the similarity score of the entire read with the whole branch, we only compute the similarity score between a part of a read which is not aligned yet with the corresponding node (the last node in that branch). Therefore, the similarity score of an aligning a read to a branch in the graph can be obtained by the summing over the similarity scores between parts of reads and their corresponding nodes in the graph. This avoids the quadratic number of computation of similarity score.

Markov Model:

To drive an n-Markov Table (1<n<=maxorder), we first align the reads to graph one by one. In the second step, for the reads that align to more than two nodes, the node chains are extracted. To avoid propagating potential errors to the model that arise as a result of a wrong alignment, only use seeds, i.e., perfect matches between (a part of) a read and the sequence implied by the path in the de Bruijn graph. Each seed is a chain of m (m>0) nodes. However, chains with less than 3 nodes are discarded. For the rest, all subchains in different size are extracted. Each subchain consists of a head (the rightmost node), the current node ( the second rightmost node) and tail ( rest of that subchain). An n-Markov transition table comprises N entries where N is the number of nodes in the graph. The current entry of this table keeps all the sub-chains with the same current node. in each entry, the frequency of observing each head is stored per tail. The sum of all frequencies is equal to the observed coverage of that chain. In BrownieAligner we are merely interested to find the transitions in this table that are wrong. These transitions are due to the wrong alignment in the first step.  For example, imagine two paths with the same tail and current node, the frequency of seeing a particular head node is 1 and the frequency of seeing another head node is 100. We can maybe conclude that seeing the second head is wrong and should be avoided. The decision is made based on the coverage, length of that chain and the frequency of seeing a particular head. Please look at the paper for more information.

# Report bugs 
Please report bugs to : Mahdi.Heydari@UGent.be
 

# Citation
Please cite our paper :
@Article{Heydari2018,
author="Heydari, Mahdi and Miclotte, Giles and Van de Peer, Yves and Fostier, Jan",
title="BrownieAligner: accurate alignment of Illumina sequencing data to de Bruijn graphs",
journal="BMC Bioinformatics",
year="2018",
month="Sep",
day="04",
volume="19",
number="1",
pages="311",
}

