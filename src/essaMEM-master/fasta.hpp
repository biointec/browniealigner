#ifndef __FASTA_HPP__
#define __FASTA_HPP__

#include <string>
#include <vector>

using namespace std;

void reverse_complement(string &seq_rc, bool nucleotides_only);
void trim(string &line, long &start, long &end);
void load_fasta(string filename, string &S, vector<string> &descr, vector<long> &startpos);

#endif // __FASTA_HPP__

