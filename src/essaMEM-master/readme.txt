essaMEM
-------------------
  
Installing:
	1 go to the folder that contains the essaMEM.tgz file
	type the following:
	tar -xvzf essaMEM.tgz
	cd essaMEM
	make
  
Usage:  ./mummer [options] <reference-file> <query-file1, query-file2, ...>

Implemented MUMmer v3 options:
  -mum           	compute maximal matches that are unique in both sequences
  -mumreference  	compute maximal matches that are unique in the reference-
                 		sequence but not necessarily in the query-sequence (default)
  -mumcand       	same as -mumreference
  -maxmatch      	compute all maximal matches regardless of their uniqueness
  -l             		set the minimum length of a match
                 		if not set, the default value is 20
  -b             		compute forward and reverse complement matches
  -F             		force 4 column output format regardless of the number of
                 		reference sequence inputs
  -n             		match only the characters a, c, g, or t
  -L             		print length of query sequence in header of matches
  -r             		compute only reverse complement matches
  -s             		print first 53 characters of the matching substring
  -c             		report the query position of a reverse complement match relative to the forward strand of the query sequence
  
  Additional options:
  -k             		sampled suffix positions (one by default)
  -threads       		number of threads to use for -maxmatch, only valid k > 1 
  -qthreads      	number of threads to use for queries 
  -suflink       		use suffix links (1=yes or 0=no) in the index and during search [auto]
  -child         		use child table (1=yes or 0=no) in the index and during search [auto]
  -skip                        sparsify the MEM-finding algorithm even more, performing jumps of skip*k [auto (l-10)/k]
  			this is a performance parameter that trade-offs SA traversal with checking of right-maximal MEMs
  
  Example usage:
  
  ./mummer -maxmatch -l 20 -b -n -k 3 -threads 3 ref.fa query.fa
  Find all maximal matches on forward and reverse strands
  of length 20 or greater, matching only a, c, t, or g.
  Index every 3rd position in the ref.fa and use 3 threads to find MEMs.
  Fastest method for one long query sequence.
  
  ./mummer -maxmatch -l 20 -b -n -k 3 -qthreads 3 ref.fa query.fa
  Same as above, but now use a single thread for every query sequence in
  query.fa. Fastest for many small query sequences.