#!/bin/bash

#CLASS CAMB 7140
#DATE: Feb 12, 2025
###Lab 5: Dumpster diving in RNA-seq data

##Task 1: 

:<<'TASK1'
Since you have no idea where these unmapped reads are coming from, you decide that a broad and unbiased approach is your best bet to identify potential non-human reads in your sample. 
One strategy is to use Sourmash to create a minHash ‘sketch’ of your fastq file and compare this sketch against a reference set of sketches. 
To begin this task, you will need to download a set of about 90,000 microbial genomes here (don’t worry, it’s a small file). 
Here’s some example code to get you started. Be sure you run this in the appropriate Conda environment.
TASK1

#ran all code in terminal 

#time = ~2min to sketch ~9M reads
sourmash sketch dna -p scaled=10000,k=31,abund SRR8668774_dehosted.fastq.gz --name-from-first

#outputs a file with the same name as the input with .sig ending

# time = ~2min
sourmash gather -k 31 SRR8668774_dehosted.fastq.gz.sig genbank-k31.lca.json.gz
# once this is done, try rerunning with an additional argument to relax the threshold used for classification: '--threshold-bp 100'

:<<'OUTPUT'
== This is sourmash version 4.8.14. ==
== Please cite Irber et. al (2024), doi:10.21105/joss.06830. ==

selecting specified query k=31
loaded query: SRR8668774.33 33/1... (k=31, DNA)
--
loaded 93249 total signatures from 1 locations.
after selecting signatures compatible with search, 93249 remain.

Starting prefetch sweep across databases.
Prefetch found 17 signatures with overlap >= 50.0 kbp.
Doing gather to generate minimum metagenome cover.

overlap     p_query p_match avg_abund
---------   ------- ------- ---------
80.0 kbp       0.0%    0.4%       1.1    FR798975.1 Leishmania braziliensis MHOM/BR/75/M2904 complete genome, chromosome 1
70.0 kbp       0.3%    1.2%      13.9    FLZI01001256.1 Escherichia coli isolate CG05C.C2 genome assembly, contig: 20151013_NOD...
found less than 50.0 kbp in common. => exiting

found 2 matches total;
the recovered matches hit 0.3% of the abundance-weighted query.
the recovered matches hit 0.5% of the query k-mers (unweighted).
OUTPUT