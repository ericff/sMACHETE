# sMACHETE
Code used for running the sMACHETE, for the portion after running the
MACHETE on the Seven Bridges Cancer Genomics Cloud (CGC), for the paper

Precise, pan-cancer discovery of gene fusions reveals a signature of
selection in primary tumors, by Donald Eric Freeman, Gillian Lee Hsieh,
Jonathan Michael Howard, Erik Lehnert, and Julia Salzman.

There are two key scripts that run the other scripts:
(i) setup_SBT_queries.pl
After one downloads outputs from many MACHETE runs, this script is
run. It generates fasta files to be used as inputs for running Sequence
Bloom Tree (SBT) queries on the CGC.

(ii) process_SBT_query_output.pl	
After one runs SBT queries and downloads the output from the queries,
this script calls other scripts to process the output to generate many
of the results and figures in the paper.
