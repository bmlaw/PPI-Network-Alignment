# sana_testing

## Overview
The goal of this project is to test the Simulated Annealing Network Aligner (SANA) for Protein Protein Interaction Networks (PPINs) by examining how it utilizes graph topology and sequence similarity under different inputs of ideal Edge Correctness (EC) and Sequence Simililarity (Seq). 

So far, we have written a network constructor (get_el.py), a protein sequence file consrtuctor (PPI_get_FASTA.py), and a experimental data viewer (sana_read_out_files.py). We are now planning to run SANA multiples times at different Seq inputs and check the output sequence similarity scores. We will also list the Reciprocal Best Hits (RBHs) for each SANA output alignment file.

## File preparation
### Network constructor: get_el.py
The program constructs network files based on the BIOGRID mitab files.

1st input: the directory of the BIOGRID mitab file. In this research, we used data from BIOGRID-ORGANISM-4.4.210.mitab. 

1st output: a preliminary version of the network file (.el): each line of which consists of both the name and the gene id of each protein involved in the interaction in the form of name_geneID (eg. Bcar1_25414	Pik3r1_25513). 

2nd output: a gene-id-only version of the first output. 

For example, if we want to get the network of rat (RNorvegicus), we run the following command on the terminal:
`python get_el.py RNorvegicus.mitab RNorvegicus.el RNorvegicus_n.el`

Note:
When constructing each network, the program applies several filters so that 
1. self-loops (interaction between the same types of proteins) will be excluded
2. interspecies interactions will be excluded
3. only physical interactions will be included
4. only experimentally detected interactions will be included



### FASTA file constructor: PPI_get_FASTA.py
The program will try to find and record the sequences of all proteins in the PPIN. Before running this program, the network constrcutor must be run. If the sequence of either protein involved in an interaction could not be found in the NCBI database, the interaction will not be included in the updated "accurate" network files. If sequence similarity parameter is not going to be 0 in a SANA run, use the network file updated by this program!

1st input: the directory of the .el file written by the network constructor program.

1st output: directory of the output .fasta file (contains the protein sequences found online).

2nd output: directory of the updated network file (gene-id-only version)

3rd output: directory of the updated network file (contains both protein name and gene id)

For example, if we want to get the fasta file and the updated network of rat (RNorvegicus), we run the following command on the terminal:
`python PPI_get_FASTA.py RNorvegicus.el RNorvegicus.fasta RNorvegicus_accurate_n.el RNorvegicus_accurate.el`

### BLAST files
To compare the sequence similarity between proteins in the two networks to be aligned, SANA uses the BLAST output. Hence, if the "ideal" sequence similarity input of a SANA run is not 0, then the relevant blast file must be prepared in advance. In this project, BLAST was run in parallel on the IWU cluster.

Below is an example of the shell script (.sh) to run BLAST:

```
#!/usr/bin/env bash

#
# Execute the job and place stdout and stderr files in the current working directory
#$ -cwd
#
# Max runtime requested
#$ -l h_rt=30:30:00
#
# Job name
#$ -N HSapiens_blast1


./blast/makeblastdb -in FASTA/HSapiens.fasta -dbtype prot -out HSapiens

time ./blast/blastp -query FASTA/MMusculus.fasta -db HSapiens -out scores/MMusculus_HSapiens_blast.out -outfmt 6 -word_size 6 -evalue 10 -num_threads 16 -max_hsps 1
time ./blast/blastp -query FASTA/RNorvegicus.fasta -db HSapiens -out scores/RNorvegicus_HSapiens_blast.out -outfmt 6 -word_size 6 -evalue 10 -num_threads 16 -max_hsps 1
time ./blast/blastp -query FASTA/SPombe.fasta -db HSapiens -out scores/SPombe_HSapiens_blast.out -outfmt 6 -word_size 6 -evalue 10 -num_threads 16 -max_hsps 14
```

In this example, we made a database based on the human (HSapiens) fasta file using makeblastdb. Then, we run blast 3 times by taking mouse (MMusculus), rat (RNorvegicus), and yeast (SPombe) as queries. 

To submit the shell scripts so that things can be run in parallel, use the qsub command. For example:
`qsub run_blast_HSapiens.sh`

To check if your submitted job is actually running, use the following command:
`qhost -j`

Just a personal naming protocol: the species name in the name of the shell script is the name of the database species.

## Run SANA
SANA can be downloaded from Prof.Wayne Hayes's github repo: [Wayne Hayes SANA](https://github.com/waynebhayes/SANA)

For details not covered in this file, refer to the tutorial included in Wayne's repo.

Before running SANA, make sure the file preparation processes described have been completed.

Here is an example to run SANA (run it under SANA directory):
`./sana -fg1 RNorvegicus.el -fg2 HSapiens.el -ec 0.3 -seq 0.7 -t 100`

Note:

1. -fg1 is the edge list file of the query species, while -fg2 is the el file for the database species. The size (number of interactions and/or proteins involved ???) of the query network must be smaller or equal to the database network. Otherwise, SANA will report error.

2. Given an "ideal" sequence similarity input -seq that is greater than 0, SANA will look up for sequence similarity files automatically in the sequence/scores directory. Hence, a good naming protocol is important! For example, the query network is RNorvegicus.el while the database is HSapiens.el, so SANA will go to sequence/scores for RNorvegicus_HSapiens_blast.out

3. Usually, make sure ec + seq = 1

4. -t is in minutes. As stated by Wayne Hayes, for networks as large as HSapiens,  around 64 min of runtime will attain a high alignment precision, while this growth will increase by a few percent if we increase the runtime to 120 min. 

