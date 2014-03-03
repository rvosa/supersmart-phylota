SUPERSMART example case: Primate phylogeny
==================================

Here we describe a test case of how to use the SUPERSMART
pipeline to infer a phylogeny for the taxonomic order Primates.
First, an example is given of how to invoke the entire pipeline
to build a primate tree using a single command. Below, we will describe
the input and output of each component of the pipeline in 
detail.

Running the primate example 
==================================

The script 'run.sh' is located in the top-level
directory of the SUPERSMART source and combines the invokation of
all intermediate steps in the pipeline for the primate example. 
Note that the environment variable SUPERSMART_HOME must be set 
to the directory which contains 'run.sh' and the 'examples/' 
subdirectory. From there, simply call

$ sh run.sh

After all steps of the pipeline have finished running, the file
'megatree.dnd' is written to the directory containing the primate example 
data ($SUPERSMART_HOME/examples/primates). This file contains the 
derived tree and can be viewed with common analysis tool in phylogeny
(e.g. FigTree).

Documentation of all intermediate steps
==================================
Input for SUPERSMART is 
a text file containing the names of all species for which the 
final tree will be derived. For the primate example, this
is the file 'names.txt' located in the primate example directory, 
which contains the species names of 217 primates. For all intermediate
steps in the pipeline, there exist scripts in the 'scripts/supersmart' 
directory of the project. Thus, each step can also be run individually.
Below, all individual steps are documented. 

===
Step 1 
===
The first step is to search the given species names in the phylota
database, to extract all taxonmic ranks for all species and write 
them into a table. 
This is accomplished by running the script 'mpi_write_taxa_table.pl'
located in the 'scripts/supersmart/' directory. The script writes
a table with all NCBI taxon ids for all taxonomic ranks for the species to STDOUT.
Below we show the table for the first three species in the example.

name    species genus   family  order   class   phylum  kingdom
Allenopithecus nigroviridis     54135   54134   9527    9443    40674   7711    33208
Allocebus trichotis     122248  122247  30615   9443    40674   7711    33208
Alouatta belzebul       30590   9499    378855  9443    40674   7711    33208

===
Step 2
===
Next, a 'common tree' according to the NCBI taxonomy is created from all taxa present 
in the table generated in step 1. The tree is written by invoking the script 
'write_common_tree.pl' with the species table as input. By default, the common tree
is written to STDOUT in the newick format.

===
Step 3
===  
For the taxons present in the species table generated in step 1, markers 
that contain phylogenetically relevant information are automatically selected. These 
marker sequences are collected from a local copy of the PhyLoTa database 
(http://phylota.net/). 

The PhyLoTa database contains publicly available nucleotide
sequences that are taxonomically organized and, at strategically chosen nodes in the
taxonomy, all-versus-all BLAST searches were conducted for representative sequences
for the specific node. The BLAST hits were then clustered around a representative
'seed sequence'. 
 
The script 'mpi_write_alignments' takes as input the species table generated in 
step 1 and queries SUPERSMART's local copy of the PhyLoTa database for relevant
sequence clusters for available markers. The sequences in each selected cluster are
then aligned and the multiple sequence alignment is written to a FASTA file named
according to the GI code of the seed sequence of the PhyLoTa cluster. 
In the example described here, 3059 files are written to the 'examples/primates'
directory. The file names of the produced FASTA files are written to STDOUT.

===
Step 4
===

    

 





