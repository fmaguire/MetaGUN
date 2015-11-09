Rehosted
========

Codebase rehosted by Finlay Maguire 

Links in paper are now broken so rehosting to prevent
orphaning


```
README file of MetaGUN Release Version 1.0
======================================================================
      Copyright Deparment of BioMedical Engineering,
      College of Engineering, Peking University, Beijing
======================================================================

INTRODUCTION:

 - MetaGUN is a gene prediction protocol for metagenomic fragments 
   based on a machine learning approach of SVM. It implements in a 
   three-stage strategy to predict genes. Firstly, it classifies 
   input fragments into phylogenetic groups by a k-mer based sequence 
   binning method. Then, protein-coding sequences are identified for 
   each group independently with SVM classifiers that integrate entropy 
   density profiles (EDP) of codon usage, translation initiation (TIS) 
   scores and open reading frame (ORF) length as input patterns. Finally, 
   the TISs are adjusted by employing a modified version of MetaTISA.
 
 - This package includes the following files and directories:
   + src -- The C++ source codes
   + bin -- The executable files
   + dat -- The parameter files 
     + binmodel -- The parameter files for fragments classification 
     + cdsmodel -- The SVM classifiers for CDS identification 
     + tismodel -- The supervised TIS parameter files 
     + Cdd -- The BLAST database of conserved domains
   + scripts -- The Python scripts
     + metagun.py -- The overall MetaGUN prediction pipeline
     + bin-seqs.py -- A script for running the subprocess of fragments 
       classification separately
     + domain-search.py -- A script for running the subprocess of domain
       searches separately
   + example -- An example FASTA sequence file for testing

======================================================================

INSTLLATION:

 - Requirements
   We built the MetaGUN pipeline by using Python scripts. So, the 
   installation of Python programming language is required before 
   running MetaGUN program. A relese version of Python2 after 
   Python2.6 is needed. You can get Python from its official website:
       http://www.python.org/
       
 - LINUX
   Run ./build script
   
 - WINDOWS
   Uncompress the win32_binary.rar and copy all the executable files 
   to the 'bin' directory, then run convert2dos.cdsmodels.bat to convert
   the SVM classifiers into DOS format.

======================================================================

IMPLEMENTATION:

 - Generally, simply run metagun.py to predict genes for metagenomic
   fragments by defining the project name and the path of the input
   sequences. Like running on the example data we provided:
       ./scripts/metagun.py -n example -s example/example.seq
   Prediction will be implemented on example.seq, and the final result
   saves in example.metagun
 
 - For convenience, the subproccesses of the fragment classification 
   and the domain searches can be run separately. We provide two scripts
   namely bin-seqs.py and domain-search.py in the 'scripts' directory.
   Taking the provided 'example.seq' as an example, you can firstly 
   obtain the fragment classification results,
       ./scripts/bin-seqs.py example/example.seq example.bin
   and the domain search results,
       ./scripts/domain-search.py example/example.seq example.blastcdd
   separately. Then, do gene prediction by specifying the fragment
   classification results with -b option, and the domain search results
   with -B option,
       ./scripts/metagun.py -n example -s example/example.seq 
                            -b example.bin -B example.blastcdd
   Prediction will be saved in file example.metagun
   
 - The prediction results is the 'MED' format, like,
   >Read01
       1     144   +
     791     868   +
     202     765   -
   >Read02
      88     870   +
       2      82   -
   >Read03
       3     869   -
```
