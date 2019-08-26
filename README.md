# HybridRead

This directory contains python scripts mentioned in the internship report: 'A new algorithm for characterization of hybrid reads in NGS data of HLA genes'.

The script AlignReads.py is for the pre-processing of SAM files and the output text file should be used as input for the script SelectHybridReads.py which is the main algorithm. The input for the script ProcessHybridRead.py are the three files for HLA-A, B and C that contain (1 switch) hybrid read data.

The UnitTests directory contains all unit tests for the main algorithm. Several examples of in- and output files that are used or created by the python scripts can be found in the ExampleInputAndOutputFiles directory.
