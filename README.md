This is an ongoing project that aims to perform intron and exon identification using a stochastic process.

Scripts Folder: contains the scripts to run the model, as well as various fasta files to use for training and testing:

  - genie.py - the script containing the different class variables
  - test.py - a tester illustrating the functionalities of genie
  - main.py - uses genie to perform the identification
  - fitness.py - calculates the accuracy of the model given the gene name (ie. ch.1_0)
  - optimize.py - a genetic algorithm used to optimize model parameters
  - /models - contain pwm and intron length data used

smallgenes: contains many .gff3 and .fa files used for training and testing
