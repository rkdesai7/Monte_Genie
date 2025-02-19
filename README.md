This is an ongoing project that aims to perform intron and exon identification using a stochastic process.

The data folder contains different isoform information, as well as position weight matrices that are used to train the model.

The scripts folder contains the scripts to run the model, as well as various fasta files to use for training and testing:

  - genie.py - the script containing the different class variables
  - test.py - a tester illustrating the functionalities of genie
  - main.py - uses genie to perform the identification
  - fitness.py - calculates the accuracy of the model
  - optimize.py - a genetic algorithm used to optimize model parameters
  - /smallgenes - contains .gff3 and .fa files for training and testing
