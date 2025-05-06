This is an ongoing project that aims to perform intron and exon identification using a monte-carlo simulation.

Scripts Folder:

- models folder: contains intron and exon length data, and acceptor and donor site pwms (for c. elegans). These are the defaults used in Monte Genie.

- fitness.py: calculates model fitness against existing annotations (.gff3 files) using Manhattan distance

- gen_isoformie.py: runs the simulation using genie and generates isoforms. outputs a textfile named after the gene used

- genie.py: outdated

- genie2.py: class that simulates snRNP and DNA sequence behaviour

- optimize.py: performs optimization of model using fasta and .gff files of desired organism

- test.py: illustrates functionality of genie classes and functions
	
Smallgenes Folder: Contains samples data from C. elegans (fasta and gff files for each gene).
