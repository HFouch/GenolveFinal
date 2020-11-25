Genolve is a Python program that identifies the most parsimonious series of genetic events that can describe the transition from genome A to genome B.  It can thus be used to investigate the evolutionary history that relates two distinct genomes.  Genolve uses an ammended version of the double-cut-and-join model [1], and extends this approach with a tree searching algorithm to set an upper bound to the minimum number of recombination events required for the observed genomic difference.  Genolve ignores small evolutionary differences such as SNPs and indels, which can effectively be analysed with existing tools.  Genolve provides a list of all the equally most parsimonious solutions.  As input, Genolve requires a file with the chromosome number and the order of all identified synonymous sequence blocks.
The files necessary to run Genolve include Class_extremities_and_adjacencies.py, Class_wrDCJ_Node.py, New_network_wrDCJ.py and commandline_script.py In addition two .txt files containing the the input genomes and a .txt file containing the weighting ratios will be required.
Genolve takes as input two lists of synteny blocks representing the two input genomes.
genA0.txt and genB0.txt are examples of the format in which the input genomes need to be. The synteny blocks confined to separate chromosomes occur on different lines in the file. Synteny blocks on the same chromosome are separated by commas. Synteny blocks in the inverse orientation should be preceded by a negative ('-') sign.
The Weight_ratios.txt file can be altered in accordance with the cost that the user wants to assign to each of the different rearrangements. The order of the weight ratios are inversions, transpositions, balanced translocations, unbalanced translocation, fissions, fusions
Python3 is required to run the commandline_script.py file. The following arguments need to be passed to the program: -target your_target_genome.txt -source your_source_genome.txt -list_of_rearrangement_ratios the_weights_you_want_to_assing_to_the_rearrangements.txt -output name_of_output_file.txt
Example command line:
genolve -target genome_fileB.txt -source genome_fileA.txt -list_of_rearrangements weights.txt –output A_to_B_output.txt

[1]  Yancopoulos S, Attie O, Friedberg R. Efficient sorting of genomic permutations by translocation, inversion and block interchange. Bioinformatics 2005; 21:3340–3346
