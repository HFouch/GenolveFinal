# GenolveFinal
The files nessisary to run Genolve include Class_extremities_and_adjacencies.py, Class_wrDCJ_Node.py, New_network_wrDCJ.py and commandline_script.py

Genolve takes as input two lists of synteny blocks representing the two input genomes. 

genA0.txt and genB0.txt are examples of the format in which the input genomes need to be. 
The synteny blocks confined to seperate chromsomes occur on different lines in the file. Synteny blocks on the same chromosome are seperated by commas. 
Synteny blocks in the inverse orientation should be preceded by a negative ('-') sign.

The Weight_ratios.txt file can be altered in accordance with the cost that the user wants to assign to each of the different rearrangements
  the order of the weight ratios are inversions, transpositions, balanced translocations, unbalanced translocation, fissions, fusions
  
Python3 is required to the commandline_script.py file. The following arguments need to be passed to the program:
-target your_target_genome.txt -source your_source_genome.txt -list_of_rearrangement_ratios the_weights_you_want_to_assing_to_the_rearrangements.txt -output name_of_output_file.txt
