# TNW-Reuma.pt

Scripts of the two experiments presented in the work IIEEC_Report_75249.pdf

Experiment 1 - Biological therapy sequences: encoder_sequence.py and temporal_reuma.py

Experiment 2 - DAS28 sequences: encoder_sequence_das.py and temporal_reuma_das.py

For both experiments the first python file consists on the encoding step to transform data from a csv file into temporal sequences and the second python file performs the alignments using the TNW algorithm and the Hierarchical Clustering step using the distance matrix

The user-defined parameters such as the gap penalty, the scoring schema, temporal penalty function, clustering distance measure and the method of transforming similarity matrix into distance matrix are all defined manually in the scripts
