# jss

In certain applications, such as document deduplication, finding similar users, etc. calculating all pairwise jaccard similarities can be computationally expensive. For example, calculating all pairwise Jaccard similarities for 100.000 documents would require (100.000*100.000) / 2 = 5.000.000.000 operations.
