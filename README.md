## jss - Fast Jaccard Similariy Search Using MinHashing and Locality Sensitive Hashing

In certain applications such as document deduplication, finding similar users, etc. the need arises to calculate all pairwise Jaccard similarities between the available sets/objects. This can be quite computationally expensive. For example, calculating all pairwise similarities between 100.000 documents would require (100.000\*100.000) / 2 = 5.000.000.000 operations. Assuming that one operation requires roughly 130 μs, then we would need 130\*5.000.000.000 = 650.000 sec = 7.5 days.   

A solution is to combine two techniques, namely, MinHashing and Locality Sensitive Hashing (LSH). In order to apply them, it is useful to represent the sets in a term-incidence matrix. For a collection of documents, the term-incidence matrix contains in its rows the k-shingles or terms of the document collection, and the columns correspond to the actual documents. 

MinHashing represents this large and sparse matrix as a "signature", a new smaller matrix, that can easily fit into memory. The signature matrix has the property that for any pair of sets, the probability that the corresponding columns are identical is equal to the Jaccard similarity. LSH is then used to process this signature matrix in a fast way in order to search for candidate sets that might be similar. The Jaccard similarity of those candidate sets is calculated and returned, if the sets are similar.



