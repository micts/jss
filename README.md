## jss - Fast Jaccard Similariy Search Using MinHashing and Locality Sensitive Hashing

In certain applications such as document deduplication, finding similar users, etc. the need arises to calculate all pairwise Jaccard similarities between the available sets/objects. This can be quite computationally expensive. For example, calculating all pairwise similarities between 100.000 documents would require (100.000\*100.000) / 2 = 5.000.000.000 operations. Assuming that one operation requires roughly 130 μs, then we would need 130\*5.000.000.000 = 650.000 sec = 7.5 days.   

A solution is to combine two techniques, namely, MinHashing and Locality Sensitive Hashing (LSH). In order to apply them, it is useful to think of the set collection as a term-incidence matrix. For a collection of documents, the term-incidence matrix contains in its rows the k-shingles or terms of the document collection, and the columns correspond to the actual documents. The matrix contains 1 if the k-shingle/term occurs in the document, and 0 otherwise.

MinHashing represents this large and sparse matrix as a "signature", a new smaller matrix, that can easily fit into memory. The signature matrix has the property that for any pair of sets, the probability that the corresponding columns are identical is equal to the Jaccard similarity. LSH is then used to process this signature matrix in a fast way in order to search for candidate sets that might be similar. Finally, the Jaccard similarity of those candidate sets is calculated. For a detailed explanation of MinHashing and LSH you can refer to the book [Mining of Massive Datasets](http://www.mmds.org/).

## How to use

Requirements: ```numpy``` and ```Python 3.x```

Clone the repository:
```git clone https://github.com/micts/jss.git```

## Examples

Refer to the ```jss_example.py``` for a short example.

## TO DO

* Check compatibility with python 2.X
* Complete the example in ```jss_example.py```

