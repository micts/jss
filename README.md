## jss - Fast Jaccard Similariy Search Using MinHashing and Locality Sensitive Hashing

In certain applications such as document deduplication, finding similar users, etc. the need arises to calculate all pairwise Jaccard similarities between the available sets/objects. This can be quite computationally expensive. For example, calculating all pairwise similarities between 100.000 documents would require (100.000\*100.000) / 2 = 5.000.000.000 operations. Assuming that one operation requires roughly 130 μs, then we would need 130\*5.000.000.000 = 650.000 sec = 7.5 days.
