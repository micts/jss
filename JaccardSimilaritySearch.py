import numpy as np
from sys import stdout

class MinHashingLSH:

    def __init__(self, data, threshold=0.5, n_bands=22, n_rows=6, min_set_size=None, random_state=None, verbose=0):
        
        """
        Parameters
        ----------
        data: ndarray 
              Array of two columns representing the object and item in the first and second column, respectively,
              encoded as integers. Zero-based encoding is a requirement. The objects correspond to the different 
              sets or vectors for which the similarity is to be computed. The items are the elements of those sets 
              or vectors. The ndarray could alternatively be the row and column indices of a sparse item-object matrix.
        
        threshold: float, default: 0.5
                   Objects with Jaccard Similarity greater than the defined threshold are considered similar.
        
        n_bands: int, default: 22
                 The number of bands, where n_bands * n_rows = n_hash_functions.
                         
        n_rows: int, default: 6
                 The number of rows per band, where n_bands * n_rows = n_hash_functions.
        
        min_set_size: int, default: None
                 Discard all objects with nunber of items less than min_set_size. That is, discard all sets with number 
                 of elements less than min_set_size.
        
        random_state: int or None, default: None
                      The seed of the pseudo random number generator to use when creating the signature matrix 
                      in MinHashing. If int, random_state is the seed used by the random number generator. 
                      If None, the random number generator is the RandomState instance used by np.random.              
        
        verbose: int, default: 0
                 Any number greater than 0 for verbosity.
                 
        Attributes
        ----------

        data: ndarray
              Same as in parameters. If min_set_size is not None, then the data attribute does not contain the
              objects with less number of items than the specified min_set_size.
              
        n_hash_functions: int
              The number of hash functions use in MinHashing, where n_hash_functions = n_bands * n_rows.
        """
        
        self.data = data
        self.threshold = threshold
        self.n_bands = n_bands
        self.n_rows = n_rows
        self.min_set_size = min_set_size
        self.n_hash_functions = self.n_rows * self.n_bands
        self._random_state = random_state
        self._verbose = verbose      
          
        if self._verbose > 0:
            print('Parameters',
                  '----------',
                  'Nr. of Hash Functions: ' + str(self.n_hash_functions),
                  'Nr. of Bands: ' + str(self.n_bands),
                  'Nr. of Rows per Band: ' + str(self.n_rows),
                  'Threshold: ' + str(self.threshold),
                  '----------',
                  sep='\n')
        
        self._n_items = np.max(data[:, 0] + 1)      
        self._n_properties = np.max(data[:, 1] + 1)
        item_bins = np.bincount(data[:, 0])
        self._item_inds = np.cumsum(item_bins)
        self._item_inds = np.hstack((0, self._item_inds))
        self.data = data[data[:, 0].argsort()]
        
        if min_set_size is not None:
            items_to_keep = np.where(item_bins > self.min_set_size)[0]
            inds_to_keep = np.in1d(self.data[:, 0], items_to_keep)
            self.data = self.data[inds_to_keep , :]
            item_bins = np.bincount(self.data[:, 0])
            self._item_inds = np.cumsum(item_bins)
            self._item_inds = np.hstack((0, self._item_inds))
        
        self._unique_items = np.unique(self.data[:, 0])
        
    def _hash_to_buckets(self, band):
        k = self.n_rows * band
        buckets = {}
        for item in self._unique_items:
            bucketID = tuple(self.signature_matrix[k:(k + self.n_rows), item])  
            if bucketID in buckets:  
                buckets[bucketID].append(item)
            else:
                buckets[bucketID] = []  
                buckets[bucketID].append(item) 
        return(buckets)        
    
    def jaccard_similarity(self, set1, set2):
        return(len(set1.intersection(set2)) / len(set1.union(set2)))
    
    def _compare_candidates(self, buckets, bucketID):
        sig_matrix_candidates = self.signature_matrix[:, buckets[bucketID]]
        n_candidates = sig_matrix_candidates.shape[1]
        bucket_array = np.asarray(buckets[bucketID]) 
    
        for candidate in range(n_candidates - 1):
            temp = sig_matrix_candidates[:, (candidate + 1):n_candidates]
            temp_bucket = bucket_array[(candidate + 1):n_candidates]
            item_signature = np.tile(sig_matrix_candidates[:, candidate], (temp.shape[1], 1)).T
            estimated_sims = np.sum(item_signature == temp, axis = 0) / self.n_hash_functions
            candidates = np.where(estimated_sims >= self.threshold)[0]
            if len(candidates) > 0:
                candidates = temp_bucket[candidates]
                item = bucket_array[candidate]
                item_set = set(self.data[self._item_inds[item]:self._item_inds[item + 1]][:, 1])
                for candidate in candidates:
                    if candidate in self._similar_items_dict[item]:
                        pass
                    else:
                        candidate_set = set(self.data[self._item_inds[candidate]:self._item_inds[candidate + 1]][:, 1])
                        jaccard_similarity = self.jaccard_similarity(item_set, candidate_set)
                        if jaccard_similarity >= self.threshold:
                            self._similar_items_dict[item].append(candidate)
                            self.jaccard_similarity_score.append(jaccard_similarity)
                            self.similar_items.append([item, candidate])

    def MinHashing(self, low_memory=True):
        
        """
        Parameters
        ----------
        low_memory: ndarray, default: True 
                    When low_memory is False, the creation of signature matrix is considerably faster, but at the same 
                    time memory expensive since a hash matrix of shape (np.max(items), n_hash_functions) needs to be
                    stored. If this hash matrix fits into memory, then use low_memory=True.
                    
        Attributes
        ----------
        signature_matrix: ndarray
                          Signature matrix of shape (n_hash_functions, np.max(objects))            
        """
        
        if self._verbose > 0:
            print('\nMinHashing...')
         
        ii32 = np.iinfo(np.int32)
        if self._n_properties < ii32.max:
            self.signature_matrix = np.zeros(shape=(self.n_hash_functions, self._n_items), dtype=np.int32)
            dtype=np.int32
        else:    
            self.signature_matrix = np.zeros(shape=(self.n_hash_functions, self._n_items), dtype=np.int64)
            dtype=np.int64
        
        if low_memory == True:
            hash_values = np.arange(self._n_properties)
            rs = np.random.RandomState(self._random_state)
            for ind in range(self.n_hash_functions):
                rs.shuffle(hash_values)  
                for item in range(len(self._item_inds) - 1):
                    item_properties = self.data[self._item_inds[item]:self._item_inds[item + 1], 1]
                    if len(item_properties) > 0: 
                        item_hash_values = hash_values[item_properties]
                        self.signature_matrix[ind, item] = np.min(item_hash_values)
        else:
            hash_values = np.arange(self._n_properties)
            hash_matrix = np.zeros(shape=(self._n_properties, self.n_hash_functions), dtype=dtype)
            rs = np.random.RandomState(self._random_state)
            for ind in range(self.n_hash_functions):
                rs.shuffle(hash_values)  
                hash_matrix[:, ind] = hash_values
        
            # Fill the signature matrix
            for item in range(len(self._item_inds) - 1):
                item_properties = self.data[self._item_inds[item]:self._item_inds[item + 1], 1] 
                if len(item_properties) > 0:
                    item_hash_matrix = hash_matrix[item_properties, :]  
                    min_hashes = np.amin(item_hash_matrix, axis=0)  
                    self.signature_matrix[:, item] = min_hashes
        if self._verbose > 0:
            print('Signature Matrix Created. MinHashing Completed.\n')                                                            

    def LSH(self, signature_matrix=None):        
        
        """
        Parameters
        ----------
        singature_matrix: ndarray, default: None
                          A signature matrix produced by the MinHashing() method. If None, then the LSH()
                          method uses the signature_matrix produced by the last call to the MinHashing() 
                          method. If a signature_matrix is specified, then this matrix is used by the LSH()
                          method, and any matrix produced by previous calls to MinHashing is ignored.
                    
        Attributes
        ----------
        similar_items: ndarray
                       Array of two columns containing the similar pairs.
                     
        jaccard_similarity: ndarray
                       Array containing the similarities that correspond to each row of similar_items.                
        """
        
        if signature_matrix is not None:
            self.signature_matrix = signature_matrix

        if self._verbose > 0:
            print('Locality Sensitive Hashing...')
        self.similar_items = []
        self._similar_items_dict = dict.fromkeys(self._unique_items, [])
        self.jaccard_similarity_score = []
        for band in range(self.n_bands):
            if self._verbose > 0:
                stdout.write("\rProcessing Band: " + str(band + 1) + "/" + str(self.n_bands))
            buckets = self._hash_to_buckets(band)
            
            for idx, bucketID in enumerate(buckets):
                if len(buckets[bucketID]) > 1:
                    self._compare_candidates(buckets, bucketID)
        
        self.similar_items = np.asarray(self.similar_items)
        self.jaccard_similarity_score = np.asarray(self.jaccard_similarity_score)
        
        if self._verbose > 0:
            print('\nLocality Sensitive Hashing Completed.',
                  '\n' + str(len(self.similar_items)) + ' Similar Pairs Found.',
                  sep='\n')
        
