import numpy as np
import JaccardSimilaritySearch as jss

data = np.load('user_movie20M.npy')

#signature_matrix = np.load('signature_matrix.npy')

mhlsh = jss.MinHashingLSH(data, threshold=0.5,min_set_size=10, random_state=10, verbose=1)
mhlsh.MinHashing(low_memory=False)
mhlsh.LSH()
