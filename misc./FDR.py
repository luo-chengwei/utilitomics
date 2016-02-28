import sys

# input a pvalue list and alpha level, return a qvalue list
def bh(pValues):
	# append an index identifier to each p-value
    indexed_pval = []
    for i, p in zip(range(len(pValues)), pValues):
    	indexed_pval.append((i, p))
    # sort p-values in descending order
    sorted_pval = sorted(indexed_pval, key = lambda x: x[1], reverse = True)
    
    # determine significant features
    N = len(pValues)
    modifier = N
    self.numSignFeatures = 0
    qvals = []
    for i in range(N): qvals.append(0.)
    for r, (i, p) in enumerate(sorted_pval):
    	q = min(1., float(p) * N/(N-r))
    	qvals[r] = q
    return qvals