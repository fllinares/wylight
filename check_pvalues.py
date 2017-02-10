import sys
import numpy as np
from scipy.stats import hypergeom

filename_pvals = sys.argv[1]
filename_labels = sys.argv[2]
filename_output = sys.argv[3]

out = np.loadtxt(filename_pvals, delimiter=',', skiprows=1)  # Load pvalues and margins a, x
N_out = out.shape[0]  # Number of significant itemsets found
y = np.loadtxt(filename_labels)  # Load class labels
n = y.sum(); N = len(y) # Compute margins n, N of 2x2 table
if n > N/2:  # Ensure n is the number of samples of the minority class
    n = N-n
pvals = np.zeros((N_out,))
for i in range(N_out):
    a = out[i,0]; x = out[i,1]  # Retrieve from output file margins a, x for i-th itemset
    hypergeom_pmf = hypergeom.pmf(np.arange(N+1),N,n,x)  # Get PDF of Hypergeom RV of proper margins
    pvals[i] = hypergeom_pmf[hypergeom_pmf <= hypergeom_pmf[a]].sum()  # Definition of two-tail p-val
np.savetxt(filename_output,np.column_stack((out,pvals)),fmt='%d,%d,%e,%e')
