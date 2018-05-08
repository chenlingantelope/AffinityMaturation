from numpy.random import normal
from numpy.random import exponential
import numpy as np

# Undefined parameters
Ea = 0.1
M = 20
N = 18
# Germinal Center 
# Three cells that passes the threshold

affinity = normal(loc=0, scale=Ea/1.644853626951, size=30)
affinity = [x for x in affinity if x>Ea ]

# replicating

s = [1]*M

H = [[]] * len(affinity)
for i in range(len(affinity)): 
	alpha = affinity[i]
	h = exponential(size = 10)
	h_mean = np.mean(h)
	h = [x/h_mean * alpha for x in h]
	H[i] = h


for i in range(9):
	affinity = affinity + affinity
