from numpy.random import uniform
from numpy.random import exponential
import numpy as np

# Undefined parameters
Ea = 10.8  #kcal per mole
N = 3 # number of strains
K = 46 # number of residues on the strain (and possible binding sites)
Kv = 40 # number of variable residues
kT = 0.0019872041 * 310  #kcal per mole 
# C = 500*1000*1000 /(6.022140857*10**23) # number of moles per liter, https://academic.oup.com/jid/article/181/4/1280/852328
C = 0.5
C = [C]*N
# Germinal Center 
# Three cells that passes the threshold

affinity = uniform(low=Ea-kT,high=Ea+kT,size=3)

# replicating
S = np.ones([N,K])
# add mutation in Sk
M=len(H)
H = np.zeros([3,N,K])

for i in range(len(affinity)): 
	for j in range(N):
		alpha = affinity[i]
		h = exponential(size = K)
		h_mean = np.mean(h)
		h = [x/h_mean * alpha for x in h]
		H[i,j] = h

for i in range(9):
	H = H + H

E = np.zeros([M, N] )

# Calculating Binding Strengh
Pi= np.zeros([M])
Ph= np.zeros([M])
for i in range(M):
	temp1 = np.zeros(N)
	temp2 = np.zeros(N)
	for j in range(N):
		E[i,j] = np.sum([H[i,j,k] * S[k,j] for k in range(Kv)]) + np.sum(H[i][k] for k in range(Kv+1,K))
		temp1[j] = C[j] * np.exp((E[i,j] - Ea)/kT)
		temp2[j] = np.exp((E[i,j])/kT)
	Pi[i] = np.sum(temp1)/(1+np.sum(temp1))
	Ph[i] = np.sum(temp2)


# Selection 


# Replication + Mutation

# Termination 
