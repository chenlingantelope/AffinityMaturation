from numpy.random import normal
from numpy.random import exponential
import numpy as np

# Undefined parameters
Ea = 10.8  #kcal per mole
N = 3 # number of strains
K = 20 # number of residues on the strain (and possible binding sites)
Kv = 18 # number of variable residues
kT = 0.0019872041 * 310  #kcal per mole 
# C = 500*1000*1000 /(6.022140857*10**23) # number of moles per liter, https://academic.oup.com/jid/article/181/4/1280/852328
C = 0.5
C = [C]*N
# Germinal Center 
# Three cells that passes the threshold

affinity = normal(loc=0, scale=(Ea/1.644853626951/K), size=30)
affinity = [x for x in affinity if x>(Ea/K)]

# replicating
S = np.ones([N,K])

for i in range(len(affinity)): 
	alpha = affinity[i]
	h = exponential(size = M)
	h_mean = np.mean(h)
	h = [x/h_mean * alpha for x in h]
	H[i] = h


for i in range(9):
	H = H + H

M=len(H)
E = np.zeros([M, N] )
H = np.zeros([M,N,K])


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
# probability of internalizing 
# mutation

