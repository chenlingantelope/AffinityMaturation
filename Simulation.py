from numpy.random import uniform
from numpy.random import exponential
from numpy.random import binomial
from numpy.random import choice
import numpy as np
from Functions import *
# Undefined parameters
Ea = 10.8  #kcal per mole
N = 3 # number of strains
K = 46 # number of residues on the strain (and possible binding sites)
Kv = 28 # number of variable residues
Ks = 6
kT = 0.0019872041 * 310  #kcal per mole 
# C = 500*1000*1000 /(6.022140857*10**23) # number of moles per liter, https://academic.oup.com/jid/article/181/4/1280/852328
C = 0.5
C = [C]*N
# Germinal Center 
# Three cells that passes the threshold

affinity = uniform(low=Ea-kT,high=Ea+kT,size=3)

S = np.ones([N,K])
# add mutation in Sk
for i in range(len(S)):
	mut = choice(range(Kv),size=3,replace=False)
	for j in mut:
		S[i,j] = -1

# replicating
H = np.zeros([3,K])

for i in range(len(affinity)): 
	for j in range(N):
		alpha = affinity[i]
		h = exponential(size = K)
		h_sum = np.sum(h)
		h = [x/h_sum * alpha for x in h]
		H[i,:] = h


for i in range(9):
	H = np.concatenate([H,H],0)

M=len(H)

tmax = 240
for t in range(tmax):
	# Calculating Binding Strengh
	Ph,Pi = BindingStrength(M, N, H, S, kT,Kv,K)
	# Selection 
	# for each cell, the probability of survival is Ph * Pi
	survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
	temp = [H[i,:] for i in range(len(survival)) if survival[i]==1]
	temp = np.rollaxis(np.dstack(temp),-1)
	H = temp
	# Replication + Mutation
	H = np.concatenate([H,H,H,H],0)
	M = len(H)
	# Termination 
	if(M>1536):
		break 
	if(M<2):
		break

