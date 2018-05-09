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
delta_E = {(-6.4, -6): 0.025, (-6, -5.6): 0.005, (-5.6, -5.2): 0.035,
(-5.2, -4.8): 0.06, (-4.8, -4.4): 0.06, (-4.4, -4): 0.15,  (-4, -3.6): 0.15, 
(-3.6, -3.2): 0.19, (-3.2, -2.8): 0.23, (-2.8, -2.4): 0.23, (-2.4, -2): 0.28, 
(-2, -1.6): .32, (-1.6, -1.2): .32, (-1.2, -.8): .55, (-.8, -.4): 0.42, (-0.4, 0): 0.7, (0, 0.4): 0.67, 
(0.4, 0.8): 0.35, (0.8, 1.2): .13, (1.2, 1.6): 0.08}
total = sum(delta_E.itervalues(), 0.0)
delta_E = {k: v / total for k, v in delta_E.items()}

values, probabilities = zip(*delta_E.items())

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
		h = uniform(-.18, 0.9, size = K)
		h_sum = np.sum(h)
		h = [x/h_sum * alpha for x in h]
		H[i,:] = h


for i in range(9):
	H = np.concatenate([H,H],0)

M=len(H)

tmax = 240
for t in range(tmax):
	# Calculating Binding Strengh
	Ph,Pi = BindingStrength(M, N, H, S, kT,Kv,K, C, Ea)
	# Selection 
	# for each cell, the probability of survival is Ph * Pi
	survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
	temp = [H[i,:] for i in range(len(survival)) if survival[i]==1]
	#temp = np.rollaxis(np.dstack(temp),-1)
	H = temp
	print len(H)
	# Replication + Mutation
	for __ in range(0,2):
		temp = []
		for item in H:
			for _ in range(0,2):
				# Check if there is a mutation
				if uniform(0, 1) < 0.14:
					# Type of mutation
					cell_fate = uniform(0,1)
					if cell_fate < 0.2:
						uniform_sampling_range = values[np.random.choice(range(0,len(probabilities)), p=probabilities)]
						dE = uniform(uniform_sampling_range[0], uniform_sampling_range[1])
						item = item.copy()
						index = np.random.choice(range(0, len(item)))
						s_k = 1
						if index < Kv:
							s_k = -1
						item[index] += s_k * dE

						temp.append(item)
					elif cell_fate < 0.7:
						temp.append(item.copy())

				else:
					temp.append(item.copy())
		H = np.array(temp)
	#H = np.concatenate([H,H,H,H],0)
	M = len(H)

	# Termination
	if(M>1536):
		break 
	if(M<2):
		break

