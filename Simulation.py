from numpy.random import uniform
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
for i in [1,2]:
	mut = choice(range(Kv),size=5,replace=False)
	for j in mut:
		S[i,j] = -1

# replicating


H=generate_H1(S, kT,Kv,K, C, Ea,affinity,values)

BindingStrength(H, S[j], kT,Kv,K, C, Ea)