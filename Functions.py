import numpy as np
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import choice

def BindingStrength(H, S, kT,Kv,K, C, Ea):
	M = len(H)
	N = len(S)
	Pi= np.zeros([M])
	Ph= np.zeros([M])
	E = np.zeros([M,N])
	for i in range(M):
		temp1 = np.zeros(N)
		temp2 = np.zeros(N)
		for j in range(N):
			E[i,j] = np.sum([H[i,k] * S[j,k] for k in range(Kv)]) + np.sum(H[i,k] for k in range(Kv+1,K))
			temp1[j] = C[j] * np.exp((E[i,j] - Ea)/kT)
			temp2[j] = np.exp((E[i,j])/kT)
		Pi[i] = np.sum(temp1)/(1+np.sum(temp1))
		Ph[i] = np.sum(temp2)
	Ctot = np.sum(C)
	Ph = [Ph[i] / (Ph[i] +  1.0/(len(Ph) * Ctot) * (np.sum(Ph[0:(i - 1)]) + np.sum(Ph[(i + 1):len(Ph)]))) for i in
		  range(len(Ph))]
	return Ph, Pi

# def BindingStrength1d(j, H, S, kT,Kv,K, C, Ea):
# 	S = S[j]
# 	C = C[j]
# 	M = len(H)
# 	Pi= np.zeros([M])
# 	Ph= np.zeros([M])
# 	E = np.zeros([M])
# 	for i in range(M):
# 		temp1 = np.zeros(N)
# 		temp2 = np.zeros(N)
# 		E[i] = np.sum([H[i,k] * S[k] for k in range(Kv)]) + np.sum(H[i,k] for k in range(Kv+1,K))
# 		temp = C* np.exp((E[i] - Ea)/kT)
# 		Pi[i] = temp/(1+temp)
# 		Ph[i] = np.exp((E[i])/kT)
# 	Ph = [Ph[i] / (Ph[i] +  1.0/(len(Ph) * C) * (np.sum(Ph[0:(i - 1)]) + np.sum(Ph[(i + 1):len(Ph)]))) for i in
# 		  range(len(Ph))]
# 	return Ph, Pi


def Replicate(H,values,probabilities):
	H = np.array(H)
	H = H[np.random.choice(H.shape[0], int(H.shape[0]*.9), replace=False), :]
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
						#if index < Kv:
							#s_k = 1
						item[index] += s_k * dE
						temp.append(item)
					elif cell_fate < 0.7:
						temp.append(item.copy())
				else:
					temp.append(item.copy())
		H = np.array(temp)
	return H

def generate_H1(S, kT,Kv,K, C, Ea,affinity,values,probabilities):
	H = np.zeros([3,K])
	N = len(S)
	for i in range(3):
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
		Ph,Pi = BindingStrength(H, S, kT,Kv,K, C, Ea)
		# Selection
		# for each cell, the probability of survival is Ph * Pi
		survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
		H = [H[i,:] for i in range(len(survival)) if survival[i]==1]
		# Replication + Mutation
		H=Replicate(H,values,probabilities)
		M = len(H)
		print(len(H))
		# Termination
		if(M>1536):
			break
		if(M<2):
			break
	return(H)


def generate_H2(S, kT,Kv,K, C, Ea,affinity,values,probabilities):
	H = np.zeros([3,K])
	N = len(S)
	for i in range(3):
		for j in range(N):
			alpha = affinity[i]
			h = uniform(-.18, 0.9, size = K)
			h_sum = np.sum(h)
			h = [x/h_sum * alpha for x in h]
			H[i,:] = h
	for i in range(9):
		H = np.concatenate([H,H],0)
	M=len(H)
	tmax = 80
	for t in range(tmax):
		# Calculating Binding Strengh
		Ph,Pi = BindingStrength(H, np.array([S[0]]*3), kT,Kv,K, C, Ea)
		# Selection
		# for each cell, the probability of survival is Ph * Pi
		survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
		H = [H[i,:] for i in range(len(survival)) if survival[i]==1]
		# Replication + Mutation
		H=Replicate(H,values,probabilities)
		M = len(H)
		print(len(H))
		# Termination
		if(M>1536):
			break
		if(M<2):
			break
	tmax = 160
	for t in range(tmax):
		# Calculating Binding Strengh
		Ph,Pi = BindingStrength(H, np.array([S[1],S[2],S[2]]), kT,Kv,K, C, Ea)
		# Selection
		# for each cell, the probability of survival is Ph * Pi
		survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
		H = [H[i,:] for i in range(len(survival)) if survival[i]==1]
		# Replication + Mutation
		H=Replicate(H,values,probabilities)
		M = len(H)
		print(len(H))
		# Termination
		if(M>1536):
			break
		if(M<2):
			break
	return(H)



def generate_H3(S, kT,Kv,K, C, Ea,affinity,values,probabilities):
	H = np.zeros([3,K])
	N = len(S)
	for i in range(3):
		for j in range(N):
			alpha = affinity[i]
			h = uniform(-.18, 0.9, size = K)
			h_sum = np.sum(h)
			h = [x/h_sum * alpha for x in h]
			H[i,:] = h
	for i in range(9):
		H = np.concatenate([H,H],0)
	M=len(H)
	tmax = 80
	for t in range(tmax):
		# Calculating Binding Strengh
		Ph,Pi = BindingStrength(H, np.array([S[0]]*3), kT,Kv,K, C, Ea)
		# Selection
		# for each cell, the probability of survival is Ph * Pi
		survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
		H = [H[i,:] for i in range(len(survival)) if survival[i]==1]
		# Replication + Mutation
		H=Replicate(H,values,probabilities)
		M = len(H)
		print(len(H))
		# Termination
		if(M>1536):
			break
		if(M<2):
			break
	tmax = 80
	for t in range(tmax):
		# Calculating Binding Strengh
		Ph,Pi = BindingStrength(H, np.array([S[1]]*3), kT,Kv,K, C, Ea)
		# Selection
		# for each cell, the probability of survival is Ph * Pi
		survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
		H = [H[i,:] for i in range(len(survival)) if survival[i]==1]
		# Replication + Mutation
		H=Replicate(H,values,probabilities)
		M = len(H)
		print(len(H))
		# Termination
		if(M>1536):
			break
		if(M<2):
			break
	tmax = 80
	for t in range(tmax):
		# Calculating Binding Strengh
		Ph,Pi = BindingStrength(H, np.array([S[2]]*3), kT,Kv,K, C, Ea)
		# Selection
		# for each cell, the probability of survival is Ph * Pi
		survival = [binomial(n=1,p=Ph[i]*Pi[i]) for i in range(len(Ph))]
		H = [H[i,:] for i in range(len(survival)) if survival[i]==1]
		# Replication + Mutation
		H=Replicate(H,values,probabilities)
		M = len(H)
		print(len(H))
		# Termination
		if(M>1536):
			break
		if(M<2):
			break
	return(H)

