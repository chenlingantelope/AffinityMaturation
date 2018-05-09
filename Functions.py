import numpy as np
def BindingStrength(M, N, H, S, kT,Kv,K, C, Ea):
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
	#Ph = [ Ph[i]/(Ph[i] + Ctot/len(Ph)*(np.sum(Ph[0:(i-1)])+np.sum(Ph[(i+1):len(Ph)])))for i in range(len(Ph))]
	Ph = [Ph[i] / (Ph[i] +  1.0/(len(Ph) * Ctot) * (np.sum(Ph[0:(i - 1)]) + np.sum(Ph[(i + 1):len(Ph)]))) for i in
		  range(len(Ph))]
	return Ph, Pi

