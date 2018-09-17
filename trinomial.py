import numpy as np
import math as m

class american_pricer:

    def __init__(self,S=1,T=60,K=1,r=0,q=0,sigma=.1, call = True):
        self.S = S
        self.T = T
        self.K = K
        self.r = r
        self.q = q
        self.sigma = sigma
        self.call  = call


    def price(self,N):
        dt = float(self.T/N)
        
        #mu is r-q - (sigma^2)/2
        mu = self.r-self.q-(self.sigma**2)/2.0
        #set sigma max for stability requirements
        smax = 2 * abs(mu) * dt**.5
        smax = max(smax, self.sigma * (2**.5))
        if smax ==0:
            return -9999
        #set up arrays to keep track of steps
        #dimension M
        M = int(5 * (N**.5))
        C_ = np.empty(2*M+1, dtype=np.float64)
        pC_ = np.empty(2*M+1, dtype=np.float64)
        S_ = np.empty(2*M+1, dtype=np.float64,)
        #set probs up, down, and same
        p = float(0.5 * (self.sigma**2) )/ (smax **2)
        p_u = p + 0.5 * mu * dt**.5 / float(smax)
        p_m = 1 - 2 * p
        p_d = p - 0.5 * mu * dt**.5 / float(smax)
        #init payoff
        D = 1.0 / (1 + self.r * dt)
        E = m.exp(smax * dt**.5)
        
        for j in range(0,len(S_)):
            if j ==0:
                S_[j] = self.S * m.exp(-M * smax * dt**.5)
            else:
                S_[j] = S_[j - 1] * E
            if self.call ==True:
                C_[j] = max(S_[j] - self.K, 0)
            else:
                C_[j] = max(self.K-S_[j], 0)

        for k in range(0,N):
            for j in range(1,2 * M):
                pC_[j] = (p_u * C_[j + 1] + p_m * C_[j] + p_d * C_[j - 1])*D
            #set boundaries
            pC_[0] = 2 * pC_[1] - pC_[2]
            pC_[2 * M] = 2 * pC_[2 * M -1] - pC_[2 * M - 2]
            
            for n in range(0,2 * M+1):
                if self.call ==True:
                    C_[n] = max(pC_[n],max(S_[n]-self.K,0))
                else:
                    C_[n] = max(pC_[n],max(self.K-S_[n],0))
        ret = C_[M]
        return ret