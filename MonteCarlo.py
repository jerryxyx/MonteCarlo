import numpy as np
import scipy as sp
import scipy.stats
from cvxopt import matrix, solvers
import pandas as pd


class MonteCarlo:
    def __init__(self,S0,K,T,r,q,sigma,kappa=0,theta=0,xi=0,rho=0,V0=0,underlying_process="geometric brownian motion"):
        self.underlying_process = underlying_process
        self.S0 = S0
        self.K = K
        self.T = T
        self.r = r
        self.q = q
        self.sigma = sigma
        self.kappa = kappa
        self.theta = theta
        self.rho = rho
        self.V0 = V0
        self.xi = xi
        
        self.value_results = None
    
    # view antithetic variates as a option of simulation method to reduce the variance    
    def simulate(self, n_trails, n_steps, antitheticVariates=False, boundaryScheme="Higham and Mao"):
        
        dt = self.T/n_steps
        mu = self.r - self.q
        self.n_trails = n_trails
        self.n_steps = n_steps
        self.boundaryScheme = boundaryScheme
        
        if(self.underlying_process=="geometric brownian motion"):
#             first_step_prices = np.ones((n_trails,1))*np.log(self.S0)
            log_price_matrix = np.zeros((n_trails,n_steps))
            normal_matrix = np.random.normal(size=(n_trails,n_steps))
            if(antitheticVariates==True):
                n_trails *= 2
                self.n_trails = n_trails
                normal_matrix = np.concatenate((normal_matrix,-normal_matrix),axis=0)
            cumsum_normal_matrix = normal_matrix.cumsum(axis=1)
#             log_price_matrix = np.concatenate((first_step_prices,log_price_matrix),axis=1)
            deviation_matrix = cumsum_normal_matrix*self.sigma*np.sqrt(dt) + \
    (mu-self.sigma**2/2)*dt*np.arange(1,n_steps+1)
            log_price_matrix = deviation_matrix+np.log(self.S0)
            price_matrix = np.exp(log_price_matrix)
            price_zero = (np.ones(n_trails)*self.S0)[:,np.newaxis]
            price_matrix = np.concatenate((price_zero,price_matrix),axis=1)
            self.price_matrix = price_matrix
        
        elif(self.underlying_process=="CIR model"):
            # generate correlated random variables
            randn_matrix_v = np.random.normal(size=(n_trails,n_steps))
            if(antitheticVariates==True):
                n_trails *= 2
                self.n_trails = n_trails
                randn_matrix_v = np.concatenate(( randn_matrix_v, -randn_matrix_v),axis=0)

            # boundary scheme fuctions
            if(boundaryScheme=="absorption"):
                f1=f2=f3=lambda x: np.maximum(x,0)
            elif(boundaryScheme=="reflection"):
                f1=f2=f3=np.absolute
            elif(boundaryScheme=="Higham and Mao"):
                f1=f2=lambda x: x
                f3 = np.absolute
            elif(boundaryScheme=="partial truncation"):
                f1=f2=lambda x: x
                f3=lambda x: np.maximum(x,0)
            elif(boundaryScheme=="full truncation"):
                f1 = lambda x: x
                f2=f3= lambda x: np.maximum(x,0)
            
            # simulate CIR process
            V_matrix = np.zeros((n_trails,n_steps+1))
            V_matrix[:,0] = self.S0

            for j in range(self.n_steps):
                V_matrix[:,j+1] = f1(V_matrix[:,j]) - self.kappa*dt*(f2(V_matrix[:,j])-self.theta) +\
                    self.xi*np.sqrt(f3(V_matrix[:,j]))*np.sqrt(dt)*randn_matrix_v[:,j]
                V_matrix[:,j+1] = f3(V_matrix[:,j+1])
                
            price_matrix = V_matrix
            self.price_matrix = price_matrix
            
        
        elif(self.underlying_process=="Heston model"):
            # generate correlated random variables
            randn_matrix_1 = np.random.normal(size=(n_trails,n_steps))
            randn_matrix_2 = np.random.normal(size=(n_trails,n_steps))
            randn_matrix_v = randn_matrix_1
            randn_matrix_S = self.rho*randn_matrix_1 + np.sqrt(1-self.rho**2)*randn_matrix_2
            if(antitheticVariates==True):
                n_trails *= 2
                self.n_trails = n_trails
                randn_matrix_v = np.concatenate(( randn_matrix_v, +randn_matrix_v),axis=0)
                randn_matrix_S = np.concatenate(( randn_matrix_S, -randn_matrix_S),axis=0)

            # boundary scheme fuctions
            if(boundaryScheme=="absorption"):
                f1=f2=f3=lambda x: np.maximum(x,0)
            elif(boundaryScheme=="reflection"):
                f1=f2=f3=np.absolute
            elif(boundaryScheme=="Higham and Mao"):
                f1=f2=lambda x: x
                f3 = np.absolute
            elif(boundaryScheme=="partial truncation"):
                f1=f2=lambda x: x
                f3=lambda x: np.maximum(x,0)
            elif(boundaryScheme=="full truncation"):
                f1 = lambda x: x
                f2=f3= lambda x: np.maximum(x,0)
            
            # simulate stochastic volatility process
            V_matrix = np.zeros((n_trails,n_steps+1))
            V_matrix[:,0] = self.V0
            log_price_matrix = np.zeros((n_trails,n_steps+1))
            log_price_matrix[:,0] = np.log(self.S0)
            for j in range(self.n_steps):
#                 V_matrix[:,j+1] = self.kappa*self.theta*dt + (1-self.kappa*dt)*V_matrix[:,j] +\
#                     self.xi*np.sqrt(V_matrix[:,j]*dt)*randn_matrix_v[:,j]
                V_matrix[:,j+1] = f1(V_matrix[:,j]) - self.kappa*dt*(f2(V_matrix[:,j])-self.theta) +\
                    self.xi*np.sqrt(f3(V_matrix[:,j]))*np.sqrt(dt)*randn_matrix_v[:,j]
                V_matrix[:,j+1] = f3(V_matrix[:,j+1])
                log_price_matrix[:,j+1] = log_price_matrix[:,j] + (mu - V_matrix[:,j]/2)*dt +\
                    np.sqrt(V_matrix[:,j]*dt)*randn_matrix_S[:,j]
            price_matrix = np.exp(log_price_matrix)
            self.price_matrix = price_matrix
            
        return price_matrix
    
    def BlackScholesPricer(self,option_type='c'):
        S = self.S0
        K = self.K
        T = self.T
        r = self.r
        q = self.q
        sigma = self.sigma
        d1 = (np.log(S/K)+(r-q)*T +0.5*sigma**2*T)/(sigma*np.sqrt(T))
        d2 = d1 - sigma*np.sqrt(T)
        N = lambda x: sp.stats.norm.cdf(x)
        call = np.exp(-q*T) * S * N(d1) - np.exp(-r*T) * K * N(d2)
        put = call - np.exp(-q*T) * S + K*np.exp(-r*T)
        if(option_type=="c"):
            return call
        elif(option_type=="p"):
            return put
        else:
            print("please enter the option type: (c/p)")
        pass
    
    def MCPricer(self,option_type='c'):
        price_matrix = self.price_matrix
        # k = n_steps
        dt = self.T/self.n_steps
        df = np.exp(- self.r*dt)
        n_basis = len(func_list)
        n_trails = self.n_trails
        n_steps = self.n_steps
        
        if(option_type=="c"):
            payoff = (price_matrix[:,n_steps] - strike)
        elif(option_type=="p"):
            payoff = (strike - price_matrix[:,n_steps])
        else:
            print("please enter the option type: (c/p)")
            return
        
        payoff = matrix(np.where(payoff<0,0,payoff))
#         vk = payoff*df
        value_results = payoff*np.exp(-risk_free_rate*time_to_maturity)
        regular_mc_price = np.average(payoff*np.exp(-risk_free_rate*time_to_maturity))
        self.payoff = payoff
        self.mc_price = regular_mc_price
        self.value_results = value_results
        return regular_mc_price
    
    def BSDeltaHedgedPricer(self,option_type="c"):
        
        regular_mc_price = self.MCPricer(option_type=option_type)
        dt = self.T/self.n_steps
        df = np.exp(- self.r*dt)
        df2 = np.exp(-(self.r-self.q)*dt)
        
        # Delta hedged cash flow
        def Delta_fun(x,tau,option_type):
            d1 = (np.log(x/self.K) + (self.r-self.q)*tau + self.sigma**2*tau/2)/(self.sigma*np.sqrt(tau))
            if(option_type=='c'):
                return sp.stats.norm.cdf(d1)
            elif(option_type=='p'):
                return -sp.stats.norm.cdf(-d1)
        
        discounted_hedge_cash_flow = np.zeros(self.n_trails)
        for i in range(self.n_trails):
            Sk_array = self.price_matrix[i,:]
            bi_diag_matrix = np.diag([-1]*(n_steps),0) + np.diag([df2]*(n_steps-1),1)
            # (Sk+1 exp(-r dt) - Sk) exp(-r*(tk-t0))
            discounted_stock_price_change = np.dot(bi_diag_matrix,Sk_array[:-1])
            discounted_stock_price_change[-1] += Sk_array[-1]*df2
            discounted_stock_price_change *= np.exp(-self.r*np.arange(n_steps)*dt)
            tau_array = dt*np.arange(self.n_steps,0,-1)
            Delta_array = np.array([Delta_fun(Sk,tau,option_type) for Sk,tau in zip(Sk_array[:-1],tau_array)])
            discounted_hedge_cash_flow[i] = np.dot(Delta_array,discounted_stock_price_change)
        
        BSDeltaBased_mc_price = regular_mc_price - discounted_hedge_cash_flow.mean()
#         print("The average discounted hedge cash flow: {}".format(discounted_hedge_cash_flow.mean()))
        
        value_results = self.payoff*np.exp(-self.r*self.T) - discounted_hedge_cash_flow
#         print("Sanity check {} = {}".format(value_results.mean(),BSDeltaBased_mc_price))
        self.value_results = value_results
        
        return BSDeltaBased_mc_price
    
    def OHMCPricer(self,option_type='c', func_list=[lambda x: x**0, lambda x: x]):
        def _calculate_Q_matrix(S_k,S_kp1,df,df2,func_list):
            dS = df2*S_kp1 - S_k
            A = np.array([func(S_k) for func in func_list]).T
            B = (np.array([func(S_k) for func in func_list])*dS).T
            return np.concatenate((-A,B),axis=1)
        
        price_matrix = self.price_matrix
        # k = n_steps
        dt = self.T/self.n_steps
        df = np.exp(- self.r*dt)
        df2 = np.exp(-(self.r-self.q)*dt)
        n_basis = len(func_list)
        n_trails = self.n_trails
        n_steps = self.n_steps
        
        if(option_type=="c"):
            payoff = (price_matrix[:,n_steps] - strike)
        elif(option_type=="p"):
            payoff = (strike - price_matrix[:,n_steps])
        else:
            print("please enter the option type: (c/p)")
            return
        
        payoff = matrix(np.where(payoff<0,0,payoff))
        vk = payoff*df
#         print("regular MC price",regular_mc_price)
    
        # k = 1,...,n_steps-1
        for k in range(n_steps-1,0,-1):
            Sk = price_matrix[:,k]
            Skp1 = price_matrix[:,k+1]
            Qk = matrix(_calculate_Q_matrix(Sk,Skp1,df,df2,func_list))
            P = Qk.T * Qk
            q = Qk.T * vk
            A = matrix(np.ones(n_trails,dtype=np.float64)).T * Qk
            b = - matrix(np.ones(n_trails,dtype=np.float64)).T * vk
            sol = solvers.coneqp(P=P,q=q,A=A,b=b)
            ak = sol["x"][:n_basis]
            bk = sol["x"][n_basis:]
            vk = matrix(np.array([func(price_matrix[:,k]) for func in func_list])).T*ak*df
        
        # k = 0
        v0 = vk
        S0 = price_matrix[:,0]
        S1 = price_matrix[:,1]
        dS0 = df2*S1 - S0
        Q0 = np.concatenate((-np.ones(n_trails)[:,np.newaxis],dS0[:,np.newaxis]),axis=1)
        Q0 = matrix(Q0)
        P = Q0.T*Q0
        q = Q0.T*v0
        A = matrix(np.ones(n_trails,dtype=np.float64)).T * Q0
        b = - matrix(np.ones(n_trails,dtype=np.float64)).T * v0
        C1 = matrix(ak).T * np.array([func(S1) for func in func_list]).T
        sol = solvers.coneqp(P=P,q=q,A=A,b=b)
        self.sol = sol
        residual_risk = (v0.T*v0 + 2*sol["primal objective"])/n_trails
        self.residual_risk = residual_risk[0]    # the value of unit matrix
        
        return sol["x"][0]
    
    def standard_error(self):
        # can not apply to the OHMC since its result is not obtained by averaging
        # sample variance
        sample_var = np.var(self.value_results,ddof=1)
        std_estimate = np.sqrt(sample_var)
        standard_err = std_estimate/np.sqrt(n_trails)
        return standard_err
        
    def pricing(self, option_type='c', func_list=[lambda x: x**0, lambda x: x]):
        OHMC_price = self.OHMCPricer(option_type=option_type,func_list=func_list)
        regular_mc_price = self.MCPricer(option_type=option_type)
        black_sholes_price = self.BlackScholesPricer(option_type)
        return({"OHMC": OHMC_price,"regular MC": regular_mc_price,"Black-Scholes":black_sholes_price})
    
    def hedging(self):
        S = self.S0
        K = self.K
        T = self.T
        r = self.r
        q = self.q
        sigma = self.sigma
        d1 = (np.log(S/K)+(r-q)*T +0.5*sigma**2*T)/(sigma*np.sqrt(T))
        d2 = d1 - sigma*np.sqrt(T)
        N = lambda x: sp.stats.norm.cdf(x)
        return({"OHMC optimal hedge": -self.sol["x"][1],"Black-Scholes delta hedge":N(d1),"OHMC residual risk":self.residual_risk})
        