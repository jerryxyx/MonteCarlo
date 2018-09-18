import numpy as np
import scipy as sp
import scipy.stats
from cvxopt import matrix, solvers
from sklearn.linear_model import LinearRegression


class MonteCarlo:
    def __init__(self, S0, K, T, r, q, sigma, kappa=0, theta=0, xi=0, rho=0, V0=0,
                 underlying_process="geometric brownian motion"):
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
    def simulate(self, n_trials, n_steps, antitheticVariates=False, boundaryScheme="Higham and Mao"):

        dt = self.T / n_steps
        mu = self.r - self.q
        self.n_trials = n_trials
        self.n_steps = n_steps
        self.boundaryScheme = boundaryScheme

        if (self.underlying_process == "geometric brownian motion"):
            #             first_step_prices = np.ones((n_trials,1))*np.log(self.S0)
            log_price_matrix = np.zeros((n_trials, n_steps))
            normal_matrix = np.random.normal(size=(n_trials, n_steps))
            if (antitheticVariates == True):
                n_trials *= 2
                self.n_trials = n_trials
                normal_matrix = np.concatenate((normal_matrix, -normal_matrix), axis=0)
            cumsum_normal_matrix = normal_matrix.cumsum(axis=1)
            #             log_price_matrix = np.concatenate((first_step_prices,log_price_matrix),axis=1)
            deviation_matrix = cumsum_normal_matrix * self.sigma * np.sqrt(dt) + \
                               (mu - self.sigma ** 2 / 2) * dt * np.arange(1, n_steps + 1)
            log_price_matrix = deviation_matrix + np.log(self.S0)
            price_matrix = np.exp(log_price_matrix)
            price_zero = (np.ones(n_trials) * self.S0)[:, np.newaxis]
            price_matrix = np.concatenate((price_zero, price_matrix), axis=1)
            self.price_matrix = price_matrix

        elif (self.underlying_process == "CIR model"):
            # generate correlated random variables
            randn_matrix_v = np.random.normal(size=(n_trials, n_steps))
            if (antitheticVariates == True):
                n_trials *= 2
                self.n_trials = n_trials
                randn_matrix_v = np.concatenate((randn_matrix_v, -randn_matrix_v), axis=0)

            # boundary scheme fuctions
            if (boundaryScheme == "absorption"):
                f1 = f2 = f3 = lambda x: np.maximum(x, 0)
            elif (boundaryScheme == "reflection"):
                f1 = f2 = f3 = np.absolute
            elif (boundaryScheme == "Higham and Mao"):
                f1 = f2 = lambda x: x
                f3 = np.absolute
            elif (boundaryScheme == "partial truncation"):
                f1 = f2 = lambda x: x
                f3 = lambda x: np.maximum(x, 0)
            elif (boundaryScheme == "full truncation"):
                f1 = lambda x: x
                f2 = f3 = lambda x: np.maximum(x, 0)

            # simulate CIR process
            V_matrix = np.zeros((n_trials, n_steps + 1))
            V_matrix[:, 0] = self.S0

            for j in range(self.n_steps):
                V_matrix[:, j + 1] = f1(V_matrix[:, j]) - self.kappa * dt * (f2(V_matrix[:, j]) - self.theta) + \
                                     self.xi * np.sqrt(f3(V_matrix[:, j])) * np.sqrt(dt) * randn_matrix_v[:, j]
                V_matrix[:, j + 1] = f3(V_matrix[:, j + 1])

            price_matrix = V_matrix
            self.price_matrix = price_matrix


        elif (self.underlying_process == "Heston model"):
            # generate correlated random variables
            randn_matrix_1 = np.random.normal(size=(n_trials, n_steps))
            randn_matrix_2 = np.random.normal(size=(n_trials, n_steps))
            randn_matrix_v = randn_matrix_1
            randn_matrix_S = self.rho * randn_matrix_1 + np.sqrt(1 - self.rho ** 2) * randn_matrix_2
            if (antitheticVariates == True):
                n_trials *= 2
                self.n_trials = n_trials
                randn_matrix_v = np.concatenate((randn_matrix_v, +randn_matrix_v), axis=0)
                randn_matrix_S = np.concatenate((randn_matrix_S, -randn_matrix_S), axis=0)

            # boundary scheme fuctions
            if (boundaryScheme == "absorption"):
                f1 = f2 = f3 = lambda x: np.maximum(x, 0)
            elif (boundaryScheme == "reflection"):
                f1 = f2 = f3 = np.absolute
            elif (boundaryScheme == "Higham and Mao"):
                f1 = f2 = lambda x: x
                f3 = np.absolute
            elif (boundaryScheme == "partial truncation"):
                f1 = f2 = lambda x: x
                f3 = lambda x: np.maximum(x, 0)
            elif (boundaryScheme == "full truncation"):
                f1 = lambda x: x
                f2 = f3 = lambda x: np.maximum(x, 0)

            # simulate stochastic volatility process
            V_matrix = np.zeros((n_trials, n_steps + 1))
            V_matrix[:, 0] = self.V0
            log_price_matrix = np.zeros((n_trials, n_steps + 1))
            log_price_matrix[:, 0] = np.log(self.S0)
            for j in range(self.n_steps):
                #                 V_matrix[:,j+1] = self.kappa*self.theta*dt + (1-self.kappa*dt)*V_matrix[:,j] +\
                #                     self.xi*np.sqrt(V_matrix[:,j]*dt)*randn_matrix_v[:,j]
                V_matrix[:, j + 1] = f1(V_matrix[:, j]) - self.kappa * dt * (f2(V_matrix[:, j]) - self.theta) + \
                                     self.xi * np.sqrt(f3(V_matrix[:, j])) * np.sqrt(dt) * randn_matrix_v[:, j]
                V_matrix[:, j + 1] = f3(V_matrix[:, j + 1])
                log_price_matrix[:, j + 1] = log_price_matrix[:, j] + (mu - V_matrix[:, j] / 2) * dt + \
                                             np.sqrt(V_matrix[:, j] * dt) * randn_matrix_S[:, j]
            price_matrix = np.exp(log_price_matrix)
            self.price_matrix = price_matrix

        return price_matrix

    def LSM(self, option_type="c", func_list=[lambda x: x ** 0, lambda x: x]):
        dt = self.T / self.n_steps
        df = np.exp(-self.r * dt)
        K = self.K
        price_matrix = self.price_matrix
        n_trials = self.n_trials
        n_steps = self.n_steps
        exercise_matrix = np.zeros(price_matrix.shape,dtype=bool)
        # american_values_matrix = np.zeros(price_matrix.shape)

        if (option_type == "c"):
            payoff_fun = lambda x: np.maximum(x - K, 0)
        elif (option_type == "p"):
            payoff_fun = lambda x: np.maximum(K - x, 0)

        # when contract is at the maturity
        stock_prices_t = price_matrix[:, -1]
        exercise_values_t = payoff_fun(stock_prices_t)
        holding_values_t = exercise_values_t
        # american_values_matrix[:,-1] = exercise_values_t

        # before maturaty
        for i in np.arange(n_steps)[:0:-1]:
            holding_values_tp1 = holding_values_t
            exercise_values_tp1 = exercise_values_t
            stock_prices_t = price_matrix[:, i]
            ITM_filter = payoff_fun(stock_prices_t) > 0  # ITM
            A_matrix = np.array([func(stock_prices_t) for func in func_list]).T
            b_matrix = np.maximum(holding_values_tp1,exercise_values_tp1)[:, np.newaxis] * df
            A_prime_matrix = A_matrix[ITM_filter, :]
            b_prime_matrix = b_matrix[ITM_filter, :]
            lr = LinearRegression(fit_intercept=False)
            lr.fit(A_prime_matrix, b_prime_matrix)
            holding_values_t = np.dot(A_matrix, lr.coef_.T)[:, 0]
            exercise_values_t = payoff_fun(stock_prices_t)
            # american_values_matrix[:,i] = np.maximum(holding_values_t,exercise_values_t)
            exercise_filter = (exercise_values_t > holding_values_t) & ITM_filter
            exercise_matrix[exercise_filter, i] = 1

        # i=0
        # holding_values_tp1 = holding_values_t
        #
        # stock_price_0 = price_matrix[0, 0]
        # holding_value_0 = np.mean(holding_values_tp1 * df)
        # exercise_value_0 = payoff_fun(stock_price_0)
        # exercise_matrix[:, 0] = exercise_value_0 > holding_value_0

        # redefine the exercise matrix since we should only exercise once
        holding_matrix = np.zeros(exercise_matrix.shape, dtype=bool)
        for i in np.arange(n_trials):
            exercise_row = exercise_matrix[i, :]
            if (exercise_row.any()):
                exercise_idx = np.where(exercise_row == 1)[0][0]
                exercise_row[exercise_idx + 1:] = 0
                holding_matrix[i,:exercise_idx+1] = 1
            else:
                exercise_row[-1] = 1
                holding_matrix[i,:] = 1

        self.holding_matrix = holding_matrix
        self.exercise_matrix = exercise_matrix
        return exercise_matrix
    
    def LSM2(self, option_type="c", func_list=[lambda x: x ** 0, lambda x: x]):
        dt = self.T / self.n_steps
        df = np.exp(-self.r * dt)
        df2 = np.exp(-(self.r - self.q) * dt)
        K = self.K
        price_matrix = self.price_matrix
        n_trials = self.n_trials
        n_steps = self.n_steps
        exercise_matrix = np.zeros(price_matrix.shape,dtype=bool)
        american_values_matrix = np.zeros(price_matrix.shape)
        # american_values_matrix = np.zeros(price_matrix.shape)

        if (option_type == "c"):
            payoff_fun = lambda x: np.maximum(x - K, 0)
        elif (option_type == "p"):
            payoff_fun = lambda x: np.maximum(K - x, 0)

        # when contract is at the maturity
        stock_prices_t = price_matrix[:, -1]
        exercise_values_t = payoff_fun(stock_prices_t)
        holding_values_t = exercise_values_t
        american_values_matrix[:,-1] = exercise_values_t

        # before maturaty
        for i in np.arange(n_steps)[:0:-1]:
            # A1
            holding_values_tp1 = holding_values_t
            exercise_values_tp1 = exercise_values_t
            stock_prices_t = price_matrix[:, i]
            ITM_filter = payoff_fun(stock_prices_t) > 0  # ITM
            A_matrix = np.array([func(stock_prices_t) for func in func_list]).T
            b_matrix = np.maximum(holding_values_tp1,exercise_values_tp1)[:, np.newaxis] * df
            A_prime_matrix = A_matrix[ITM_filter, :]
            b_prime_matrix = b_matrix[ITM_filter, :]
            lr = LinearRegression(fit_intercept=False)
            lr.fit(A_prime_matrix, b_prime_matrix)
            holding_values_t = np.dot(A_matrix, lr.coef_.T)[:, 0]
            exercise_values_t = payoff_fun(stock_prices_t)
            exercise_filter = (exercise_values_t > holding_values_t) & ITM_filter
            exercise_matrix[exercise_filter, i] = 1

            # A2: there is no way to construct a portfolio
            american_values_tp1 = american_values_matrix[:,i+1]
            lr.fit(A_matrix, american_values_tp1[:,np.newaxis]*df)
            holding_values_t2 = np.dot(A_matrix, lr.coef_.T)[:, 0]
            american_values_t = np.maximum(holding_values_t2,exercise_values_t)
            american_values_matrix[:,i] = american_values_t
            
        
        # i=0
        # regular martingale pricing
        american_value1 = american_values_matrix[:,1].mean() * df
        # with control variates
        v0 = matrix((american_values_matrix[:,1] * df)[:,np.newaxis])
        S0 = price_matrix[:, 0]
        S1 = price_matrix[:, 1]
        dS0 = df * S1 - S0
        Q0 = np.concatenate((-np.ones(n_trials)[:, np.newaxis], dS0[:, np.newaxis]), axis=1)
        Q0 = matrix(Q0)
        P = Q0.T * Q0
        q = Q0.T * v0
        # print(q)
        A = matrix(np.ones(n_trials, dtype=np.float64)).T * Q0
        b = - matrix(np.ones(n_trials, dtype=np.float64)).T * v0
        # C1 = matrix(ak).T * np.array([func(S1) for func in func_list]).T
        sol = solvers.coneqp(P=P, q=q, A=A, b=b)
        self.sol = sol
        residual_risk = (v0.T * v0 + 2 * sol["primal objective"]) / n_trials
        self.residual_risk = residual_risk[0]  # the value of unit matrix
        american_value2 = sol["x"][0]
        delta_hedge = sol["x"][1]
        print(american_value1,american_value2,delta_hedge)
        american_values_matrix[:,0] = american_value2
        
        # holding_values_tp1 = holding_values_t
        #
        # stock_price_0 = price_matrix[0, 0]
        # holding_value_0 = np.mean(holding_values_tp1 * df)
        # exercise_value_0 = payoff_fun(stock_price_0)
        # exercise_matrix[:, 0] = exercise_value_0 > holding_value_0

        # redefine the exercise matrix since we should only exercise once
        holding_matrix = np.zeros(exercise_matrix.shape, dtype=bool)
        for i in np.arange(n_trials):
            exercise_row = exercise_matrix[i, :]
            if (exercise_row.any()):
                exercise_idx = np.where(exercise_row == 1)[0][0]
                exercise_row[exercise_idx + 1:] = 0
                holding_matrix[i,:exercise_idx+1] = 1
            else:
                exercise_row[-1] = 1
                holding_matrix[i,:] = 1

        self.holding_matrix = holding_matrix
        self.exercise_matrix = exercise_matrix
        self.american_values_matrix = american_values_matrix
        return exercise_matrix

    def BlackScholesPricer(self, option_type='c'):
        S = self.S0
        K = self.K
        T = self.T
        r = self.r
        q = self.q
        sigma = self.sigma
        d1 = (np.log(S / K) + (r - q) * T + 0.5 * sigma ** 2 * T) / (sigma * np.sqrt(T))
        d2 = d1 - sigma * np.sqrt(T)
        N = lambda x: sp.stats.norm.cdf(x)
        call = np.exp(-q * T) * S * N(d1) - np.exp(-r * T) * K * N(d2)
        put = call - np.exp(-q * T) * S + K * np.exp(-r * T)
        if (option_type == "c"):
            return call
        elif (option_type == "p"):
            return put
        else:
            print("please enter the option type: (c/p)")
        pass

    def MCPricer(self, option_type='c', isAmerican=False):
        price_matrix = self.price_matrix
        n_steps = self.n_steps
        n_trials = self.n_trials
        strike = self.K
        risk_free_rate = self.r
        time_to_maturity = self.T
        dt = time_to_maturity / n_steps
        if (option_type == "c"):
            payoff_fun = lambda x: np.maximum(x-strike,0)
        elif (option_type == "p"):
            payoff_fun = lambda x: np.maximum(strike-x, 0)
        else:
            print("please enter the option type: (c/p)")
            return

        if (isAmerican == False):

            payoff = payoff_fun(price_matrix[:, n_steps])
            #         vk = payoff*df
            value_results = payoff * np.exp(-risk_free_rate * time_to_maturity)
            regular_mc_price = np.average(value_results)
            self.payoff = payoff
            self.mc_price = regular_mc_price
            self.value_results = value_results
        else:
            exercise_matrix = self.exercise_matrix
            t_exercise_array = dt * np.where(exercise_matrix == 1)[1]
            value_results = payoff_fun(price_matrix[np.where(exercise_matrix == 1)]) * np.exp(-risk_free_rate * t_exercise_array)
            regular_mc_price = np.average(value_results)
            self.mc_price = regular_mc_price
            self.value_results = value_results
        return regular_mc_price

    def BSDeltaHedgedPricer(self, option_type="c"):

        regular_mc_price = self.MCPricer(option_type=option_type)
        dt = self.T / self.n_steps
        df2 = np.exp(-(self.r - self.q) * dt)

        # Delta hedged cash flow
        def Delta_fun(x, tau, option_type):
            d1 = (np.log(x / self.K) + (self.r - self.q) * tau + self.sigma ** 2 * tau / 2) / (
                    self.sigma * np.sqrt(tau))
            if (option_type == 'c'):
                return sp.stats.norm.cdf(d1)
            elif (option_type == 'p'):
                return -sp.stats.norm.cdf(-d1)

        discounted_hedge_cash_flow = np.zeros(self.n_trials)
        for i in range(self.n_trials):
            Sk_array = self.price_matrix[i, :]
            bi_diag_matrix = np.diag([-1] * (self.n_steps), 0) + np.diag([df2] * (self.n_steps - 1), 1)
            # (Sk+1 exp(-r dt) - Sk) exp(-r*(tk-t0))
            discounted_stock_price_change = np.dot(bi_diag_matrix, Sk_array[:-1])
            discounted_stock_price_change[-1] += Sk_array[-1] * df2
            discounted_stock_price_change *= np.exp(-self.r * np.arange(self.n_steps) * dt)
            tau_array = dt * np.arange(self.n_steps, 0, -1)
            Delta_array = np.array([Delta_fun(Sk, tau, option_type) for Sk, tau in zip(Sk_array[:-1], tau_array)])
            discounted_hedge_cash_flow[i] = np.dot(Delta_array, discounted_stock_price_change)

        BSDeltaBased_mc_price = regular_mc_price - discounted_hedge_cash_flow.mean()
        #         print("The average discounted hedge cash flow: {}".format(discounted_hedge_cash_flow.mean()))

        value_results = self.payoff * np.exp(-self.r * self.T) - discounted_hedge_cash_flow
        #         print("Sanity check {} = {}".format(value_results.mean(),BSDeltaBased_mc_price))
        self.value_results = value_results

        return BSDeltaBased_mc_price

    def OHMCPricer(self, option_type='c', isAmerican=False, func_list=[lambda x: x ** 0, lambda x: x]):
        def _calculate_Q_matrix(S_k, S_kp1, df, df2, func_list):
            dS = df2 * S_kp1 - S_k
            A = np.array([func(S_k) for func in func_list]).T
            B = (np.array([func(S_k) for func in func_list]) * dS).T
            return np.concatenate((-A, B), axis=1)

        price_matrix = self.price_matrix
        # k = n_steps
        dt = self.T / self.n_steps
        df = np.exp(- self.r * dt)
        df2 = np.exp(-(self.r - self.q) * dt)
        n_basis = len(func_list)
        n_trials = self.n_trials
        n_steps = self.n_steps
        strike = self.K

        if (option_type == "c"):
            payoff_fun = lambda x: np.maximum(x-strike,0)
            # payoff = (price_matrix[:, n_steps] - strike)
        elif (option_type == "p"):
            payoff_fun = lambda x: np.maximum(strike-x,0)
            # payoff = (strike - price_matrix[:, n_steps])
        else:
            print("please enter the option type: (c/p)")
            return

        if isAmerican is True:
            holding_matrix = self.holding_matrix
        else:
            holding_matrix = np.ones(price_matrix.shape,dtype=bool)

        # At maturity
        holding_filter_k = holding_matrix[:, n_steps]
        payoff = matrix(payoff_fun(price_matrix[holding_filter_k,n_steps]))
        vk = payoff * df
        Sk = price_matrix[holding_filter_k,n_steps]
        #         print("regular MC price",regular_mc_price)

        # k = n_steps-1,...,1
        for k in range(n_steps - 1, 0, -1):

            holding_filter_kp1 = holding_filter_k
            holding_filter_k = holding_matrix[:, k]
            Skp1 = price_matrix[holding_filter_kp1, k+1]
            Sk = price_matrix[holding_filter_kp1, k]
            Qk = matrix(_calculate_Q_matrix(Sk, Skp1, df, df2, func_list))
            P = Qk.T * Qk
            q = Qk.T * vk
            A = matrix(np.ones(holding_filter_kp1.sum(), dtype=np.float64)).T * Qk
            b = - matrix(np.ones(holding_filter_kp1.sum(), dtype=np.float64)).T * vk
            # print(Sk)
            # print(Skp1)

            sol = solvers.coneqp(P=P, q=q, A=A, b=b)
            ak = sol["x"][:n_basis]
            bk = sol["x"][n_basis:]
            vk = matrix(np.array([func(price_matrix[holding_filter_k, k]) for func in func_list])).T * ak * df
            # break

        # k = 0
        v0 = vk
        holding_filter_1 = holding_filter_k
        holding_filter_0 = holding_matrix[:, 0]
        S0 = price_matrix[holding_filter_1, 0]
        S1 = price_matrix[holding_filter_1, 1]
        dS0 = df2 * S1 - S0
        Q0 = np.concatenate((-np.ones(holding_filter_1.sum())[:, np.newaxis], dS0[:, np.newaxis]), axis=1)
        Q0 = matrix(Q0)
        P = Q0.T * Q0
        q = Q0.T * v0
        A = matrix(np.ones(holding_filter_1.sum(), dtype=np.float64)).T * Q0
        b = - matrix(np.ones(holding_filter_1.sum(), dtype=np.float64)).T * v0
        C1 = matrix(ak).T * np.array([func(S1) for func in func_list]).T
        sol = solvers.coneqp(P=P, q=q, A=A, b=b)
        self.sol = sol
        residual_risk = (v0.T * v0 + 2 * sol["primal objective"]) / holding_filter_1.sum()
        self.residual_risk = residual_risk[0]  # the value of unit matrix

        return sol["x"][0]

    def standard_error(self):
        # can not apply to the OHMC since its result is not obtained by averaging
        # sample variance
        sample_var = np.var(self.value_results, ddof=1)
        std_estimate = np.sqrt(sample_var)
        standard_err = std_estimate / np.sqrt(self.n_trials)
        return standard_err

    def pricing(self, option_type='c', func_list=[lambda x: x ** 0, lambda x: x]):
        OHMC_price = self.OHMCPricer(option_type=option_type, func_list=func_list)
        regular_mc_price = self.MCPricer(option_type=option_type)
        black_sholes_price = self.BlackScholesPricer(option_type)
        return ({"OHMC": OHMC_price, "regular MC": regular_mc_price, "Black-Scholes": black_sholes_price})

    def hedging(self):
        S = self.S0
        K = self.K
        T = self.T
        r = self.r
        q = self.q
        sigma = self.sigma
        d1 = (np.log(S / K) + (r - q) * T + 0.5 * sigma ** 2 * T) / (sigma * np.sqrt(T))
        d2 = d1 - sigma * np.sqrt(T)
        N = lambda x: sp.stats.norm.cdf(x)
        return ({"OHMC optimal hedge": -self.sol["x"][1], "Black-Scholes delta hedge": N(d1),
                 "OHMC residual risk": self.residual_risk})
