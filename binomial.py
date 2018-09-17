import numpy as np
import math


def call(S, K, sigma, r, t, div=0, n=100, am=False):
    """
    Price a call option using the Binomial Options Pricing model
    S: initial spot price of stock
    K: strick price of option
    sigma: volatility
    r: risk-free interest rate
    t: time to maturity (in years)
    div: dividend yield (continuous compounding)
    n: binomial steps
    am: True for American option, False for European option
    %timeit results: 0 loops, best of 3: 49.6 ms per loop
    """
    return option(S, K, sigma, r, t, div, 1, n, am)


def put(S, K, sigma, r, t, div=0, n=100, am=False):
    """
    Price a put option using the Binomial Options Pricing model
    S: initial spot price of stock
    K: strick price of option
    sigma: volatility
    r: risk-free interest rate
    t: time to maturity (in years)
    div: dividend yield (continuous compounding)
    n: binomial steps
    am: True for American option, False for European option
    %timeit results: 10 loops, best of 3: 50.5 ms per loop
    """
    return option(S, K, sigma, r, t, div, -1, n, am)


def option(S, K, sigma, r, t, div=0, call=1, n=100, am=False):
    """
    Price an option using the Binomial Options Pricing model
    S: initial spot price of stock
    K: strick price of option
    sigma: volatility
    r: risk-free interest rate
    t: time to maturity (in years)
    div: dividend yield (continuous compounding)
    call: 1 if call option, -1 if put
    n: binomial steps
    am: True for American option, False for European option
    """
    delta = float(t) / n
    u = math.exp(sigma * math.sqrt(delta))
    d = float(1) / u
    q = float((math.exp((r - div) * delta) - d)) / (u - d)  # Prob. of up step
    stock_val = np.zeros((n + 1, n + 1))
    opt_val = np.zeros((n + 1, n + 1))

    # Calculate stock value at maturity
    stock_val[0, 0] = S
    for i in range(1, n + 1):
        stock_val[i, 0] = stock_val[i - 1, 0] * u
        for j in range(1, i + 1):
            stock_val[i, j] = stock_val[i - 1, j - 1] * d

    # Recursion for option price
    for j in range(n + 1):
        opt_val[n, j] = max(0, call*(stock_val[n, j] - K))
    for i in range(n - 1, -1, -1):
        for j in range(i + 1):
            opt_val[i, j] = \
                (q * opt_val[i + 1, j] + (1 - q) * opt_val[i + 1, j + 1]) \
                / math.exp(r * delta)
            if am:
                opt_val[i, j] = max(opt_val[i, j], call*(stock_val[i, j] - K))
    return opt_val[0, 0]


def call2(S, K, sigma, r, t, div=0, n=100, am=False):
    """
    Price a call option using the Binomial Options Pricing model
    S: initial spot price of stock
    K: strick price of option
    sigma: volatility
    r: risk-free interest rate
    t: time to maturity (in years)
    div: dividend yield (continuous compounding)
    n: binomial steps
    %timeit results: 10000 loops, best of 3: 139 us per loop
    """
    return option2(S, K, sigma, r, t, div, 1, n, am)


def put2(S, K, sigma, r, t, div=0, n=100, am=False):
    """
    Price a put option using the Binomial Options Pricing model
    S: initial spot price of stock
    K: strick price of option
    sigma: volatility
    r: risk-free interest rate
    t: time to maturity (in years)
    div: dividend yield (continuous compounding)
    n: binomial steps
    %timeit results: 10000 loops, best of 3: 136 us per loop
    """
    return option2(S, K, sigma, r, t, div, -1, n, am)


def option2(S, K, sigma, r, t, div=0, call=1, n=100, am=False):
    """
    Price an option using the Binomial Options Pricing model
    S: initial spot price of stock
    K: strick price of option
    sigma: volatility
    r: risk-free interest rate
    t: time to maturity (in years)
    div: dividend yield (continuous compounding)
    n: binomial steps
    """
    delta = float(t) / n
    u = math.exp(sigma * math.sqrt(delta))
    d = float(1) / u
    pu = float((math.exp((r - div) * delta) - d)) / (u - d)  # Prob. of up step
    pd = 1 - pu  # Prob. of down step
    u_squared = u * u
    S = S * pow(d, n)  # stock price at bottom node at last date
    prob = pow(pd, n)  # prob. of bottom node at last date
    opt_val = prob * max(0, call*(S - K))
    for i in range(1, n):
        S = S * u_squared
        prob = prob * (float(pu) / pd) * (n - i + 1) / i
        opt_val = opt_val + prob * max(call*(S - K), 0)
    return math.exp(-r * t) * opt_val