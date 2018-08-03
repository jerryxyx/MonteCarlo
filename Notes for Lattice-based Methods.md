# FE621 Notes

## Lattice-based Methods

Tree methods equations are derived by moment matching for mean and variance. Note that real world probabilities do not play any role in valuing an option in a tree method.

### Binomial Trees

* Additive tree ($log S_t$)

  * Trigeogis tree: $\Delta x_u = \Delta x_d$
    $$
    p_u \Delta x - p_d \Delta x = \mu \Delta t\\
    p_u \Delta x^2 + p_d \Delta x ^2 - (\mu \Delta t)^2 = \sigma^2 \Delta t\\
    p_u + p_d = 1\\
    \mu = r-\frac{\sigma^2}{2}
    $$

  * Jarrow Rudd tree: $p_u = p_d$

* Multiplicative tree ($S_t$)

  * Cox-Ross-Rubinstein (CRR) tree:  similar to Trigeogis tree, but different in the first moment  and involves approximation. $u = e^{\sigma \sqrt{\Delta t}}$
    $$
    p_u u + p_d d = e^{r \Delta t}\\
    p_u \Delta x^2 + p_d \Delta x ^2 - (\mu \Delta t)^2 = \sigma^2 \Delta t\\
    u = \frac{1}{d}=e^{\Delta x}\\
    p_u + p_d = 1\\
    \mu = r-\frac{\sigma^2}{2}
    $$










### Trinomial Trees

$$
p_u \Delta x - p_d \Delta x = \mu \Delta t\\
p_u \Delta x^2 + p_d \Delta x ^2 - (\mu \Delta t)^2 = \sigma^2 \Delta t\\
p_u +p_m+ p_d = 1\\
\mu = r-\frac{\sigma^2}{2}
$$

In order to make probabilities to be numbers between 0 and 1, we have a sufficient condition: 
$$
\Delta x > \sigma \sqrt{3 \Delta t}
$$
Any $\Delta x$ with this property produces a convergent tree.

### Finite Difference Method

More general to the addictive trinomial tree method.

### Convergence Comparison

| method             | rate of convergence                                          | convergence condition                    |
| :----------------- | ------------------------------------------------------------ | ---------------------------------------- |
| binomial tree      | $O((\Delta x)^2 + \Delta t)$                                 | NA                                       |
| trinomial tree     | $O((\Delta x)^2 + \Delta t)$                                 | $\Delta x \geq \sigma \sqrt{3 \Delta t}$ |
| explicit FDM       | $O((\Delta x)^2 + \Delta t)$                                 | $\Delta x \geq \sigma \sqrt{3 \Delta t}$ |
| implicit FDM       | $O((\Delta x)^2 + \Delta t)$                                 | stable                                   |
| Crank-Nicolson FDM | $O((\Delta x)^2 + (\frac{\Delta t}{2})^2)$                   | stable                                   |
| Monte Carlo        | $O\left(max\left(\Delta t, \frac{\sigma}{\sqrt{N_x}}\right)\right)$ | stable                                   |



## Variance Reduction

|      method        |              explanation                    |
| :-----------------:| :-----------------------------------------: |
| antithetic variates| the payoff of the antithetic pair $(X_1,X_2)$, $f(X_1), f(X_2)$ is negative correlated. A sufficient condition is to make payoff function monotone |
| delta-based control variates |                                                              |





## Risk-Neutral Measure

Risk-neutral measures make it easy to express the value of a derivative in a formula. 

$$H_0 = P(0,T) E_Q[H_T]$$

where the risk-neutral measure is denoted by Q. This can be re-stated in terms of the physical measure P as

$$H_0 = P(0,T)E_P[\frac{dQ}{dP} H_T]$$

Another name for the risk-neutral measure is the equivalent martingale measure. If there is just one unique risk-neutral measure in the market, then there is a unique arbitrage-free price for each asset in the market. This is the fundamental theorem of arbitrage-free pricing.If there are more such measures, then in an interval of prices no arbitrage is possible. If no equivalent martingale measure exists, arbitrage opportunities do. 

Suppose our economy consists of 2 assets, a stock and a risk-free bond, and that we use Black-Scholes model. In the model the evolution of the stock price can be described by Geometric Brownian Motion:

$$dS_t = \alpha S_t dt + \sigma S_t dW_t$$

where $W_t$ is a standard Brownian motion with respect to the physical measure. If we define

$$\tilde{W}_t = W_t + \frac{\alpha-r}{\sigma}t$$

Girsanov's theorem states that there exists a measure $Q$ under which ${\displaystyle {\tilde {W}}_{t}}$ is a Brownian motion. Put this back in the original equation:

$$dS_t = r S_t dt + \sigma S_t d\tilde{W}_t$$

Then, the discounted stock price $\tilde{S}_t$is a $Q$-martingale.

$$d \tilde{S}_t = \sigma \tilde{S}_t d\tilde{W}_t$$

Note that risk neutral measure is powerful because you don't need to replicate a portfolio in order to be arbitrage-free compared to risk-free bond. Under such measure, expected value(first moment) of securities are equal to rolling it into a deposit account (inflate it to the maturity date at the riskless rate of interest). 