Report – eWHORM simulations
================

## Transformation of the data

**Step 1:**

- Let’s consider $x_0, x_1$ the variable denoting microfilaria at
  baseline and at twelve months and $sd_0$ the common standard deviation
  for a specific treatment.

- $\mu_0 = mean(x_0)$, and $\mu_1 = (1 - r_0) \mu_0$

**Step 2:**

- we consider the $log(x_1)$ and $log(x_0)$ to calculate the difference:
  $\mu = log(x_1) - log(x_0)$ with

the following assumptions:

- $\tilde{x_i} = log(x_i)$ with the estimated mean $m_i$ and the
  variance $\sigma_{\tilde{x_i}}^2$

- $(\tilde{x_0}, \tilde{x_1})$ following bivariate normal distribution
  with mean

$m = (m_0, m_1)$ and variance-covariance matrix with elements
$\sigma_{\tilde{x_0}}^2, \sigma_{\tilde{x_1}}^2, and \quad cov(\tilde{x_1}, \tilde{x_0}) =cov(\tilde{x_0}, \tilde{x_1}) = \rho \sigma_{\tilde{x_0}} \sigma_{\tilde{x_1}}$

- From the known formula we obtained
  $m_i = log(\frac{\mu_i^2}{ \sqrt{\mu_i^2 + mean(sd_i)^2} })$,
  $\sigma_{\tilde{x_i}} = \sqrt{log(1 + \frac{sd_i^2}{\mu_i^2} )})$,
  $sd_0 = sd_i$, $\forall i$

- The difference of the log is following the normal distribution with
  mean = $m_1 - m_0$, var =
  $\sigma_{\tilde{x_0}}^2 + \sigma_{\tilde{x_1}}^2 - 2 \rho \sigma_{\tilde{x_1}} \sigma_{\tilde{x_0}} = \sigma_{\tilde{x_0}}^2$
  ($\rho = \frac{1}{2}, \sigma_{\tilde{x_0}}^2 = \sigma_{\tilde{x_i}}^2)$

## Values for the simulations

baseline means, sd, correlations, …

## Considered parameters for the design

Values

alpha1, alpha, power, sample sizes, …

- code examples

Let’s consider $mean(x_0) = 650, sd_0 = 575$
