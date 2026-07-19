# Volatility Surface Exploration

Working through stochastic volatility modeling from the underlying stochastic calculus
to a calibrated implied volatility surface, following Gatheral's *The Volatility
Surface*. Each model is derived in full and then implemented, so the code and the math
can be checked against each other.

## Heston (1993)

### Derivation

`Heston/heston_derivation.md` is a complete derivation, from the model definition to a
Fourier pricing formula:

1. Motivation — why constant volatility and local volatility both fail
2. The model under the risk-neutral measure
3. The CIR variance process and its stationary behavior
4. The Feller condition
5. Log-price dynamics
6. The Kolmogorov backward equation
7. Affine structure and the exponential-affine ansatz
8. Deriving the Riccati ODEs
9. Solving the Riccati equation
10. The characteristic function
11. Fourier pricing
12. Parameter intuition on the implied volatility surface
13. Connections and limitations

The through-line is why Heston is tractable at all: log-price and variance form an
affine system, so the characteristic function has an exponential-affine form whose
coefficients solve Riccati ODEs in closed form. That gives semi-analytic option prices
by Fourier inversion instead of Monte Carlo.

The Feller condition, `2κθ > ξ²`, gets its own section because it decides whether
variance can reach zero — the boundary behavior that determines whether a discretization
scheme is even valid.

### Implementation

`Heston/heston_implemented.ipynb` implements the characteristic function, prices options
by numerical Fourier inversion (`scipy.integrate.quad`), inverts prices back to Black-
Scholes implied volatility, and calibrates `(κ, θ, ξ, ρ, v₀)` by numerical optimization.

Output figures:

| Figure | Shows |
| --- | --- |
| `heston_smiles.png` | Implied volatility smiles by maturity |
| `heston_surface.png` | The fitted implied volatility surface |
| `heston_params.png` | Effect of each parameter on surface shape |

`ρ` controls the skew (negative for equities, where downside protection is bid up), `ξ`
controls the curvature of the smile, and `κ` controls how quickly the smile flattens with
maturity.

## Reading the derivation

The math renders as LaTeX. GitHub's Markdown viewer handles most of it; for the full
rendering, open the file in a Markdown editor with math support.

## Requirements

```bash
pip install numpy scipy matplotlib
```
