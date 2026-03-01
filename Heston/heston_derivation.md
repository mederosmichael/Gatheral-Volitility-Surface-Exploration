# The Heston Stochastic Volatility Model: Complete Derivation

## 1. Motivation

Black-Scholes assumes constant volatility σ. This is empirically wrong — implied volatility varies by strike and maturity (the volatility surface). Local volatility (Dupire) fixes the marginal distributions but makes volatility deterministic given the spot path, producing unrealistic dynamics. Heston (1993) introduces genuine randomness in variance itself.

---

## 2. The Model

Under the risk-neutral measure:

$$dS_t = rS_t\,dt + \sqrt{v_t}\,S_t\,dW_t^{(1)}$$

$$dv_t = \kappa(\theta - v_t)\,dt + \xi\sqrt{v_t}\,dW_t^{(2)}$$

$$dW^{(1)}_t\,dW^{(2)}_t = \rho\,dt$$

**Parameters:**
- θ: long-run variance (where v_t mean-reverts to)
- κ: speed of mean reversion
- ξ: vol of vol (volatility of the variance process)
- ρ: correlation between spot and variance (typically negative for equities)
- v₀: initial variance

---

## 3. The CIR Variance Process

The variance follows a Cox-Ingersoll-Ross process:

$$dv_t = \kappa(\theta - v_t)\,dt + \xi\sqrt{v_t}\,dW_t$$

### 3.1 Mean Reversion

Without noise, the ODE dv/dt = κ(θ − v) has the solution:

$$v_t = \theta + (v_0 - \theta)e^{-\kappa t}$$

The gap to the long-run level decays exponentially with half-life ln(2)/κ. Adding noise creates a tug-of-war: noise knocks v_t away from θ, the drift pulls it back. The stationary variance is σ²/(2κ).

### 3.2 The √v Diffusion Scaling

Why √v_t and not a constant or v_t?

- **Constant diffusion (σ):** Noise doesn't vanish near zero → process goes negative. Unacceptable for variance.
- **Linear diffusion (σv):** Noise vanishes too fast near zero → process can be absorbed at zero. Also too aggressive away from zero.
- **Square root (σ√v):** Near zero, √v vanishes slower than v (since √v ≫ v for small v). The mean-reverting drift κθ can win against the shrinking noise, pushing the process back. But the noise does shrink, preventing the process from being blasted through zero.

The √v scaling is the Goldilocks choice: vanishes fast enough to prevent negativity, slow enough that the process doesn't get stuck at zero.

---

## 4. The Feller Condition

**Statement:** If 2κθ ≥ ξ², the variance process v_t never reaches zero (zero is unattainable). If 2κθ < ξ², v_t can touch zero but immediately reflects back.

### 4.1 The Scale Function

**Goal:** Transform v_t into a process with no drift (a local martingale) to analyze boundary behavior.

We seek s(v) such that s(v_t) is a local martingale. Apply Itô's formula to s(v_t):

$$ds(v_t) = \frac{ds}{dv}\,dv_t + \frac{1}{2}\frac{d^2s}{dv^2}(dv_t)^2$$

Substitute dv_t = κ(θ − v)dt + ξ√v dW:

**First term:** (ds/dv)[κ(θ − v)dt + ξ√v dW]

**Second term:** (1/2)(d²s/dv²) · ξ²v dt   [since (dv)² = ξ²v dt by Itô rules]

Collecting dt and dW terms:

$$ds(v_t) = \left[\kappa(\theta - v)\frac{ds}{dv} + \frac{1}{2}\xi^2 v\frac{d^2s}{dv^2}\right]dt + \xi\sqrt{v}\frac{ds}{dv}\,dW$$

For s(v_t) to be a local martingale, the dt coefficient must vanish:

$$\kappa(\theta - v)\frac{ds}{dv} + \frac{1}{2}\xi^2 v\frac{d^2s}{dv^2} = 0$$

### 4.2 Solving the Scale Function ODE

Substitute p(v) = ds/dv to get a first-order ODE:

$$\kappa(\theta - v)\,p + \frac{1}{2}\xi^2 v\,\frac{dp}{dv} = 0$$

Rearrange:

$$\frac{dp}{dv} = -\frac{2\kappa(\theta - v)}{\xi^2 v}\,p$$

This is separable. Divide both sides by p:

$$\frac{1}{p}\frac{dp}{dv} = -\frac{2\kappa(\theta - v)}{\xi^2 v}$$

The left side is (d/dv)(ln p). The right side splits:

$$= -\frac{2\kappa\theta}{\xi^2}\cdot\frac{1}{v} + \frac{2\kappa}{\xi^2}$$

Define α = 2κθ/ξ² and β = 2κ/ξ². Integrate:

$$\ln p(v) = -\alpha\ln v + \beta v + \text{const}$$

Exponentiate:

$$p(v) = \frac{ds}{dv} = v^{-\alpha}\,e^{\beta v}$$

### 4.3 Boundary Classification via Hitting Probabilities

**Setup:** s(v_t) is a local martingale. For an interval (a, b), let τ = min(τ_a, τ_b). The stopped process s(v_{t∧τ}) is bounded on [a, b], hence a true martingale. Apply the optional stopping theorem:

$$\mathbb{E}[s(v_\tau)] = s(v_0)$$

Since v_τ equals either a or b:

$$s(a)(1 - p) + s(b)\cdot p = s(v_0)$$

where p = P(hit b before a). Solving:

$$P(\text{hit } a \text{ before } b) = \frac{s(b) - s(v_0)}{s(b) - s(a)}$$

**The Optional Stopping Theorem** says E[M_τ] = M₀ for a martingale M stopped at a bounded stopping time τ. The proof decomposes M_τ − M₀ into a telescoping sum of martingale increments, each with zero conditional expectation. The key fact is that the indicator 1_{τ≥k} is F_{k−1}-measurable (since τ is a stopping time, the decision to keep playing at step k depends only on information up to time k−1), so it factors out of conditional expectations.

### 4.4 Testing Whether Zero is Reachable

Let a → 0⁺ in the hitting probability formula. The numerator s(b) − s(v₀) stays finite. The denominator s(b) − s(a) depends on s(0⁺).

Since s(a) = −∫_a^c (ds/dv) dv → −∫_0^c v^{−α}e^{βv} dv as a → 0:

**If ∫_0^c v^{−α}e^{βv} dv = ∞:** Then s(0⁺) = −∞, the denominator → ∞, and P(hit 0) → 0. Zero is unreachable — it is infinitely far away in the intrinsic geometry of the process.

**If ∫_0^c v^{−α}e^{βv} dv < ∞:** Then s(0⁺) is finite, and P(hit 0) > 0. Zero is reachable.

**The integral's meaning:** The scale function s warps the coordinate system to reflect how the process actually experiences distance. If s(0⁺) = −∞, zero looks close in v-space but is infinitely far in the process's natural geometry, accounting for the drift pushing away and the noise shrinking.

### 4.5 The Convergence Test

Near v = 0, e^{βv} → 1, so the integrand behaves like v^{−α}. The integral ∫_0^c v^{−α} dv has antiderivative v^{1−α}/(1−α):

- **α < 1:** v^{1−α} → 0 as v → 0. Integral converges. Zero is reachable.
- **α = 1:** ∫v^{−1} dv = ln v → −∞. Diverges. Zero is unreachable.
- **α > 1:** v^{1−α} → ∞ as v → 0. Diverges. Zero is unreachable.

Converges iff α < 1, i.e., 2κθ/ξ² < 1, i.e., 2κθ < ξ².

**Feller condition:** Zero is unreachable iff 2κθ ≥ ξ².

### 4.6 Boundary Type When Feller is Violated

The speed measure density m(v) = 1/(σ²(v) · s'(v)) = v^{α−1}/(ξ²e^{βv}). Since ∫_0^c v^{α−1} dv converges for α > 0 (always true), zero is a **regular boundary** when attainable — the process reaches it and we can impose reflection. In Heston, variance touches zero and immediately bounces back.

When Feller holds, zero is an **entrance boundary** — the process can start there but can never reach it from the interior.

---

## 5. Log-Price Dynamics

Apply Itô's formula to ln(S_t) with dS = rS dt + √v S dW⁽¹⁾:

$$d\ln S = \frac{1}{S}dS + \frac{1}{2}\left(-\frac{1}{S^2}\right)(dS)^2$$

$$= \frac{1}{S}(rS\,dt + \sqrt{v}S\,dW^{(1)}) + \frac{1}{2}\left(-\frac{1}{S^2}\right)v S^2\,dt$$

$$= \left(r - \frac{v}{2}\right)dt + \sqrt{v}\,dW^{(1)}$$

The −v/2 term is the Itô correction: the same reason geometric Brownian motion's log has drift r − σ²/2, except σ² is replaced by the stochastic v_t.

Define x_t = ln S_t. The system becomes:

$$dx_t = \left(r - \frac{v_t}{2}\right)dt + \sqrt{v_t}\,dW_t^{(1)}$$

$$dv_t = \kappa(\theta - v_t)\,dt + \xi\sqrt{v_t}\,dW_t^{(2)}$$

---

## 6. The Kolmogorov Backward Equation

### 6.1 Why Conditional Expectations Satisfy PDEs

The characteristic function φ(ω, x, v, t) = E[e^{iωx_T} | x_t = x, v_t = v] is a conditional expectation. By the tower property, conditional expectations are martingales:

$$\mathbb{E}[\mathbb{E}[g(X_T)|\mathcal{F}_t] | \mathcal{F}_s] = \mathbb{E}[g(X_T)|\mathcal{F}_s] \quad \text{for } s < t$$

This is because the tower property says conditioning on coarse information (F_s) after conditioning on fine information (F_t) gives the coarser conditioning. Hence E[M_t | F_s] = M_s — the martingale property.

Since φ(x_t, v_t, t) is a martingale, apply Itô and set the drift to zero — the same principle used for the scale function and Black-Scholes.

### 6.2 The Two-Dimensional Backward PDE

Apply two-dimensional Itô to φ(x, v, t). The second-order Taylor expansion in stochastic calculus keeps all terms up to order dt:

- (dx)² = v dt [survives]
- (dv)² = ξ²v dt [survives]  
- (dx)(dv) = √v · ξ√v · ρ dt = ρξv dt [survives because dW⁽¹⁾dW⁽²⁾ = ρ dt]

Setting the drift to zero gives:

$$\frac{\partial\phi}{\partial t} + \left(r - \frac{v}{2}\right)\frac{\partial\phi}{\partial x} + \kappa(\theta - v)\frac{\partial\phi}{\partial v} + \frac{1}{2}v\frac{\partial^2\phi}{\partial x^2} + \frac{1}{2}\xi^2 v\frac{\partial^2\phi}{\partial v^2} + \rho\xi v\frac{\partial^2\phi}{\partial x\,\partial v} = 0$$

Terminal condition: φ(ω, x, v, T) = e^{iωx}.

This is the Black-Scholes PDE generalized to two dimensions. Each coefficient comes directly from the SDE: drifts give first-order terms, squared diffusions give second-order terms, correlation gives the cross term.

---

## 7. The Affine Structure and Ansatz

### 7.1 What "Affine" Means

The model is affine because all drift and squared-diffusion coefficients are linear (a + bv) in the state variable v:

- Drift of x: r − v/2 (affine in v)
- Drift of v: κθ − κv (affine in v)
- (Diffusion of x)²: v (affine in v)
- (Diffusion of v)²: ξ²v (affine in v)
- Cross: ρξv (affine in v)

No v², no √v in drift, nothing nonlinear. This structural property guarantees that the characteristic function is exponential-affine in the state variables.

### 7.2 The Ansatz

We guess:

$$\phi(\omega, x, v, t) = \exp\left(C(\omega, \tau) + D(\omega, \tau)\,v + i\omega\,x\right)$$

where τ = T − t, with C(0) = 0 and D(0) = 0.

This is not a blind guess — for affine models, the Duffie-Pan-Singleton (2000) theorem guarantees this form works. The exponent is affine in the state variables (x, v), and C, D satisfy ODEs (not PDEs) because the affine structure allows the state variables to separate.

---

## 8. Deriving the Riccati ODEs

### 8.1 Partial Derivatives of the Ansatz

Let f = C(τ) + D(τ)v + iωx, so φ = e^f. Since τ = T − t, ∂τ/∂t = −1.

| Derivative | Result |
|-----------|--------|
| ∂φ/∂t | (−dC/dτ − (dD/dτ)v) φ |
| ∂φ/∂x | iω φ |
| ∂φ/∂v | D φ |
| ∂²φ/∂x² | −ω² φ |
| ∂²φ/∂v² | D² φ |
| ∂²φ/∂x∂v | iωD φ |

### 8.2 Substitution into the Backward PDE

Substitute, then divide by φ (which appears in every term):

$$-\frac{dC}{d\tau} - \frac{dD}{d\tau}v + \left(r - \frac{v}{2}\right)(i\omega) + \kappa(\theta - v)D + \frac{1}{2}v(-\omega^2) + \frac{1}{2}\xi^2 v D^2 + \rho\xi v(i\omega D) = 0$$

Expand:

$$-\frac{dC}{d\tau} - \frac{dD}{d\tau}v + i\omega r - \frac{i\omega v}{2} + \kappa\theta D - \kappa v D - \frac{\omega^2 v}{2} + \frac{\xi^2 D^2 v}{2} + \rho\xi i\omega D v = 0$$

### 8.3 Separation

This must hold for all v. Since the equation is A + Bv = 0 for all v, both A = 0 and B = 0:

**Constant terms (ODE for C):**

$$\frac{dC}{d\tau} = i\omega r + \kappa\theta D$$

**Coefficient of v (ODE for D):**

$$\frac{dD}{d\tau} = \frac{\xi^2}{2}D^2 + (\rho\xi i\omega - \kappa)D - \frac{1}{2}(i\omega + \omega^2)$$

The separation works precisely because every coefficient in the PDE was constant or linear in v (the affine property). A non-affine model would produce v² terms that can't be separated.

---

## 9. Solving the Riccati Equation

### 9.1 The Riccati → Linear ODE Reduction

The D equation has the form dD/dτ = aD² + bD + c where:
- a = ξ²/2
- b = ρξiω − κ
- c = −(iω + ω²)/2

Substitute D = −(1/a)(1/u)(du/dτ). This logarithmic derivative substitution linearizes the quadratic ODE.

**Computing dD/dτ:** By the quotient rule on u'/u:

$$\frac{dD}{d\tau} = -\frac{1}{a}\left(\frac{u''}{u} - \frac{(u')^2}{u^2}\right)$$

**Substituting into aD² + bD + c:** After cancellation of the (u')²/u² terms from both sides:

$$-\frac{1}{a}\frac{u''}{u} = -\frac{b}{a}\frac{u'}{u} + c$$

Multiply by −au:

$$u'' - bu' + acu = 0$$

This is a constant-coefficient linear second-order ODE.

### 9.2 Solving the Linear ODE

Characteristic equation: λ² − bλ + ac = 0

$$\lambda = \frac{b \pm d}{2}$$

where d = √(b² − 4ac). Define λ₁ = (b + d)/2 and λ₂ = (b − d)/2.

General solution: u(τ) = Ae^{λ₁τ} + Be^{λ₂τ}

### 9.3 Applying D(0) = 0

Since D = −(1/a)(u'/u), the condition D(0) = 0 requires u'(0) = 0.

u'(0) = Aλ₁ + Bλ₂ = 0, giving A = −(λ₂/λ₁)B.

Setting B = λ₁:

u(τ) = −λ₂e^{λ₁τ} + λ₁e^{λ₂τ}

### 9.4 The Solution for D

After computing D = −(1/a)(u'/u) and simplifying with g = (b − d)/(b + d):

$$D(\tau) = \frac{b - d}{\xi^2}\cdot\frac{1 - e^{-d\tau}}{1 - g\,e^{-d\tau}}$$

where:

$$d = \sqrt{(\kappa - \rho\xi i\omega)^2 + \xi^2(i\omega + \omega^2)}$$

### 9.5 The Solution for C

Integrate dC/dτ = iωr + κθD:

$$C(\tau) = i\omega r\tau + \frac{\kappa\theta}{\xi^2}\left[(b - d)\tau - 2\ln\left(\frac{1 - ge^{-d\tau}}{1 - g}\right)\right]$$

---

## 10. The Characteristic Function

Combining everything:

$$\phi(\omega) = \exp\left(C(\omega, \tau) + D(\omega, \tau)\,v_0 + i\omega\,x_0\right)$$

with C and D given above. This is a closed-form expression for E[e^{iω ln S_T}], encoding the entire distribution of the log-price under the Heston model.

---

## 11. Fourier Pricing

### 11.1 Decomposing the Call Price

$$C_{\text{call}} = e^{-rT}\mathbb{E}[\max(S_T - K, 0)] = e^{-rT}\mathbb{E}[S_T\mathbf{1}_{\{S_T > K\}}] - Ke^{-rT}\mathbb{P}(S_T > K)$$

The second term needs P(ln S_T > ln K) = P₂.

The first term uses a change of measure: define a tilted probability P̃ where P̃(A) = E[S_T 1_A]/E[S_T]. Under this measure, the S_T weighting becomes a plain probability:

$$e^{-rT}\mathbb{E}[S_T\mathbf{1}_{\{S_T > K\}}] = S_0\,\widetilde{\mathbb{P}}(S_T > K) = S_0 P_1$$

### 11.2 The Gil-Pelaez Inversion Formula

Both P₁ and P₂ are tail probabilities recoverable from characteristic functions:

$$P_j = \frac{1}{2} + \frac{1}{\pi}\int_0^\infty \text{Re}\left[\frac{e^{-i\omega\ln K}\,\phi_j(\omega)}{i\omega}\right]d\omega$$

**Derivation sketch:**
1. Write P(X > k) = ∫_k^∞ f(x) dx
2. Replace f with Fourier inversion: f(x) = (1/2π) ∫ e^{−iωx} φ(ω) dω
3. Swap integration order (Fubini)
4. Evaluate the inner integral ∫_k^∞ e^{−iωx} dx = e^{−iωk}/(iω) (using a damping factor for convergence)
5. Handle the pole at ω = 0 (contributes the 1/2)
6. Use conjugate symmetry φ(−ω) = φ(ω)* to fold the integral from (−∞,∞) to (0,∞)

Here φ₂ = φ (the original characteristic function) and φ₁ = φ(ω − i)/(e^{rT}φ(−i)) (the characteristic function under the tilted measure).

### 11.3 The Final Formula

$$C_{\text{call}} = S_0 P_1 - Ke^{-rT}P_2$$

This has the same structure as Black-Scholes: S₀N(d₁) − Ke^{−rT}N(d₂). In BS, the probabilities are normal CDFs (because log-price is Gaussian). In Heston, the distribution is non-Gaussian, so the CDFs are replaced by Fourier integrals.

Numerically, each integral converges fast (φ decays at high frequencies) and can be evaluated with ~50-100 quadrature points. This makes Heston calibration fast.

---

## 12. Parameter Intuition on the Implied Volatility Surface

### ρ (Correlation) — Controls the Skew

ρ < 0: spot drops → variance increases → large downside moves amplified → OTM puts expensive → downside skew. More negative ρ = steeper skew. At ρ = 0, the smile is symmetric.

**ρ answers: "When the stock moves, does vol move with it or against it?"**

### ξ (Vol of Vol) — Controls the Curvature

High ξ → variance fluctuates widely → fat tails in both directions → both OTM puts and calls expensive → convex smile. At ξ = 0, variance is deterministic and the smile is flat (Black-Scholes).

**ξ answers: "How uncertain is future volatility?"**

### κ (Mean Reversion Speed) — Controls the Term Structure

High κ → vol shocks die quickly → long-dated options insensitive to current v₀ → smile flattens with maturity faster. ATM vol converges to √θ as T → ∞.

**κ answers: "How persistent are vol shocks?"**

### θ (Long-Run Variance) — Controls the Long-End Level

θ sets the asymptotic ATM implied vol for long-dated options. Short-dated ATM vol is driven by v₀, long-dated by θ.

**θ answers: "What is the market's long-run expectation of variance?"**

### v₀ (Initial Variance) — Controls the Short-End Level

v₀ > θ: ATM term structure slopes downward (vol is elevated, will revert). v₀ < θ: slopes upward. v₀ = θ: roughly flat.

---

## 13. Connections and Limitations

### Black-Scholes as a Special Case

Set ξ = 0 (no vol randomness). Then v_t = θ + (v₀ − θ)e^{−κt} is deterministic. In the limit κ → ∞, v_t = θ = σ² (constant), and Heston reduces to Black-Scholes with flat implied vol.

### Why Heston Beats Local Vol

Local volatility (Dupire) produces the correct marginal distributions but wrong dynamics: the forward smile flattens unrealistically, and vol is deterministic given the spot path. Heston's stochastic vol generates genuine randomness in vol, producing more realistic smile dynamics and better hedging.

### What Heston Can't Do

Heston can't independently control short-dated and long-dated smiles (same five parameters govern the entire surface). Short-dated equity smiles have more convexity than Heston produces. This motivates:
- **Bates model:** Heston + jumps in spot
- **Double Heston:** Two variance factors
- **Rough volatility:** Fractional Brownian motion driving variance
- **SVI parameterization (Gatheral):** More flexible fitting of individual maturity slices

---

## Appendix: Key Results Used

### Itô's Formula (1D)
For f(X_t) where dX = μ dt + σ dW:
$$df = \frac{df}{dx}dX + \frac{1}{2}\frac{d^2f}{dx^2}(dX)^2$$
The second-order term survives because (dW)² = dt (Brownian motion is rough enough that the quadratic variation is non-zero).

### Itô's Formula (2D)
For u(X_t, Y_t):
$$du = \frac{\partial u}{\partial x}dX + \frac{\partial u}{\partial y}dY + \frac{1}{2}\frac{\partial^2 u}{\partial x^2}(dX)^2 + \frac{1}{2}\frac{\partial^2 u}{\partial y^2}(dY)^2 + \frac{\partial^2 u}{\partial x\partial y}(dX)(dY)$$
The cross term (dX)(dY) survives when the driving Brownian motions are correlated.

### Itô Multiplication Rules
- dt · dt = 0
- dt · dW = 0
- dW · dW = dt
- dW⁽¹⁾ · dW⁽²⁾ = ρ dt

### Optional Stopping Theorem
If M_t is a martingale and τ is a bounded stopping time, then E[M_τ] = E[M₀]. A clever stopping strategy cannot beat a fair game.

### Kolmogorov Backward Equation
u(x, t) = E[g(X_T) | X_t = x] satisfies:
$$\frac{\partial u}{\partial t} + \mu(x)\frac{\partial u}{\partial x} + \frac{1}{2}\sigma^2(x)\frac{\partial^2 u}{\partial x^2} = 0$$
Derived by noting u is a martingale (tower property), applying Itô, and setting the drift to zero. The Black-Scholes PDE is a special case.

### Fourier Transform
$$\hat{f}(\omega) = \int_{-\infty}^{\infty}e^{i\omega x}f(x)\,dx \quad \iff \quad f(x) = \frac{1}{2\pi}\int_{-\infty}^{\infty}e^{-i\omega x}\hat{f}(\omega)\,d\omega$$
The characteristic function φ(ω) = E[e^{iωX}] is the Fourier transform of the density. It always exists (|e^{iωX}| = 1) and uniquely determines the distribution.
