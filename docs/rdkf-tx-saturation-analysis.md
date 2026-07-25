# RDKF vs DKF on 3D LFM data — TX saturation and RMSE deficit

Status: **diagnosed, not fixed (2026-07-25).** No repo code was modified. This
document records why RDKF transmits at 100% and still loses to DKF, what the
tuning sweeps do and do not tell us, and the checklist of next steps.

Update (2026-07-25, later): **checklist item 1 is resolved — see §6.** The
compounding deflation is faithful to the paper, so that remedy is off the table;
§6.1 replaces it with an exact closed-form feasibility bound that explains the
saturation and lets items 2–3 be computed instead of swept.

Runs analysed:

- `results/estimateSST3dLfm_T1000_N100s100_b0.05_blfm0.05_runs200_20260721-044004.mat`
  (sha `b540599`, 200 runs, T=1000, b = b_lfm = 0.05)
- `results/tuneRDKFbLfm_sst3d_T1000_N100s100_b0.05_blfm0.05_runs50_20260720-195223.mat`
  (sha `7fba732`, b sweep on LFM data, trigger held at α=10, β=0.2, δ=0.5)
- `results/tuneRDKF_sst3d_T1000_N100s100_b0.05_runs50_20260720-200410.mat`
  (sha `7fba732`, β×δ grid on **nominal** data, α=10 and b=0.05 fixed)
- `results/tuneDKF_sst3d_T1000_N100s100_b0.05_runs50_20260720-194926.mat` (reference)

All RMSE figures are the steady-state statistic (`utils/ssRmseStats.m`: mean over
the second half of the horizon, then mean across runs). TX figures are averaged
over the same window.

---

## 1. The symptom

| filter | ssRMSE | ±std | TX (ss) |
|---|---|---|---|
| CKF (centralized nominal) | 0.620 | 0.060 | — |
| **CRKF (centralized robust)** | **0.564** | 0.041 | — |
| **DKF** | **1.009** | 0.085 | **0.469** |
| **RDKF** | **1.123** | 0.104 | **1.0000** |
| RDKFLOC | 1.121 | 0.104 | 1.0000 |
| SRDKF-Open | 1.123 | 0.104 | 0.9989 |
| SRDKF-Closed | 1.379 | 0.119 | 0.345 |
| SDKF-Closed | 1.555 | 0.206 | 0.348 |
| SDKF-Open | 1.372 | 0.175 | 0.9989 |
| SRDKFLOC-Closed | 1.376 | 0.119 | 0.345 |
| SRDKFLOC-Open | 1.121 | 0.104 | 0.9989 |
| DSEACP | 4.498 | 0.225 | — |

RDKF's TX is not merely *near* 1 — it is **exactly** 1.0 at every step of every
one of the 200 runs (`frac(txRate == 1) = 1.0000` over the steady-state window).
After t=1 it never skips a single broadcast. So it pays 2.1× DKF's bandwidth and
still loses 11% on RMSE.

Two corollaries:

- RDKFLOC and SRDKF-Open are numerically indistinguishable from RDKF. With TX≈1
  all three degenerate to the same all-transmit robust filter, so those rows carry
  no independent information in this run.
- **Robustness works centrally**: CRKF beats CKF by 9% on the same LFM data. The
  defect is specific to the *distributed* robust construction, not to the LFM
  setup or the KL radius.

---

## 2. Root cause of the saturation

The RDKF trigger (`estimators/RDKF.m:139-160`) skips a broadcast iff

    eNorm = e' Ω_i e <= α        AND        Ω_i/(1+β) ⪯ Ψ̄_i ⪯ (1+δ) Ω_i

Instrumenting all three conditions per node (α=10, β=0.2, δ=0.5, steady-state
window):

| b | TX | P(eNorm≤α) | P(lower ⪯ ok) | P(upper ⪯ ok) | median λmin(Ω⁻¹Ψ̄) | θ/λmin(Ω_pred) |
|---|---|---|---|---|---|---|
| 0 | 0.505 | 0.9998 | **0.4949** | 1.0000 | 0.830 | 0 |
| 0.003 | 0.520 | 0.9996 | **0.4799** | 1.0000 | 0.797 | 0.094 |
| 0.01 | **1.0** | 0.9999 | **0.0000** | 1.0000 | 0.792 | 0.168 |
| 0.05 | **1.0** | 0.9993 | **0.0000** | 1.0000 | 0.634 | 0.337 |

### 2.1 α and δ are dead knobs

- The eNorm condition passes ~100% of the time. Median eNorm ≈ 0.4 against α=10 —
  α is ~25× too loose to ever bind.
- The upper Loewner bound passes 100% of the time: λmax ≈ 1.00 against
  1+δ = 1.5.

This is why the `tuneRDKF` grid reports *identical* RMSE (0.66067) for
δ ∈ {0.1, 0.5, 1, 2.5} at β=0.1 — δ only enters via the fusion fallback
(`RDKF.m:179-180`), which never executes when TX=1.

### 2.2 The lower Loewner bound is the only binding condition, and θI breaks it

Ψ̄ is the deflated pair, Ψ = Ω_pred − θI. At b=0.05, θ is 34% of λmin(Ω_pred),
which drags λmin(Ω⁻¹Ψ̄) down to 0.634 against a required 1/(1+β) = 0.833.

The distribution is razor-tight — p5 = 0.627, p95 ≈ 0.64. This is a
**deterministic bias, not variance**: every node fails, every step. Hence the
cliff between b=0.003 and b=0.01 in the sweep, where TX jumps 0.542 → 1.000 with
nothing in between.

The β required to admit the observed λmin at b=0.05 is **≈ 0.6**, i.e. 3× the
shipped 0.2.

### 2.3 The deflation compounds during silence — a hard TX floor

`updateGlobalPriors` (`RDKF.m:211-233`) re-deflates the *already-deflated* Ψ̄ on
every step a node stays silent: line 221 sets `Omega_i_check = Psi_bar(:,:,i)`,
and θ̄ is subtracted again at line 230. Correlating λmin against consecutive
silence age (β=1, δ=0.1):

| silence age | b=0.05: median λmin | frac failing lower | b=0: median λmin | frac failing |
|---|---|---|---|---|
| 0 | 0.712 | 0.0005 | 0.931 | 0 |
| 1 | **0.396** | **0.9993** | 0.889 | 0 |
| 2 | 0.256 | 1.0 | 0.759 | 0 |
| 3 | — | — | 0.747 | 0 |

At b=0.05 each silent step costs 44% of λmin, so **a node can essentially never be
silent twice in a row** — β shifts where the cliff sits but cannot remove it. At
b=0 the decay is gentle and asymptotes (0.93 → 0.89 → 0.76 → 0.75), never failing
the bound; silence is then limited by the other conditions instead.

This is a structural TX floor of ~50% for RDKF at b=0.05, independent of tuning.
**Whether the compounding is faithful to Ghion & Zorzi or an implementation
artifact is the single highest-leverage open question** (see checklist item 1).

---

## 3. What the tuning sweeps say

### 3.1 `tuneRDKFbLfm_sst3d` — b sweep on LFM data (β=0.2, δ=0.5)

| b | 0 | 1e-4 | 3e-4 | 1e-3 | 3e-3 | 1e-2 | 3e-2 | 5e-2 | 1e-1 |
|---|---|---|---|---|---|---|---|---|---|
| SsRMSE | 0.994 | **0.991** | 0.993 | 0.996 | 1.013 | 1.089 | 1.111 | 1.137 | 1.188 |
| TX | 0.502 | 0.511 | 0.516 | 0.530 | 0.542 | 1.0 | 1.0 | 1.0 | 1.0 |

- The reported "best" b=1e-4 beats the b=0 anchor by 0.2% against a std of 0.064.
  That is noise, not a signal. **The robust layer buys nothing at any b tested.**
- The rows b ∈ {0.01 … 0.1} are *all* at TX=1.0, so they form an accidental
  matched-bandwidth comparison: with the trigger fully neutralized, RMSE still
  rises monotonically 1.089 → 1.188 as b grows. **More robustness is strictly
  worse even when given every measurement.** This is the cleanest damning number
  in the set.
- The sweep is confounded by design: b moves robustness *and* transmission rate
  at once (bigger b → smaller λmin → more TX → better RMSE), so the low-b rows
  conflate the two effects. It never tests "meaningful b with a healthy trigger".

### 3.2 `tuneRDKF_sst3d` — β×δ grid, nominal data, b=0.05

Pareto set: β=0.1/any δ (TX 1.0, 0.661) · β=0.5/δ=0.1 (TX 0.710, 0.729) ·
β=1/δ=0.1 (TX 0.506, 0.788) · β=2.5/δ=0.1 (TX 0.339, 0.890).

- **δ=0.1 dominates at every β.** The fusion fallback discount 1/(1+δ) is
  expensive: at β=1, δ=0.1 → 0.788 vs δ=2.5 → 1.662. The shipped δ=0.5 is a
  straight loss.
- β must be ≥0.5 for the trigger to fire at all at b=0.05 — independently
  confirming §2.2 from a completely different measurement.

### 3.3 `tuneDKF_sst3d` reference

DKF reaches 0.531 @ TX 0.353 (β=0.2, δ=0.1) and 0.496 @ TX 0.525 (β=0.1, δ=0.1).
DKF's whole frontier sits well below RDKF's on nominal data.

---

## 4. Direct probe of the retuned parameters

5 LFM trajectories, T=400 — **small sample, directional only, not
publication-grade**:

| config | ssRMSE | TX |
|---|---|---|
| DKF β=0.2 δ=0.5 (shipped) | 1.043 | 0.470 |
| DKF β=1.0 δ=0.1 | 1.430 | 0.204 |
| RDKF β=0.2 δ=0.5 b=0.05 (shipped) | 1.173 | 1.000 |
| **RDKF β=1.0 δ=0.1 b=0.05** | **1.205** | **0.500** |
| RDKF β=0.6 δ=0.1 b=0.05 | 1.209 | 0.501 |
| RDKF β=2.5 δ=0.1 b=0.05 | 1.290 | 0.333 |
| RDKF β=1.0 δ=0.1 b=0 | 1.329 | 0.254 |

Retuning **halves RDKF's bandwidth (1.00 → 0.50) for a 2.7% RMSE cost**, which
makes the RDKF tradeoff curve reportable again. But at matched TX≈0.5, DKF is
still 13% better (1.043 vs 1.205).

**Retuning fixes the plot, not the result.**

---

## 5. Why RDKFLOC cannot rescue this

`extras.bLocal` in the estimate run: min 0.0481 / median 0.0483 / max 0.0494
against a global b = 0.05. Algorithm 2's local tolerances are within 4% of the
global one, so RDKFLOC ≡ RDKF (1.121 vs 1.123).

Reason: `utils/computeLocalTolerances.m` slices the joint covariance of
z = (x_{t+1}, y_t) down to the n=6 state rows plus node i's own p_i=3 measurement
rows. The model mismatch lives almost entirely in the state block, which *every*
node shares — so deleting the other 99 nodes' measurement rows barely reduces the
KL. Whether this is the paper's intent or a scaling error in the application of
eq. (26) needs checking.

Leading hypothesis for the overall RDKF deficit: each node deflates by the full b
against its *locally fused* Ω, which after the convex fusion is on the scale of a
single node's information — so the network collectively defends ~N× the intended
KL ball. That is exactly what per-node tolerances should correct, and here they
don't.

---

## 6. Item 1 resolved — the compounding deflation is faithful to the paper

Checked `estimators/RDKF.m` line by line against Ghion & Zorzi (2023)
(`papers/robust-distributed-kalman-event-triggered.pdf`).

**Verdict: not a bug.** `updateGlobalPriors` (`RDKF.m:211-233`) is a literal
transcription of eqs. (13)–(14):

    (13)  q̌ⁱ_t = cⁱ_t qⁱ_{t|t} + (1−cⁱ_t) q̄ⁱ_t
          Ω̌ⁱ_t = cⁱ_t Ωⁱ_{t|t} + (1−cⁱ_t) Ψ̄ⁱ_t      ← the already-deflated pair
    (14)  Ω̄ⁱ_{t+1} = Q⁻¹ − Q⁻¹A(AᵀQ⁻¹A + Ω̌ⁱ_t)⁻¹AᵀQ⁻¹
          find θ̄ⁱ_t s.t. γ(Ω̄ⁱ_{t+1}, θ̄ⁱ_t) = b
          Ψ̄ⁱ_{t+1} = Ω̄ⁱ_{t+1} − θ̄ⁱ_t I
          q̄ⁱ_{t+1} = Ψ̄ⁱ_{t+1} A (Ω̌ⁱ_t)⁻¹ q̌ⁱ_t

Feeding Ψ̄ back in when cⁱ=0 (line 221) and subtracting θ̄ again (line 230) is
precisely the second line of eq. (13) followed by eq. (14). **So §2.3's TX floor
is a property of Algorithm 1, not of our code, and the leading remedy in the
original checklist is off the table.**

Also verified faithful: correction (7) including the relay-node branch
Ωⁱ_{t|t} := Ψⁱ_{t|t−1}; the trigger (9); fusion (10)–(11) with the 1/(1+δ)
shrinkage; prediction (12); `updateOmega` = the Ω̄ recursion; and `findTheta`'s γ
matches the Section-2 definition exactly.

Sub-item (siblings): `SRDKF.m:249-268` has the identical branch (via
`predictNoTransmit`); `RDKFLOC.m` and `SRDKFLOC.m` define no
`updateGlobalPriors` at all — they subclass and substitute only bⁱ, which is
exactly what Algorithm 2 / eq. (28) prescribes. All four are consistent.

Two deliberate deviations from the paper's text, both benign:

- `calcFusionWeights` uses π_ij = (d_i+1)⁻¹; the paper writes (d_j+1)⁻¹. Ours is
  right — the paper calls Π a convex combination, which requires unit row sums.
  Already documented in the function header.
- `SingleTarget3dModel`'s Cⁱ split ({px,py} and {py,pz}) follows the paper's
  *prose* rather than its *formulas* (diag(1,0,1) / diag(0,1,1) = {px,pz} /
  {py,pz}). Immaterial by symmetry.

Nit: `findTheta`'s `theta_high = min(eig(Omega)) - 1e-10` uses an absolute
margin. Harmless at our λmin ≈ 2.6, but it would break on a badly scaled Ω.

### 6.1 What replaces it: an exact (b, β) feasibility bound

θ solves γ(Ω̄,θ) = b, and γ(Ω̄,·) is monotone increasing in θ. So the lower
Loewner condition Ω/(1+β) ⪯ Ω̄ − θI is *equivalent* to a ceiling on b. Taking the
best case Ω̄ = Ω_{t|t} (i.e. ignoring the extra information Ω_{t|t} carries from
fusion and from its own correction step):

    the trigger can fire at all   ⟺   b ≤ γ(Ω, λmin(Ω)·β/(1+β))        (*)

This is closed form, computable from the steady-state Ω alone — no Monte Carlo,
no sweeping. For this plant the steady-state Ω_{t|t} has eigenvalues

    [2.61  4.31  5.49 | 46.6  53.1  59.0]        cond = 22.6

a clean 3+3 split (velocity / position), so the effective dimension entering γ is
~3, not 6. Evaluating (*):

| β | 0.1 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.5 | 2.5 |
|---|---|---|---|---|---|---|---|---|
| max feasible b | 0.0037 | **0.0136** | 0.047 | 0.094 | 0.149 | 0.211 | 0.384 | 0.772 |
| λmin(I − θΩ⁻¹) there | 0.909 | 0.833 | 0.714 | 0.625 | 0.556 | 0.500 | 0.400 | 0.286 |

and inverted, the minimum β that admits a given b:

| b | 0.001 | 0.003 | 0.005 | 0.01 | 0.02 | 0.03 | 0.05 | 0.1 |
|---|---|---|---|---|---|---|---|---|
| β_min (necessary) | 0.05 | 0.09 | 0.12 | 0.17 | 0.25 | 0.31 | **0.41** | 0.62 |

So the shipped (β=0.2, b=0.05) is **3.7× over the hard feasibility budget** — the
trigger provably cannot fire, which is exactly the observed TX ≡ 1.0.

The bound is *necessary, not sufficient*: the real cutoff sits ~2–4× tighter
(measured: TX 0.542 at b=0.003, 1.0 at b=0.01) because Ω̄ is propagated from the
node's own pre-fusion information and is strictly less informative than Ω_{t|t}.

This reproduces from first principles both §2.2's "β needed ≈ 0.6" and §3.2's "β
must be ≥ 0.5", and it retires the guesswork in checklist items 2–3: the grids
can be **computed** rather than searched. (β_min at b=0.05 comes out 0.41–0.47
depending on the random Rⁱ permutation draw; both are far above 0.2.)

### 6.2 A new discrepancy with the paper

The paper's Section 6 uses **exactly our α=10, β=0.2, δ=0.5, b=0.05** on
**exactly this plant** (eq. 29, which `SingleTarget3dModel` replicates) and
reports **RDKF TX ≈ 0.38** (Fig. 3), with the robust filters *beating* DKF on
RMSE. By (*) that parameter combination is infeasible for any Ω with our
eigenvalue structure — it needs Ω to be nearly isotropic:

| Ω conditioning | θ/λmin at b=0.05 | β_min |
|---|---|---|
| isotropic (n_eff = 6) | 0.162 | 0.194 |
| synthetic n_eff = 3 | 0.219 | 0.281 |
| synthetic n_eff = 1 | 0.341 | 0.516 |
| **our actual Ω** | **0.34** | **0.41–0.47** |

β=0.2 admits b=0.05 only in the isotropic n_eff=6 case, which a constant-velocity
target with unmeasured velocity cannot produce. Our real Ω lands between the
n_eff=3 and n_eff=1 rows because the three small eigenvalues are themselves
spread (2.61 / 4.31 / 5.49).

Corroborating evidence that their Ω differs from ours: the paper's Figs. 5–6 show
steady-state θ ≈ 0.02–0.1, ours is θ ≈ 0.9. Since θ/λmin(Ω) is invariant to
rescaling Ω, their Ω must be ~10–40× smaller. Sweeping the measurement-noise
scale k (`P.noiseScale`; the paper writes Rⁱ = √k·P_i R₀ P_iᵀ but never defines
k, and we default to 1) reproduces their θ magnitudes at k ≈ 10³–10⁴:

| k | 1 | 10 | 100 | 10³ | 10⁴ |
|---|---|---|---|---|---|
| λmin(Ω) | 3.85 | 1.72 | 0.750 | 0.323 | 0.138 |
| cond(Ω) | 19.6 | 33.7 | 59.0 | 104 | 185 |
| θ at b=0.05 | 1.23 | 0.552 | 0.242 | 0.104 | 0.044 |
| **β_min** | 0.470 | 0.474 | 0.475 | 0.476 | 0.476 |

**k is a red herring for the trigger.** It reconciles the θ magnitudes with the
paper's figures but moves β_min by 1%, because θ/λmin depends only on Ω's
eigenvalue *ratios*, not its scale. So k cannot explain TX ≈ 0.38 either.

The two remaining unreplicated divergences from the paper's Section 6 are:

- **network**: 80 sensor nodes + 20 relay ("communication") nodes, random with 4%
  connection density. Ours is 100 sensors / 0 relays on a spatial radius graph
  (avg degree ≈ 4.5) — a deliberate choice recorded in
  `docs/sst3d-extension-plan.md`, which keeps the relay paths dormant.
- **initial condition**: x₀ ~ N(0, I), V_{0|−1} = I. Ours is
  x₀ = [25 25 25 1000 1000 1000]ᵀ, P₀ = 10³·I.

Neither should be able to move the steady-state Ω (the Ω/Ψ̄ recursions are
independent of the data and of x₀, and P₀ washes out of the second half of a
T=1000 horizon), so neither is expected to explain the gap — see §6.4.

### 6.3 Sensor geometry: one real bug, but not the cause

Prompted by a visibility heatmap that looked wrong, we checked C and the sensor
distribution.

**The plot is genuinely broken; the plant is not.** `plants/plotSystemChecks.m:53`
hardcodes two measurement rows per sensor:

```matlab
rows = (2*(sIdx-1)+1):(2*sIdx);   % 3D model has THREE rows per sensor
```

so in the 3D model it slices sensor 2 as rows 3:4 (straddling sensors 1 and 2),
and the local-observability heatmap is meaningless from sensor 2 onward. It is
called by all four estimate scripts (`estimateSST3dLfm.m:37`). **Diagnostics
only** — every filter derives `senBlock = plant.p / S` correctly, so no result in
this document is affected. Worth fixing regardless: it is the one thing in the
repo that will actively mislead about sensor coverage.

**But the plant has the wrong number of sensor types.** The paper says every
sensor measures "either two horizontal dimensions or a combination of one
horizontal dimension and the vertical dimension". With x,y horizontal and z
vertical, the second category has *two* realizations, so the design has **three**
types — all three coordinate pairs:

    {px,py}   two horizontal
    {px,pz}   one horizontal + vertical
    {py,pz}   one horizontal + vertical

The two matrices printed in the paper (`diag(1,0,1)`, `diag(0,1,1)`) are one
exemplar from each verbal category, not the complete set.
`plants/SingleTarget3dModel.m:61-64` implements only `{px,py}` and `{py,pz}`,
dropping `{px,pz}` — which is what produces a lopsided coverage:

| | px | py | pz | vel eig spread | β_min | max feasible b at β=0.2 | starved nodes |
|---|---|---|---|---|---|---|---|
| shipped (2 types) | 50 | **100** | 50 | 1.87× | 0.339 | 0.019 | **13 / 100** |
| 3 types (all pairs) | 67 | 67 | 66 | 1.17× | 0.305 | 0.023 | **3 / 100** |

`rank(Cⁱ) = 2` per sensor (each has one identically-zero row: row 3 for the
{px,py} type, row 1 for {py,pz}), `rank(C) = 3`. Velocity is unmeasured
throughout, as intended.

**Rⁱ, by contrast, is faithful.** `Rⁱ = √k · P_i R₀ P_iᵀ` with R₀ = 0.5·diag(1,4,7)
and a fresh random 3×3 permutation per sensor: verified that every Rⁱ is diagonal
and a permutation of R₀'s spectrum {0.5, 2, 3.5}, that the permutation is uniform
across row indices (36/28/36, 36/33/31, 28/39/33 over S=100; expected 33.3), that
D Dᵀ = R and both are block-diagonal, and that `R(idx,idx)` in the filters slices
the right block. It also adds **no asymmetry of its own**:

    sum_i Ci' Ri^-1 Ci, position diagonal = [49.9  97.4  50.9]
    predicted from coverage alone         = [46.4  92.9  46.4]
    off-diagonal position coupling        = 0 (exactly)

so the 2:1:1 imbalance is attributable entirely to the missing third C type; the
R₀ permutation averages out evenly across coordinates, as intended.

Small real bug alongside it: `SingleTarget3dModel.m:68-69` stores `self.Perm` with
the comment "stored in self.Perm so a re-constructed plant reproduces the same
Rⁱ", but the constructor never replays it — it always redraws. Reconstructing
under a different RNG state gives a different R (verified). Reproducibility rests
entirely on `rng` seeding; the comment promises something the code does not do.

**This does not explain the TX saturation.** Symmetrising coverage nearly removes
the velocity-eigenvalue spread (1.87× → 1.17×) but only moves β_min from 0.339 to
0.305, and the b-ceiling at β=0.2 from 0.019 to 0.023 — against the required
b=0.05, still **2.2× short**. So the sensor-type count is a genuine fidelity bug
but not the cause of §2's saturation. (For completeness, the *paper-literal*
2-type pairing {px,pz}/{py,pz} would be marginally worse than ours: β_min 0.373.)

(β_min lands at 0.29–0.47 across §6.1/§6.2/§6.3 depending on the fusion
approximation and the random Rⁱ permutation draw. The spread is why the bound
should be read as "β=0.2 is far out of reach", not as a precise threshold.)

**Where it does matter is the RMSE deficit.** On the shipped spatial radius graph,
fusion sets average 4.76 nodes (min 2), and **13 of 100 nodes have a fusion set
that misses a position coordinate entirely** — they can only recover it through
multi-hop consensus over many steps. Going to three types drops that to **3 of
100**, a 4× reduction, because an all-one-type neighbourhood becomes much rarer.
This is the most promising lead in this section for §1's RMSE gap, and it has no
counterpart in the paper's setup (random 4% graph, plus 20 relay nodes that mix
across the network rather than locally).

### 6.4 The network *does* matter — first result, with a caveat

Running RDKF at the shipped α=10, β=0.2, δ=0.5, b=0.05 under our configuration
versus the paper's Section 6 configuration (T=150, 1 LFM trajectory; TX is nearly
deterministic here so one run is informative):

| config | avg degree | TX (ss) |
|---|---|---|
| A — shipped: S=100, spatial graph, P₀=10³I | 3.76 | **1.0000** |
| D — paper: S=80 + 20 relays, ER 4% graph, P₀=I | 4.90 | ~~0.8403~~ discarded |

Config A reproduces §2 exactly (P(lower ok) = 0.0000, median λmin = 0.634), so the
short horizon is not distorting anything.

**Config D's number is discarded, not merely caveated.** It was computed while the
run was emitting a flood of singular-matrix warnings that went unnoticed because
stdout was redirected to `/dev/null`. A re-run with stdout captured produced
**622 MB of log containing 1.55 million `Matrix is singular to working precision`
warnings in ~11 minutes** — the warning I/O, not the filtering, was the dominant
cost, and the run never reached its output stage.

### 6.4.1 Root cause: the fusion weights are built over the wrong edge direction

The collapse is **not** caused by b, by the deflation, or by the relay nodes. It is
a latent bug in `calcFusionWeights` that a *directed* graph exposes:

- `networks/calcFusionWeights.m` builds π_ij = (d_i+1)⁻¹ over node i's
  **out**-neighbours, giving unit *row* sums.
- `RDKF.fusion` (`RDKF.m:167-186`) sums over node i's **in**-neighbours:
  `[~, nids] = inedges(self.G, i)`.

On any symmetric graph in-neighbours = out-neighbours, so the weights sum to
exactly 1 and nothing goes wrong. On a genuine digraph they do not:

| graph | Σ_j π_ij over the set `fusion()` actually sums | required |
|---|---|---|
| shipped spatial (symmetric) | min 1.0000, max 1.0000 | 1.0 |
| ER 4% digraph | min 0.1111, **median 0.2000**, max 0.5000 | 1.0 |

So the fusion step discards ~80% of the information every iteration. Ω decays
geometrically to zero; then, since A is invertible,
Q⁻¹ − Q⁻¹A(AᵀQ⁻¹A + Ω^F)⁻¹AᵀQ⁻¹ → 0 as Ω^F → 0, so Ω_pred → 0 and `findTheta`
inverts a singular matrix. Relay nodes hit it first (step 8, all 20 of them)
because they have no `CᵢᵀRᵢ⁻¹Cᵢ` term to replenish; sensors follow (97.5% of
sensor node-steps below 1e-10, with λmin going slightly negative at −6.9e-15).

The measured collapse is **independent of b**, which is what rules out the
deflation as the cause:

| b | TX (ss) | sensors: frac λmin(Ω_pred) < 1e-10 | relays collapsed |
|---|---|---|---|
| 0.05 | 0.8402 | 0.9750 | 20/20 by step 8 |
| 0.01 | 0.8400 | 0.9750 | 20/20 by step 8 |

Identical to four digits. A b-driven mechanism could not do that.

**Two things to record:**

1. **Repo bug**: `calcFusionWeights` and `fusion` disagree on edge direction. Any
   directed network silently breaks every distributed filter here — no error, just
   information leaking away. Latent because every graph the repo builds is
   symmetric (`createSpatialNetwork` builds an undirected `graph` and converts).
   The paper's Section 6 network is described as a random digraph, so this must be
   fixed before any replication attempt.
2. **My test harness was also wrong**: `paperRep2.m`/`relay2.m` built the ER graph
   from an asymmetric `rand(N,N) < 0.04`. So config D was invalid twice over.

Two earlier claims in this document were wrong and are retracted:

- That the relay-node collapse was an **A1** violation (b not small enough) — no,
  it is b-independent.
- §6.2's claim that the network "should not be able to move the steady-state Ω" —
  the reasoning was wrong (the network does change Ω through fusion-set
  composition and relay nodes), but the ER-graph evidence used to overturn it was
  itself invalid. **The question of whether the paper's network changes TX is
  therefore still open**, and needs a re-run with a correctly-weighted digraph.

---

## 7. Next steps (checklist)

Ordered by leverage. Items 2–3 are cheap but should not be run before item 1,
since item 1 can change what the right grids are.

- [x] **1. Check the compounding deflation in `updateGlobalPriors` against the
      paper** (`estimators/RDKF.m:211-233`, line 221 → 230). **Done — see §6.
      It is a faithful transcription of eqs. (13)–(14); not a bug.** The TX floor
      is a property of Algorithm 1. Replaced by the exact feasibility bound in
      §6.1 and the new paper discrepancy in §6.2.
  - [x] Cross-check the same question in `SRDKF.m` / `RDKFLOC.m` / `SRDKFLOC.m`.
        All consistent; the LOC variants inherit and only substitute bⁱ (§6).
- [ ] **2. Re-tune β and δ; drop α from the grid.** The current grids spend half
      their configs on knobs that provably never bind. Grid endpoints now come
      from (*) in §6.1 rather than guesswork.
  - [ ] `tuneRdkfBetaGrid` → log-spaced over **[0.5, 3]** at b=0.05. β_min = 0.41
        is a *necessary* bound, so anything below ~0.5 is dead weight (confirmed:
        the current grid's β=0.1 is a saturated config).
  - [ ] `tuneRdkfDeltaGrid` → [0.02, 0.3]; δ=0.1 dominates today, so resolve
        *below* it (current grid's 1 and 2.5 are pure loss).
  - [ ] Report α as measured-non-binding rather than silently fixing it — or
        sweep it down to ~0.5–2 where it would start to matter. Same applies to
        `tuneDkfAlphaFixed`.
- [ ] **3. Refine the b grid across the cliff** — `[0.003, 0.005, 0.007, 0.01]` —
      run at β≈1, δ=0.1 so the trigger stays alive across the whole sweep.
      Cross-check against the β_min table in §6.1: at β=1 every b up to 0.21 is
      feasible, so the whole grid should trigger.
  - [ ] Better: run the b sweep at **matched TX** by re-tuning β per b — and now
        the per-b β can be *initialised from (*)* instead of searched blind.
        Otherwise the RMSE-vs-b curve stays confounded with bandwidth (§3.1) and
        cannot support any claim about robustness.
  - [ ] Consider promoting (*) into a repo utility (e.g. `utils/feasibleB.m`) and
        asserting it in the tune scripts, so a config that cannot possibly
        trigger fails loudly instead of silently reporting TX = 1.
- [ ] **3b. Close the §6.2 discrepancy with the paper.** Their reported TX ≈ 0.38
      at (β=0.2, b=0.05) is infeasible under (*) for our Ω. Either their Ω is
      nearly isotropic (unlikely for this plant) or something in the trigger's Ω̄
      differs from our reading of eqs. (13)–(14).
  - [x] Replicate their Section 6 network: 80 sensors + 20 relay nodes, ER 4%
        graph. **Done (§6.4) — TX drops 1.0000 → 0.8403, so the network matters.**
        Still far from their 0.38, and not a usable operating point.
  - [ ] **Bisect §6.4's three simultaneous changes** (80/20 split vs ER topology
        vs P₀=I) to find which one carries the effect. Budget ~15 min per config.
  - [ ] **Fix `calcFusionWeights` / `fusion` edge-direction mismatch first — it
        blocks every directed-graph experiment** (§6.4.1). Weights are built over
        out-neighbours; `fusion` sums over in-neighbours. Symmetric graphs hide it;
        on a digraph the fusion weights sum to ~0.2 and all information leaks away
        silently. Either build π from in-neighbours, or sum over out-neighbours —
        whichever matches the paper's Π (it defines π_ij for j ∈ N_i, the
        in-neighbour set, so building π from in-degree is the likely fix).
    - [ ] Add a regression test asserting the fusion weights sum to 1 over exactly
          the set `fusion()` iterates, on a deliberately asymmetric digraph.
    - [ ] Make `findTheta` fail loudly (or clamp) on a singular Ω rather than
          returning a θ computed from `RCOND = NaN` arithmetic — it silently
          produced plausible-looking output through this entire failure.
    - [ ] Exercise the `N > S` relay path too; it has been dormant since the 3D
          extension (`sensorCount = nodeCount = 100`) and was never separable from
          the bug above in this run.
  - [ ] Suppress `MATLAB:singularMatrix` / `nearlySingularMatrix` in any batch
        diagnostic. Unsuppressed, one T=120 N=100 run wrote **622 MB / 1.55M
        warnings in 11 min** and never reached its output stage.
  - [x] Ruled out: the sensor-type distribution / C matrix (§6.3) and the
        measurement-noise scale k (§6.2).
  - [ ] If TX stays at 1.0 under their exact configuration, the paper's Fig. 3 is
        not reproducible from Algorithm 1 as written, which is worth stating
        explicitly in the writeup.
  - [ ] Decide what to do about `P.noiseScale` (k). k ≈ 10³ matches the paper's
        Fig. 5–6 θ magnitudes; it does not affect TX (§6.2), but it does change
        the RMSE scale and so affects comparability of every reported number.
- [ ] **3c. Plant/diagnostics fidelity fixes found in §6.3.** None of these change
      the TX story, but two are real bugs and one changes the RMSE baseline.
  - [ ] **Add the missing third sensor type** `{px,pz}` in
        `SingleTarget3dModel.m:61-64`. The paper's two printed Cⁱ are exemplars of
        two verbal categories, and "one horizontal + vertical" has two
        realizations, so there are three types. Fixing it balances coverage
        (50/100/50 → 67/67/66) and cuts nodes structurally missing a position
        coordinate from 13/100 to 3/100. **This changes every 3D number in the
        repo**, so do it before the item-5 re-run, not after.
  - [ ] Fix `plotSystemChecks.m:53` to use `senBlock = size(C,1)/sensorCount`
        instead of the hardcoded 2 — it currently produces a meaningless
        observability heatmap in every 3D run (it reads sensor 2 as rows 3:4).
  - [ ] Either replay `self.Perm` when reconstructing the plant, or drop the
        claim in the `SingleTarget3dModel.m:68-69` comment that it makes Rⁱ
        reproducible. It does not.
  - [ ] Consider reporting the single-type-fusion-set count as a network
        statistic (or asserting on it), since those nodes are structurally
        starved regardless of how many sensor types exist.
- [ ] **4. Chase the CRKF-beats-CKF / RDKF-loses-to-DKF gap** (§5). Test whether
      the per-node deflation over-defends by ~N×.
  - [ ] Confirm whether `computeLocalTolerances` is applying eq. (26) at the
        right scale, given that b^i ≈ b for all nodes here.
  - [ ] Sanity experiment: deflate against the *fused network* information scale
        rather than the local fused Ω, and see whether the deficit closes.
- [ ] **5. Once 1–4 settle, re-run the full `estimateSST3dLfm` Monte Carlo** (200
      runs, T=1000) and refresh `plotRMSEvsTXrate`. Any matched-rate claim about
      RDKF vs DKF should come from that run, not from the T=400/5-run probes in §4.

---

## 8. Reproducing the diagnostics

The instrumentation lives in the session scratchpad, not the repo:

- `RDKFDiag.m` — copy of `estimators/RDKF.m` logging per-node, per-step `eNorm`,
  `λmin(Ω⁻¹Ψ̄)`, `λmax`, `θ`, and `λmin(Ω_pred)`.
- `diag.m` — §2 table (condition pass rates vs b).
- `silence.m` — §2.3 table (λmin vs consecutive silence age).
- `probe.m` — §4 table (β/δ/b probe on LFM data).
- `anaEst.m`, `anaTune.m` — §1 and §3 tables.
- `feasib.m` — §6.1 feasibility tables. Cheap (seconds): it recomputes the
  steady-state Ω from the information recursion directly and evaluates (*); no
  filter run and no Monte Carlo needed.
- `kappa.m` — §6.2 noise-scale table.
- `cmat.m` — §6.3 sensor-coverage and C-variant tables. Cheap.
- `paperRep.m` / `paperRep2.m` — §6.4 config bisection against the paper's
  Section 6. Expensive (~15 min/config: N=100 nodes × 2 `findTheta` bisections
  per node per step); `paperRep2.m` is the trimmed retry that appends to
  `paperRep2.txt` after each config so partial progress survives.

If these need to survive, promote them to `tests/` or a `tuning/diag/` folder;
the scratchpad is session-scoped and will be lost.
