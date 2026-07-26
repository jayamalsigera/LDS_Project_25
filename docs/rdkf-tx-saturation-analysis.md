# RDKF vs DKF on 3D LFM data — TX saturation and RMSE deficit

Status: **diagnosed; the four fidelity bugs are now fixed (2026-07-25).** This
document records why RDKF transmits at 100% and still loses to DKF, what the
tuning sweeps do and do not tell us, and the checklist of next steps.

Update (2026-07-25, later): **checklist item 1 is resolved — see §6.** The
compounding deflation is faithful to the paper, so that remedy is off the table;
§6.1 replaces it with an exact closed-form feasibility bound that explains the
saturation and lets items 2–3 be computed instead of swept.

Update (2026-07-25, later still): **six bugs found, five fixed — see §6.5 (bugs
1–4), §6.8 (bug 6) and §6.9 (impact).** Both weight fixes leave every
symmetric-graph result bit-for-bit unchanged (verified end-to-end across all six
estimators, §6.9), but the sensor-type fix **invalidates every 3D number in
`results/`**, including the tables in §1–§4 of this document. Those tables are
retained as the record of the pre-fix state and are labelled as such. Bug 5
(sensor count) is **not** fixed — it is a shared-config conflict with the
stochastic-trigger replication; see checklist item 0.

Update (2026-07-25, final): **§6.7 finds a fifth bug — we run 100 sensor nodes,
the paper runs 20 — and with it corrected the paper's central result reproduces:
RDKF beats DKF, and the robust filter transmits more than the nominal one at the
same ratio the paper reports.** §6.6's "Fig. 3 is not reproducible" is retracted;
it was measured on the wrong configuration. §6.1's feasibility bound (\*) is also
superseded — it omitted the prediction-loss and fusion-heterogeneity terms that
actually drive the trigger. One quantity remains unexplained: their per-node Ω is
isotropic, ours has cond ≈ 22.6.

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

> **Pre-fix data.** Every number in §1–§4 comes from runs made with the 2-type C
> matrix (bug 2, fixed in §6.5), so the RMSE values are historical and will move
> when checklist item 5 re-runs the Monte Carlo. The *qualitative* findings are
> unaffected — §6.6 re-confirms TX ≡ 1.0 with the 3-type C — but no figure here
> should be reported as-is. The fusion-weight fix (bug 1) does **not** affect
> them: Π is bitwise identical on the symmetric spatial graph these runs used.

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

> Note (§6.7): this is a fact about *our shipped 100-sensor config*, not about the
> algorithm. The paper states the opposite — α binds and β, δ do not — and at the
> paper's 20-sensor network our eNorm lands on their α scale. So "α is a dead knob"
> is itself a symptom of bug 5, and the cleanest single clue in this document.

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

> Update (§6.7): this hypothesis is now **supported quantitatively**. The b that
> actually works on the paper's network is ≈0.002–0.01, i.e. 5–25× below the
> paper's 0.05 — the right order for a per-node-vs-centralized scale error. Note
> also that our bLocal (0.048–0.049) sits well above the paper's Fig. 4 values
> (0.028–0.032); with only 20 sensor nodes instead of 100, `computeLocalTolerances`
> slices a different joint covariance, so **bug 5 is a candidate explanation for
> that gap too, and re-running Algorithm 2 at S=20 is a cheap test of it.**

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

> ⚠️ **Superseded by §6.7.** The bound (\*) below is *necessary* but it measures
> only the deflation term λmin(I − θΩ⁻¹). The trigger's actual requirement is
> λmax(Ω⁻¹M) ≤ β/(1+β), which also contains the one-step prediction loss and the
> fusion heterogeneity gap. Every β_min figure in §6.1–§6.3 is therefore the wrong
> statistic, and §6.2's dismissal of k rests on it. §6.7 replaces (\*) with the
> exact Lemma-2 interval and the corrected quantity. Kept for the record.

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

**k is a red herring for the trigger.** *(Conclusion right, reasoning wrong — see
§6.7. β_min comes from (\*), which omits the prediction-loss term that k does move.
Re-tested correctly in §6.7 over a full (Q, k) grid: the true trigger quantity goes
0.7370 → 0.7799, still short of 0.8333, so k is indeed not the answer — but not for
the reason given here.)* It reconciles the θ magnitudes with the
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
should be read as "β=0.2 is far out of reach", not as a precise threshold. **All
of these β_min values are computed from (\*) and are therefore the wrong statistic
— see the banner in §6.1 and the corrected treatment in §6.7.** The *qualitative*
reading — that the third sensor type is a fidelity bug but not the cause of the
saturation — survives, and is independently confirmed by §6.7's b-sweep.)

**Where it does matter is the RMSE deficit.** On the shipped spatial radius graph,
fusion sets average 4.76 nodes (min 2), and **13 of 100 nodes have a fusion set
that misses a position coordinate entirely** — they can only recover it through
multi-hop consensus over many steps. Going to three types drops that to **3 of
100**, a 4× reduction, because an all-one-type neighbourhood becomes much rarer.
This is the most promising lead in this section for §1's RMSE gap, and it has no
counterpart in the paper's setup (random 4% graph, plus 20 relay nodes that mix
across the network rather than locally).

### 6.4 The network — first result, since retracted (superseded by §6.6)

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
  itself invalid. **Resolved in §6.6:** with the weights fixed, the network moves
  median λmin by +4.6% and TX not at all, so §6.2's *conclusion* was right even
  though its reasoning was not.

### 6.5 Bugs 1–4 are fixed

The four defects surfaced by the item-1 audit are repaired. None of them was the
cause of the TX saturation, but two were blocking further work. (Two more were
found later: **bug 5** — sensor count, §6.7, not fixed — and **bug 6** —
`calcMetropolisWeights`, §6.8, fixed. §6.9 verifies the impact of both weight
fixes on all six consuming estimators.)

| # | file | defect | fix | changes results? |
|---|---|---|---|---|
| 1 | `networks/calcFusionWeights.m` | π built over **out**-neighbours; every filter fuses over **in**-neighbours | build π from the in-degree, index `A(j,i)` | **no** — Π is bitwise identical on any symmetric graph |
| 2 | `plants/SingleTarget3dModel.m:55-78` | only 2 of the paper's 3 sensor types | three types, sensors split 34/33/33 | **yes — every 3D number** |
| 3 | `plants/plotSystemChecks.m:47-62` | hardcoded 2 measurement rows per sensor | derive `senBlock = size(C,1)/sensorCount`, assert divisibility | diagnostics only |
| 4 | `plants/SingleTarget3dModel.m:26-28, 67-71` | comment claimed `self.Perm` makes Rⁱ reproducible | comment corrected (inspection only; seed `rng` instead) | no |

**Bug 1 — fusion weights.** `calcFusionWeights` now computes
`degrees = sum(Aoff > 0, 1)'` (column sums = in-degree) and assigns
`Pi(i,j) = 1/(1+d_i)` when `A(j,i) > 0`. Verified:

| graph | Σ_j π_ij over the set `fusion()` visits — before | after |
|---|---|---|
| shipped spatial (symmetric, N=100) | 1.0000 / 1.0000 / 1.0000 | unchanged, **Π bitwise identical** |
| ER 4% digraph (N=60) | min 0.1429, median 0.3333, max 0.6667 | **1.0000 exactly, all rows** |

The bitwise-identical result on the symmetric graph is the important one: **no
existing 2D or 3D result changes because of this fix**, so it is safe to apply
retroactively. `RDKFDiag` in the scratchpad calls the same shared function, so
the diagnostics inherit the fix.

New regression test `tests/fusionWeightsUnitTest.m` asserts, for a hand-built
asymmetric digraph, an ER 4% digraph, and a shipped spatial graph, that (a) row
sums over exactly `{i} ∪ N_i` equal 1, (b) Π places no mass on nodes `fusion()`
never visits, (c) all weights are non-negative, and (d) fusing a set of
in-neighbours that all carry the same Ω returns that Ω unchanged — the convexity
property whose violation caused the geometric decay to zero. It passes, and the
old out-neighbour formula fails cases (a) and (d) as shown above.

**Bug 2 — third sensor type.** All three coordinate pairs are now generated
(`{px,py}`, `{py,pz}`, `{px,pz}`), split as evenly as possible in contiguous
blocks. At S=100 that is 34/33/33, giving position coverage **[67 67 66]**
exactly as predicted in §6.3, with `rank(Cⁱ) = 2` for every sensor, `rank(C) = 3`,
velocity columns identically zero, and `Σᵢ Cᵢ'Rᵢ⁻¹Cᵢ` position off-diagonal
exactly 0. The remaining diagonal imbalance ([67.6 64.1 53.4] for one seed) is
the random Rⁱ permutation draw, not the geometry.

**This changes every 3D number in the repo.** The tables in §1–§4 and everything
in `results/` for the 3D model predate the fix and are now historical. Re-running
the Monte Carlo (checklist item 5) is required before any 3D figure is reported.

**Bug 3 — observability heatmap.** Confirmed correct after the fix: sensor 1
reads coords {px,py}, sensor 2 reads {px,py} (previously it was sliced as rows
3:4, straddling two sensors), sensor 12 reads {px,pz}, and every sensor has
`rank(obsv(A,Cᵢ)) = 4` of 6 — two positions plus the two velocities they
integrate, which is the correct answer. The 2D path still works unchanged
(2 rows/sensor, n=4).

### 6.6 §6.4 re-run with correct fusion weights

With bug 1 fixed the directed-network experiment is unblocked, so §6.4's
discarded config D can be redone properly, and its three simultaneous changes
bisected. Every config now reports the minimum fusion-weight row sum and the
fraction of node-steps with a numerically singular Ω_pred, so a silent collapse
cannot recur unnoticed. T=150, 1 LFM trajectory, shipped α=10, β=0.2, δ=0.5,
b=0.05, and the 3-type C matrix.

| config | S / relays | graph | P₀ | avg deg | **TX (ss)** | P(lower ok) | P(upper ok) | median λmin | min row sum | frac singular Ω_pred |
|---|---|---|---|---|---|---|---|---|---|---|
| A shipped | 100 / 0 | spatial | 10³I | 3.76 | **1.0000** | 0.0000 | 1.0000 | 0.6397 | 1.0000 | 0.0000 |
| B relays only | 80 / 20 | spatial | 10³I | 3.76 | **1.0000** | 0.0000 | 1.0000 | 0.6379 | 1.0000 | 0.0000 |
| C topology only | 100 / 0 | ER 4% | 10³I | 4.90 | **1.0000** | 0.0000 | 1.0000 | 0.6693 | 1.0000 | 0.0000 |
| E init only | 100 / 0 | spatial | I | 3.76 | **1.0000** | 0.0000 | 1.0000 | 0.6397 | 1.0000 | 0.0000 |
| D paper (all three) | 80 / 20 | ER 4% | I | 4.90 | **1.0000** | 0.0000 | 1.0000 | 0.6601 | 1.0000 | 0.0000 |

**The network does not matter. §6.4's "the network *does* matter" is retracted in
full** — the entire effect was the fusion-weight bug. Under the paper's own
Section 6 configuration, correctly weighted, RDKF still transmits at **exactly
100%**, with the lower Loewner bound failing at **every** node on **every** step.

Reading the columns:

- **Ω_pred is never singular** (frac = 0.0000 in all five configs, versus 97.5%
  of sensor node-steps below 1e-10 in the buggy run), and the fusion row sums are
  1.0000 everywhere. The relay path (configs B and D, 20 relay nodes with no
  `CᵢᵀRᵢ⁻¹Cᵢ` term) now runs clean, so **the relay-node collapse was the bug, not
  the relay construction.** That closes the dormant-`N > S`-path question too.
- The ER graph's larger fusion sets do help the trigger, but only marginally:
  median λmin(Ω⁻¹Ψ̄) rises 0.6397 → 0.6693 (+4.6%) against the 0.8333 required.
  More neighbours mean a more informative Ω and hence a relatively smaller θ, but
  the effect is an order of magnitude too small — exactly what §6.1 predicts,
  since the bound depends on Ω's eigenvalue *ratios* and fusion barely changes
  them.
- P₀ is inert to four decimal places (A vs E identical), confirming it washes out.
- The 80/20 split is a slight *loss* (0.6397 → 0.6379): 20 fewer sensors means
  less total information.

**Consequence for the paper discrepancy (§6.2).** The network, the sensor/relay
split and the initial condition are all eliminated. **But the "not reproducible"
conclusion drawn here was premature — §6.7 finds a fifth real configuration bug
on our side and shows that §6.1's bound (\*) was the wrong quantity.**

### 6.7 A fifth bug (|S| = 20, not 100), and the flaw in §6.1's bound

A second, closer pass over the paper — this time reading Section 6, Figs. 4–6 and
**Lemma 2 in the Appendix** rather than only Section 3 — turned up one more real
defect on our side and invalidated the analysis that had been used to close the
question.

**The sharpest clue, which §6.2 walked past.** Section 6 says, verbatim:

> *"We also tried to increase the transmission rate by keeping fixed α = 10 and
> changing β, δ; however, we did not notice a significant growth in terms of
> transmission rate."*

and DKF2 differs from DKF1 **only** in α (10 → 0.01), moving TX 0.28 → 0.80. So
**in the paper α is the binding condition and β, δ are the dead knobs — the exact
inverse of our implementation** (§2.1: α and δ dead, β binding). Any explanation
has to flip which condition binds, not just shift a threshold.

**Bug 5: we run 100 sensor nodes; the paper runs 20.** Section 6's sentence
*"network of dimension N = 100 where there are 80 sensors and 20 sensor nodes"* is
garbled and was read as 80 sensors + 20 relays. Three independent checks say it
means **20 sensor nodes and 80 communication nodes**:

- Fig. 4 plots bⁱ "at each sensors node" and shows **~20 bars** at irregular
  indices 7, 9, 12, 13, 24, 27, 28, 30, 39, 44, 48, 53, 55, 63, 67, 70, 71, 76,
  78, 96 — mean gap 4.7, consistent with a 20-of-100 subset and impossible for an
  80-of-100 one (which would have mostly unit gaps).
- Figs. 5 and 6 split the reported θ into "across the communication nodes" and
  "across the sensor nodes", with **C := N \ S** defined explicitly in the text.
- Measured: **P(eNorm ≤ 0.01) = 0.1796 at S=20** versus 0.0029 at S=100. DKF2
  reaches TX = 0.80 with α = 0.01, which requires P(eNorm ≤ α) ≈ 0.2. Only the
  20-sensor configuration puts our eNorm on the paper's α scale.

`scripts/sst3dParams.m:17-21` sets `sensorCount = nodeCount = 100` with a comment
saying every node is kept a sensor "for now". That is a genuine fidelity bug and
it is the one thing that calibrates α.

| S (of N=100) | TX | TX sensors | TX comm | eNorm p50 | P(eNorm≤0.01) | P(lower ok) sens | P(lower ok) comm |
|---|---|---|---|---|---|---|---|
| 100 (ours) | 1.0000 | 1.0000 | — | 0.708 | 0.0029 | 0.0000 | — |
| 80 + 20 | 1.0000 | 1.0000 | 1.0000 | 0.536 | 0.0113 | 0.0000 | 0.0000 |
| 50 + 50 | 1.0000 | 1.0000 | 1.0000 | 0.185 | 0.0458 | 0.0000 | 0.0000 |
| **20 + 80 (paper)** | 1.0000 | 1.0000 | 1.0000 | 0.050 | **0.1796** | 0.0000 | 0.0000 |

So bug 5 fixes the α scale but **not** the saturation: communication nodes, which
have no `CᵢᵀRᵢ⁻¹Cᵢ` term at all, still fail the lower bound at every step with
λmin ≈ 0.65. That is the finding that kills the previous mechanism.

**§6.1's bound (\*) measured the wrong quantity.** Tracing Algorithm 1 one step:

    Ωⁱ_{t|t} = P(Ω^{i,F}_{t−1}) − θ + Mᵢ          Mᵢ = CᵢᵀRᵢ⁻¹Cᵢ,  P(X) = (A X⁻¹Aᵀ + Q)⁻¹
    Ψ̄ⁱ_t    = P(Ωⁱ_{t−1|t−1}) − θ̄

If the network were homogeneous, Ω^{i,F} = Ωⁱ, hence θ = θ̄ and **the deflation
cancels exactly**, leaving Ψ̄ = Ωⁱ_{t|t} − Mᵢ. The trigger's real requirement is
therefore

    λmin(Ω⁻¹(P(Ω) − θI)) ≥ 1/(1+β)          i.e.   λmax(Ω⁻¹M) ≤ β/(1+β)

which contains the **one-step prediction loss** and the **fusion heterogeneity
gap** — neither of which appears in (\*). (\*) evaluated only λmin(I − θΩ⁻¹). The
b-dependence enters indirectly, through how much the deflation shrinks the
steady-state Ω itself. **This is why §6.2's k-sweep "proved" k irrelevant: it
swept k through a bound that could not see the effect.**

**Lemma 2 gives the exact admissible interval.** The Appendix proves that
γ(Ω,θ) = b implies Ω − θI ⪰ μΩ with μ = λ̄⁻¹, where λ̄ > 1 solves
λ − log λ − 1 = 2b. At b = 0.05, λ̄ = 1.517 → **μ = 0.659**. The opposite extreme
is isotropic Ω, where (n/2)[t/(1−t) + log(1−t)] = b gives t\* = 0.163 →
**1 − t\* = 0.837**. So for *any* Ω in n = 6 at b = 0.05:

    0.659  ≤  λmin(I − θΩ⁻¹)  ≤  0.837

and the paper's 1/(1+β) = **0.8333 sits 0.4% below the isotropic ceiling.** β=0.2
is admissible only for an Ω that is near-perfectly isotropic *and* has no
prediction/measurement gap at all. The theoretical floor is β ≥ 0.195; the paper
chose 0.2. That margin is too tight to be accidental.

**Their Ω is isotropic; ours is not.** Figs. 5–6 report steady-state θ ≈ 0.02–0.03
for both node classes. Inverting γ(Ω,θ) = 0.05 at θ = 0.025 in the isotropic case
gives Ω ≈ 0.153·I, i.e. **cond(Ω) ≈ 1**. Our steady-state Ω is
[2.61 4.31 5.49 | 46.6 53.1 59.0], **cond = 22.6** — the velocity/position split
that constant-velocity tracking with position-only measurements necessarily
produces (velocity information is only reachable through the A⁻¹ congruence, so it
lags position by roughly (n_eff·Ts)²).

**No plant scaling closes it.** Sweeping both knobs the paper leaves ambiguous —
the process-noise scale (its text is internally inconsistent: "Wiener process with
rate of variance 0.1" at Ts = 0.1 gives Q = 0.01·I, but it prints B = √0.001·I,
i.e. Q = 0.001·I) and the measurement-noise scale k (never defined) — over
Q ∈ 10⁻⁶…10·I and k ∈ 1…10⁴, the best attainable value of the true trigger
quantity is **0.7799 against the 0.8333 required**. It is bounded away from
feasibility because θ/λmin is set by Ω's eigenvalue *ratios*, and those ratios are
fixed by A, Ts and the position-only C — not by any noise scale.

| Q = qs·10⁻³I | k=1 | k=100 | k=10⁴ |
|---|---|---|---|
| 10⁻³ | 0.7769 | 0.7789 | **0.7799** |
| 1 (shipped) | 0.7370 | 0.7618 | 0.7723 |
| 10³ | 0.3428 | 0.5729 | 0.6824 |

**Where the trigger does become feasible: b.** At the shipped plant, β = 0.2:

| b | 0.05 | 0.03 | 0.02 | **0.01** | 0.005 | 0.001 | 0 |
|---|---|---|---|---|---|---|---|
| θ | 0.449 | 0.423 | 0.398 | 0.352 | 0.301 | 0.184 | 0 |
| trigger λmin | 0.737 | 0.777 | 0.803 | **0.839** | 0.865 | 0.900 | 0.928 |
| fires? | no | no | no | **yes** | yes | yes | yes |

So b ≲ 0.015 is the crossing, reproducing §2's measured cliff (TX 0.542 at
b=0.003, 1.0 at b=0.01) and §6.1's ceiling — but now from the correct quantity.
Note the b = 0 row is 0.928, not 1.0: the residual 7% is the prediction +
heterogeneity gap, which is what leaves nominal DKF at TX ≈ 0.5 rather than lower.

**With bug 5 fixed, the paper's central claim reproduces.** Paper network (N=100,
S=20, ER 4%, P₀=I), 3-type C, paper trigger held at α=10, β=0.2, δ=0.5, sweeping
only b (T=250, 3 LFM runs):

| b | TX | ssRMSE | P(eNorm≤α) | P(lower ok) | P(upper ok) | median λmin |
|---|---|---|---|---|---|---|
| 0.05 (paper's) | 1.0000 | 1.2985 | 0.9990 | 0.0000 | 1.0000 | 0.6501 |
| 0.02 | 1.0000 | 1.2835 | 0.9993 | 0.0000 | 1.0000 | 0.7553 |
| 0.015 | 0.9602 | 1.2385 | 0.9993 | 0.0398 | 0.9999 | 0.7786 |
| 0.01 | 0.7168 | **1.0586** | 0.9987 | 0.2835 | 1.0000 | 0.7627 |
| 0.005 | 0.6390 | **1.0456** | 0.9985 | 0.3615 | 0.9998 | 0.7837 |
| **0.002** | **0.5754** | **1.0357** | 0.9984 | 0.4256 | 0.9983 | 0.8099 |
| 0 (= DKF) | 0.4416 | 1.0685 | 0.9986 | 0.5736 | 0.9708 | 0.8562 |

**RDKF now beats DKF** — 1.0357 vs 1.0685 at b = 0.002, and it wins at every
b ≤ 0.01. §3.1's "more robustness is strictly worse even when given every
measurement", called there *"the cleanest damning number in the set"*, was an
artifact of the 100-sensor configuration and the 2-type C; **it is overturned.**
Three further signatures of the paper's Section 6 now reproduce:

- **Robust transmits more than nominal**, which the paper calls out explicitly as
  the reason it had to introduce DKF2: ours 0.5754 vs 0.4416, ratio **1.30**;
  theirs 0.38 vs 0.28, ratio **1.36**.
- **RMSE magnitude**: ours 1.04–1.30, their Fig. 1 ≈ 1.2–1.5.
- **The RMSE-vs-b curve is non-monotone with an interior optimum** (1.0685 at b=0
  → 1.0357 at b=0.002 → 1.2985 at b=0.05), which is what a working robust filter
  is supposed to look like. Every previous sweep (§3.1) was monotone-worsening.

The residual offset is b, and only b: their 0.05 versus our feasible ≲0.015.

**What this leaves.** Bug 5 is real and must be fixed. But the *residual*
discrepancy is now sharply localised to a single quantity: **the paper's per-node
Ω is isotropic with λmin ≈ 0.15; ours is anisotropic with cond ≈ 22.** Everything
else has been verified line-by-line against the paper this round — eq. (9)'s
direction, Π = (dᵢ+1)⁻¹ over in-neighbours (the paper's own text says "dᵢ" and
"convex combination", and defines π_ij ≠ 0 iff (j,i) ∈ E, confirming §6.5's fix),
eqs. (12)–(14), γ, `findTheta`, `loewnerBetweenEig`, the three Cⁱ (printed
explicitly in Section 6 — confirming §6.5's fix), Q = 0.001·I from B = √0.001·I,
R₀ = 0.5·diag(1,4,7), P₀ = I, x₀ ~ N(0,I), T = 250. **The remaining question is
not "is the trigger implemented correctly" — it is "what makes their per-node
steady-state information matrix isotropic".** That is a well-posed question about
the Ω recursion, and it is the next thing to chase, not a reason to stop.

### 6.8 Does our DKF match the paper? Mostly yes — which localises the defect

**Bug 6 first: `calcMetropolisWeights` had the identical edge-direction defect.**
`DKF`, `SDKF` and `DSEACP` use `calcMetropolisWeights`, not `calcFusionWeights`,
and it built W over `A(i,j)` (out-neighbours) with `degrees = sum(A,2)`
(out-degree) while all three fuse over `inedges`. Same latent failure, same fix
(`A(j,i)`, `sum(A,1)'`). Verified bitwise identical on the symmetric spatial
graph, so no existing result moves; `tests/fusionWeightsUnitTest.m` now exercises
**both** builders against the same contract and both pass. The degree convention
(which counts the self-loop, so it is deg+1 relative to textbook Metropolis) was
deliberately left alone to keep that neutrality — worth revisiting separately.

**Our DKF is the right comparison object.** The paper's DKF1/DKF2 come from [22],
whose trigger compares against the *nominal* propagated pair Ω̄, not the deflated
Ψ̄ — the paper stresses that eq. (9) "is intrinsically different from the
transmission rule in [22]". `DKF.m` does exactly that (`exchange` tests
`loewnerBetweenEig(lower, Omega_bar, upper)`), and its fusion shrinks the stale
Ω̄ⱼ by 1/(1+δ) as [22] prescribes. So DKF ≠ RDKF-with-b=0: they differ in *which
matrix* the trigger compares against, not only in θ.

Paper network (N=100, 20 sensors + 80 comm, ER 4%, P₀=I), T=250, 5 LFM runs,
β=0.2, δ=0.5 fixed:

| filter | our TX | **paper TX** | our ssRMSE | paper Fig. 1 |
|---|---|---|---|---|
| DKF1 (α=10) | 0.3874 | **0.28** | 1.1309 | ≈1.2–1.5 |
| **DKF2 (α=0.01)** | **0.7839** | **0.80** | 1.3841 | ≈1.2–1.5 |
| RDKF (α=10, b=0.05) | 1.0000 | **0.38** | 1.3449 | ≈1.2–1.5 |
| RDKF (α=10, b=0.002) | 0.5709 | — | **1.0924** | — |
| RDKF (α=0.01, b=0.002) | 0.8688 | — | 1.1828 | — |

**DKF2 is a near-exact match: TX 0.7839 vs 0.80.** DKF1 is 0.3874 vs 0.28 — the
same ballpark and the right side of DKF2. RMSE magnitudes (1.09–1.38) sit inside
the paper's band. So the plant, the LFM, the network, C, Rⁱ, Q, the fusion, the
Metropolis weights and **the [22] trigger with its α calibration all check out
against the paper's own numbers.**

That is the most informative result in this document, because of what it excludes.
DKF and RDKF share every line except which matrix the trigger compares against —
Ω̄ (nominal) versus Ψ̄ (deflated). DKF hits the paper's numbers; RDKF saturates at
TX = 1.0. **The defect is therefore confined to the deflated comparison in
eq. (9), and specifically to the magnitude of θ̄ that γ(Ω̄,θ̄) = b = 0.05 produces
for a per-node Ω** — which is precisely the [0.659, 0.837] Lemma-2 window of §6.7
against a required 0.8333. Nothing upstream of it is implicated any more.

**One genuine anomaly, in the opposite direction from the paper.** The paper says
*"DKF2 outperforms DKF1 in the steady state because the latter is penalized by the
low transmission rate"* — more transmission, better RMSE. Ours inverts it: DKF2
transmits 2× as much (0.784 vs 0.387) and is **22% worse** (1.3841 vs 1.1309). The
same inversion shows up in the b-sweep of §6.7 and in §1's SDKF-Open/Closed pair.
More information making the estimate worse is not a tuning artifact — it points at
an over-confidence (consistency) problem in the fused pair, where the 1/(1+δ)
shrinkage applied to silent neighbours is accidentally acting as the only
regulariser. **This is now the best-defined open bug in the repo**, and it is
independent of the robust layer since DKF1 vs DKF2 differ only in α.

### 6.9 Impact of the weight fixes on every consumer — none, verified three ways

`calcMetropolisWeights` feeds **DKF, SDKF (closed and open), DSEACP**;
`calcFusionWeights` feeds **RDKF, SRDKF** and the LOC subclasses. All six iterate
`inedges(G, i)` and index `W(i,j)` / `Pi(i,j)`, so all six require row
stochasticity over the in-neighbourhood, and all six are touched by the fixes.
Three independent checks, all negative for regression:

**A — the only graph constructor in the repo is symmetric.**
`networks/createSpatialNetwork.m` is the sole graph builder; every script and
tuning file calls it, and it builds an undirected `graph` before converting to
`digraph`, so `adjacency` is symmetric and unweighted. Over 30 configurations
(N ∈ {20, 40, 100, 200} × S ∈ {N, N/2, 20} × 3 seeds), adjacency was symmetric in
every case and **both weight matrices were bitwise identical** before and after
the fix (max |ΔW| = max |ΔΠ| = 0.000e+00, `isequal` true).

**B — end-to-end, with the pre-fix matrices injected.** `W` and `Pi` are public
properties, so the old matrices can be pushed into a constructed estimator and the
actual filter outputs compared rather than just the weights. On a symmetric
N=60 graph, T=120, LFM data, with `rng` pinned before **every** run:

| estimator | TX new | ssRMSE new | TX old | ssRMSE old | ΔTX | ΔRMSE |
|---|---|---|---|---|---|---|
| DKF | 0.4500 | 1.0693 | 0.4500 | 1.0693 | 0 | 0 |
| SDKF-Closed | 0.3601 | 1.3736 | 0.3601 | 1.3736 | 0 | 0 |
| SDKF-Open | 0.0000 | 10.5875 | 0.0000 | 10.5875 | 0 | 0 |
| DSEACP | — | 1.3099 | — | 1.3099 | — | 0 |
| RDKF | 1.0000 | 1.1424 | 1.0000 | 1.1424 | 0 | 0 |
| SRDKF-Closed | 0.3664 | 1.3494 | 0.3664 | 1.3494 | 0 | 0 |

**All six bit-identical.** Pinning the RNG matters: `SDKF`/`SRDKF` use Han's
*stochastic* trigger (`rand > ν`), so an unpinned comparison shows spurious
differences of ~1% purely from RNG state — that is what an earlier version of this
test showed, and it is not a weight effect.

(Unrelated pre-existing issue visible in the table: **SDKF-Open is degenerate in
this configuration** — TX ≈ 0, ssRMSE 10.6, versus TX 0.9989 / 1.372 in §1's
N=100, T=1000 run. Identical old and new, so not caused by the fix, but the open
-loop Z = 10⁻⁶·I evidently does not survive this shorter/smaller setup. Do not
read anything into that row.)

**C — on a genuine digraph the fix is a strict repair.** Same plant, N=60 ER 8%
digraph:

| | row sum over the set `fusion()` visits | DKF | SDKF-Closed | DSEACP | RDKF |
|---|---|---|---|---|---|
| **old** | Metropolis [0.182, 0.812]; fusion [0.100, 0.500] | ssRMSE **NaN** | **NaN** | **20.77** | **18.83** |
| **new** | both exactly [1.0000, 1.0000] | 0.9775 | 1.3676 | 1.2650 | 0.9873 |

So the pre-fix code did not merely degrade on directed graphs — DKF and SDKF
returned `NaN` and DSEACP/RDKF were off by an order of magnitude, silently.

**Conclusion: zero risk to existing results, and directed networks go from broken
to working.** Nothing in `results/` moves because of either weight fix. (The
`results/` 3D files are still invalidated — by the *sensor-type* fix, §6.5, which
is a different change.) One judgement call left in place: `calcMetropolisWeights`
counts the self-loop in `degrees`, making it deg+1 relative to textbook Metropolis.
That is pre-existing, was deliberately not touched so neutrality would hold, and
is listed as a separate item below.

---

## 7. Next steps (checklist)

Ordered by leverage. Items 2–3 are cheap but should not be run before item 1,
since item 1 can change what the right grids are.

- [ ] **0c. Decide the `calcMetropolisWeights` degree convention** (§6.9). It uses
      `sum(A, ·)` on an adjacency that carries a self-loop on every node, so
      `degrees` is deg(i)+1 relative to the textbook Metropolis weight
      1/(1+max(dᵢ,dⱼ)). Left unchanged so the edge-direction fix stays bitwise
      neutral; changing it *would* move every DKF/SDKF/DSEACP result, so it needs a
      deliberate decision and a re-run, not a drive-by fix.
- [ ] **0d. `SDKF-Open` is degenerate at N=60, T=120** (TX ≈ 0, ssRMSE 10.6) while
      fine at N=100, T=1000 (§1). Pre-existing, unrelated to the weight fixes
      (bit-identical old vs new), but it means the open-loop Z = 10⁻⁶·I is not
      scale-robust and any smaller/shorter run will report nonsense for it.
- [ ] **0a. More transmission makes RMSE worse (§6.8) — the best-defined open bug.**
      DKF2 transmits 2× DKF1 (0.784 vs 0.387) and is 22% worse (1.3841 vs 1.1309);
      the paper reports the opposite ordering. DKF1 and DKF2 differ *only* in α, so
      this is independent of the robust layer, of b, and of θ. Suspect an
      over-confidence/consistency problem in the fused pair, with the 1/(1+δ)
      shrinkage on silent neighbours acting as the only regulariser.
  - [ ] Check whether Ωⁱ_{t|t} is consistent (compare the filter-reported Ωⁱ⁻¹
        against the empirical error covariance across MC runs, per node — the same
        test `tests/clsetKfUnitTest.m` already applies to CLSET-KF).
  - [ ] If it is over-confident, the circulating-information (double-counting)
        route is the first suspect even though convex fusion is supposed to
        preclude it; compare against a covariance-intersection fusion as a foil.
- [ ] **0. Fix bug 5: `scripts/sst3dParams.m:17-21` sets
      `sensorCount = nodeCount = 100`; the paper uses 20 sensor nodes + 80
      communication nodes** (§6.7). With it, RDKF beats DKF and three of the
      paper's Section 6 signatures reproduce.
  - [ ] **Not a plain error — a shared-config conflict, so do not just edit the
        constant.** `sensorCount = nodeCount` is deliberate: the stochastic-trigger
        replication (Han et al.) does not specify how a non-sensor node decides to
        transmit, so `SDKF`/`SRDKF`/`SRDKFLOC` have no defined behaviour for
        communication nodes. The 3D config currently serves two paper
        replications with incompatible network requirements. Options: (a) give the
        deterministic-trigger filters their own `sensorCount` and keep 100 for the
        stochastic ones, reporting the two families on different networks; (b) pick
        a defensible non-sensor rule for the stochastic trigger and document it as
        our own extension; (c) keep 100 everywhere and report the Ghion–Zorzi
        comparison as off-configuration. RDKF/DKF themselves handle the relay path
        correctly (§6.6: zero singular Ω_pred with 20 relays).
  - [ ] Re-derive the α grid afterwards: at S=20 our eNorm finally lands on the
        paper's α scale (P(eNorm ≤ 0.01) = 0.18 vs their ≈0.2 for DKF2), so α
        stops being a dead knob and `tuneDkfAlphaFixed` becomes meaningful.
- [ ] **0b. Close the last gap: why is their per-node Ω isotropic (λmin ≈ 0.15,
      cond ≈ 1) while ours has cond ≈ 22.6?** (§6.7.) This is now the *only*
      unexplained quantity; everything else in the trigger has been verified
      line-by-line against the paper. Ruled out already: b, network topology,
      sensor/relay split, P₀, Q scale, k, sensor types, Rⁱ, Π, γ, the Loewner
      test. Candidate directions:
  - [ ] Check whether their reported θ is averaged over *nodes and runs* in a way
        that makes the implied Ω isotropy an artifact of the averaging (Figs. 5–6
        plot θ^c and θ^s, means over C and S respectively).
  - [ ] Test whether an Ω that is *rescaled per block* (velocity vs position)
        before the trigger — e.g. a normalised/whitened comparison — would make
        β=0.2 bind the way the paper describes.
  - [ ] Derive λmin(Ω⁻¹(P(Ω)−θI)) in closed form for the CV plant to confirm
        cond(Ω) ≈ (n_eff·Ts)⁻² is forced, which would make near-isotropy
        unreachable for *any* implementation of this plant and localise the
        discrepancy to the paper's own numbers.

- [x] **1. Check the compounding deflation in `updateGlobalPriors` against the
      paper** (`estimators/RDKF.m:211-233`, line 221 → 230). **Done — see §6.
      It is a faithful transcription of eqs. (13)–(14); not a bug.** The TX floor
      is a property of Algorithm 1. Replaced by the exact feasibility bound in
      §6.1 and the new paper discrepancy in §6.2.
  - [x] Cross-check the same question in `SRDKF.m` / `RDKFLOC.m` / `SRDKFLOC.m`.
        All consistent; the LOC variants inherit and only substitute bⁱ (§6).
- [ ] **2. Re-tune β and δ; drop α from the grid.** The current grids spend half
      their configs on knobs that provably never bind. Grid endpoints now come
      from §6.7's corrected quantity (**not** (\*) in §6.1, which is superseded).
      Note §6.7 changes the priority here: at S=20 with b ≲ 0.01 the trigger already
      fires at β=0.2, so re-tuning β may be unnecessary once bug 5 is resolved.
  - [ ] `tuneRdkfBetaGrid` → log-spaced over **[0.5, 3]** at b=0.05. (The β_min=0.41
        figure that motivated this range came from (\*) and is superseded; the range
        is still right for b=0.05 because Lemma 2 caps λmin at 0.837 there, §6.7.
        At b ≲ 0.01 the useful range is instead **[0.1, 0.5]**.)
  - [ ] `tuneRdkfDeltaGrid` → [0.02, 0.3]; δ=0.1 dominates today, so resolve
        *below* it (current grid's 1 and 2.5 are pure loss).
  - [ ] Report α as measured-non-binding rather than silently fixing it — or
        sweep it down to ~0.5–2 where it would start to matter. Same applies to
        `tuneDkfAlphaFixed`.
- [ ] **3. Refine the b grid across the cliff** — `[0.003, 0.005, 0.007, 0.01]` —
      run at β≈1, δ=0.1 so the trigger stays alive across the whole sweep.
      **§6.7 has largely done this already** (b ∈ {0.05 … 0} at β=0.2 on the paper
      network), finding the crossing at b ≈ 0.015 and an RMSE optimum at b ≈ 0.002.
      What remains is to redo it at the shipped spatial network and with enough runs
      to be publication-grade.
  - [ ] Better: run the b sweep at **matched TX** by re-tuning β per b, initialised
        from §6.7's corrected quantity rather than searched blind. Otherwise the
        RMSE-vs-b curve stays confounded with bandwidth (§3.1) and cannot support
        any claim about robustness. (§6.7's curve is already non-monotone with an
        interior optimum, which is the first time that has held.)
  - [ ] Consider promoting a feasibility check into a repo utility (e.g.
        `utils/feasibleB.m`) and asserting it in the tune scripts, so a config that
        cannot possibly trigger fails loudly instead of silently reporting TX = 1.
        **Use §6.7's quantity λmax(Ω⁻¹M) ≤ β/(1+β) — or, for a cheap Ω-free
        necessary check, Lemma 2's ceiling 1 − t\*(b,n) ≥ 1/(1+β) — not (\*).**
- [ ] **3b. Close the §6.2 discrepancy with the paper.** *(Substantially advanced
      in §6.7–§6.8: bug 5 found, the paper's central claims reproduced, and the
      residual narrowed to one quantity — their per-node Ω is isotropic, ours has
      cond ≈ 22.6. See item 0b, which supersedes this. The sub-items below are
      retained for their outcomes.)*
  - [x] ~~Replicate their Section 6 network: 80 sensors + 20 relay nodes, ER 4%
        graph.~~ The first attempt (§6.4) was invalid — it hit the fusion-weight
        bug. Re-run with correct weights is in §6.6.
  - [x] **Bisected §6.4's three simultaneous changes** (80/20 split vs ER topology
        vs P₀=I). **Done — §6.6. None of them carries any effect: all five configs
        give TX = 1.0000.** §6.4's "the network does matter" is retracted; the
        whole effect was the fusion-weight bug.
  - [x] **Fixed: `calcFusionWeights` / `fusion` edge-direction mismatch** (§6.4.1,
        §6.5). π is now built from the in-degree and indexed `A(j,i)`, matching
        the set `fusion()` iterates and the paper's Π (which defines π_ij for
        j ∈ N_i, the in-neighbour set). Π is bitwise unchanged on any symmetric
        graph, so no existing result moves.
    - [x] Regression test added: `tests/fusionWeightsUnitTest.m` asserts unit row
          sums over exactly the set `fusion()` iterates, no stray mass,
          non-negativity, and fusion idempotence, on two asymmetric digraphs plus
          a shipped spatial graph. Passes; the old formula fails it.
    - [ ] Make `findTheta` fail loudly (or clamp) on a singular Ω rather than
          returning a θ computed from `RCOND = NaN` arithmetic — it silently
          produced plausible-looking output through this entire failure. Still
          open; the §6.6 harness works around it by reporting the singular-Ω
          fraction per config instead.
    - [x] Exercised the `N > S` relay path (configs B and D in §6.6, 20 relay
          nodes). Runs clean post-fix: zero singular Ω_pred. **The relay collapse
          was the fusion bug, not the relay construction.**
  - [ ] Suppress `MATLAB:singularMatrix` / `nearlySingularMatrix` in any batch
        diagnostic. Unsuppressed, one T=120 N=100 run wrote **622 MB / 1.55M
        warnings in 11 min** and never reached its output stage. (Done in the
        scratchpad harnesses; not needed in repo code.)
  - [x] Ruled out: the sensor-type distribution / C matrix (§6.3), the
        measurement-noise scale k (§6.2), Rⁱ (§6.3), and — post-fix — the network
        topology, the sensor/relay split and the initial condition (§6.6).
  - [x] ~~TX stays at 1.0 under their exact configuration, so Fig. 3 is not
        reproducible.~~ **Retracted (§6.7): the configuration was not theirs — we
        had 100 sensor nodes, they have 20. With that fixed the paper's central
        claims do reproduce and only the value of b remains off.**
  - [ ] Decide what to do about `P.noiseScale` (k). k ≈ 10³ matches the paper's
        Fig. 5–6 θ magnitudes; it does not affect TX (§6.2), but it does change
        the RMSE scale and so affects comparability of every reported number.
- [x] **3c. Plant/diagnostics fidelity fixes found in §6.3. Done — see §6.5.**
      None change the TX story, but one changes the RMSE baseline.
  - [x] **Added the missing third sensor type** `{px,pz}` in
        `SingleTarget3dModel.m`. Verified coverage 50/100/50 → **67/67/66** and
        34/33/33 sensors per type. **This changes every 3D number in the repo**,
        so the item-5 re-run is now mandatory before reporting any 3D figure.
  - [x] Fixed `plotSystemChecks.m` to derive `senBlock = size(C,1)/sensorCount`
        (with a divisibility assert) instead of the hardcoded 2. Verified on both
        the 3D (3 rows/sensor) and 2D (2 rows/sensor) plants.
  - [x] Dropped the false reproducibility claim on `self.Perm`; the comment now
        says the draws are recorded for inspection and that reproducing a given
        Rⁱ requires seeding `rng` before construction.
  - [ ] Consider reporting the single-type-fusion-set count as a network
        statistic (or asserting on it), since those nodes are structurally
        starved regardless of how many sensor types exist.
- [ ] **4. Chase the CRKF-beats-CKF / RDKF-loses-to-DKF gap** (§5). Test whether
      the per-node deflation over-defends by ~N×. **Partly resolved: §6.7 shows
      RDKF does beat DKF once bug 5 and the sensor types are corrected, at
      b ≈ 0.002–0.01. So the "over-defends by ~N×" hypothesis is consistent with
      the observed factor (b must drop ~25× below the paper's 0.05), and §5's
      reasoning stands — but the *symptom* it was written to explain is gone.**
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
  filter run and no Monte Carlo needed. **Evaluates the superseded quantity.**
- `kappa.m` — §6.2 noise-scale table. **Superseded by `qsweep.m`.**
- `cmat.m` — §6.3 sensor-coverage and C-variant tables. Cheap.
- `predloss.m` — §6.7 decomposition of the trigger into deflation vs prediction
  loss vs the true quantity, over k and the sensor fraction. Cheap (seconds);
  homogeneous-consensus fixed point, no filter run.
- `qsweep.m` — §6.7 (Q, k) feasibility grid and the b-crossing table. Cheap. This
  is the corrected replacement for `feasib.m` / `kappa.m`.
- `sensorFrac.m` — §6.7 sensor-count sweep (S ∈ {100, 80, 50, 20}) with per-group
  TX and the eNorm percentiles that calibrate α. ~4 min.
- `bfind.m` — §6.7 b-sweep on the paper network at the paper trigger. ~15 min
  (7 values × 3 runs × T=250).
- `dkfcmp.m` — §6.8 DKF1/DKF2/RDKF comparison against the paper's Fig. 1/Fig. 3
  numbers. ~5 min.
- `impact.m` / `impact2.m` — §6.9 impact analysis of the weight fixes. Part A
  (neutrality over 30 graph configs) is seconds; Parts B and C run all six
  estimators with the pre-fix matrices injected. `impact2.m` is the version with
  `rng` pinned per run, which is required for the stochastic filters.
- `verify.m` / `verify2.m` — §6.5 fix verification. Cheap.
- `metcheck.m` — §6.9 bitwise check of `calcMetropolisWeights` on a symmetric
  graph.
- `paperRep.m` / `paperRep2.m` — §6.4 config bisection against the paper's
  Section 6. Both invalid (fusion-weight bug); superseded by `paperRep3.m`.
- `paperRep3.m` — §6.6, the valid 5-config bisection. Suppresses the
  singular-matrix warnings, asserts the fusion row sums per config, reports the
  singular-Ω fraction, and appends to `paperRep3.txt` after each config so
  partial progress survives. ~2 min/config post-fix (the buggy version was slow
  because near-singular solves plus warning I/O dominated).
- `relay.m` / `relay2.m` — §6.4.1 relay-collapse measurement. Historical: the
  collapse they measured was the fusion bug and no longer reproduces.
- `verify.m` / `verify2.m` — §6.5 fix verification (old-vs-new Π row sums,
  Π bitwise equality on a symmetric graph, 3-type coverage, `plotSystemChecks`
  slicing on both the 2D and 3D plants). Cheap.

The one piece of this that is now in the repo is
`tests/fusionWeightsUnitTest.m` (§6.5). The rest is session-scoped and will be
lost; if it needs to survive, promote it to `tests/` or a `tuning/diag/` folder.
