# Biological Invariants Under a Singular-Set-Mediated ECL

## Computational supplement

Companion code for **Biological Invariants Under a Singular-Set-Mediated ECL: A SEQUENTION Account of Canalization, Convergence, and Operational Time**, by **Henry Arellano-Peña**, manuscript Version 6, dated 15 October 2026.

This repository contains the Python supplement for the illustrative feasibility calculations in **Section 10.3 (Test 3)**. The manuscript examines recurrent biological organization, canalization, convergence, and ecological recoverability within TCGS–SEQUENTION. Its population-genetic application compares calibrated near-neutral responses with an additional gated-response hypothesis. Within the framework, biological histories remain operationally ordered records, while the source ontology is non-temporal.

The script addresses a narrower computational question: **how many directional substitutions or reciprocal invasion trials would be required to distinguish the specified models under idealized, known-parameter assumptions?** It separates detecting a small residual bias from distinguishing that bias from a calibrated Moran prediction.

**The outputs are deterministic analytical planning calculations, not biological observations, fitted experimental results, or evidence establishing Counterspace.** The response hypothesis is an additional biological model interpreted within SEQUENTION, not a result derived from the source axioms.

## Contents

| File | Purpose |
| --- | --- |
| [`test3_feasibility.py`](test3_feasibility.py) | Exact Moran fixation probabilities, the declared quadratic-response endpoint completion, normal-approximation planning counts, and numerical identity checks. |
| [`README.md`](README.md) | Scientific scope, execution instructions, model definitions, and interpretation of outputs. |

The manuscript, LaTeX source, bibliography, and empirical datasets are not included in this code repository. The supplement does not implement all seven test families, estimate the nine candidate biological invariants, or analyze cecal anatomy or recovery data.

## Requirements and execution

**Python 3.10 or later; standard library only.** No third-party packages, network connection, input dataset, or random seed are required to run the calculations. The published script was checked with Python 3.13.5.

Clone the repository and run the built-in checks:

```bash
git clone https://github.com/hyrucanji77/Biological-Invariants-Under-a-Singular-Set-Mediated-ECL.git
cd Biological-Invariants-Under-a-Singular-Set-Mediated-ECL
python test3_feasibility.py --check
```

Use `python3` instead of `python` when that is the Python 3 executable on your system. On Windows, `py -3` can be used instead.

Print the planning table:

```bash
python test3_feasibility.py
```

Export full-precision values and assumptions as JSON:

```bash
python test3_feasibility.py --json > test3_feasibility_results.json
```

Show command-line help:

```bash
python test3_feasibility.py --help
```

Every normal command-line execution runs the built-in checks before returning results. Run without Python's `-O` or `-OO` options, because the checks use assertions. The script prints to standard output; the JSON file above is created by shell redirection.

## Model and notation

The planning scenarios specify population size $n$, log relative fitness $s$, and the dimensionless gradient ratio $y=g/a^\dagger$. Here $g$ is the magnitude of the independently calibrated biological-potential gradient, and $a^\dagger$ is the proposed response crossover scale. The script accepts $y$ as an input; it does not infer $g$, $a^\dagger$, or the potential from data.

For the declared single-founder, well-mixed Moran benchmark,

$$
u_n(s)=\frac{1-e^{-s}}{1-e^{-ns}},\qquad u_n(0)=\frac{1}{n}.
$$

The implementation uses `expm1` to evaluate the expression stably near neutrality. Define the reciprocal benchmark probabilities and their mean by

$$
\rho_E^+=u_n(s),\qquad \rho_E^-=u_n(-s),\qquad \bar\rho_E=\frac{\rho_E^++\rho_E^-}{2}.
$$

The selected unregularized mobility is

$$
\mu_0(y)=\frac{y}{\sqrt{1+y^2}}.
$$

The illustrative quadratic-response endpoint completion preserves the benchmark mean while imposing the declared reciprocal contrast:

$$
\Delta_E=\rho_E^+-\rho_E^-,\qquad
\Delta_Q=\frac{n-1}{n}s\,\mu_0(y),\qquad
\rho_Q^\pm=\bar\rho_E\pm\frac{\Delta_Q}{2}.
$$

This completion is checked for probability admissibility. It is not a fitted biological mechanism or the exact fixation formula for every possible gated birth–death process. The signed relative biases are

$$
\delta_E=\frac{\Delta_E}{2\bar\rho_E},\qquad
\delta_Q=\frac{\Delta_Q}{2\bar\rho_E}.
$$

Uppercase `Delta` denotes an absolute fixation-probability contrast; lowercase `delta` denotes a dimensionless relative directional bias. They are not interchangeable. The supplementary calculations use the Moran benchmark; the manuscript's separate Kimura comparison is not implemented in this script.

For a contrast with gap $D$, variance coefficients $v_0,v_1$, two-sided significance level $\alpha$, and nominal power $\pi$, the code evaluates the equal-allocation Gaussian planning expression

$$
N_{\mathrm{plan}}=
\left\lceil\frac{\left[z_{1-\alpha/2}\sqrt{v_0}+z_\pi\sqrt{v_1}\right]^2}{D^2}\right\rceil.
$$

Here $z_p$ is the standard-normal quantile. For invasion probabilities, the variance coefficient is $V=\rho^+(1-\rho^+)+\rho^-(1-\rho^-)$. For relative directional biases in the declared count approximation, it is $(1-\delta^2)/2$. These are planning approximations, not exact finite-sample power calculations.

## Reproduced planning values

Defaults are **two-sided $\alpha=0.05$**, **nominal power $0.80$**, and **$y=0.2$**.

| $n$ | $s$ | $y$ | $\delta_E$ | $\delta_Q$ | `m_EQ` | `R_EQ` | `R_0Q` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 0.0200 | 0.200 | 0.08976 | 0.01761 | 750 | 6,763 | 113,614 |
| 20 | 0.0100 | 0.200 | 0.09472 | 0.01858 | 673 | 12,821 | 215,390 |
| 50 | 0.0040 | 0.200 | 0.09769 | 0.01916 | 633 | 31,080 | 522,186 |
| 20 | 0.0020 | 0.200 | 0.01900 | 0.00373 | 16,823 | 319,658 | 5,370,865 |

**`m_EQ`** is the average expected substitution count **per direction** for distinguishing the benchmark and quadratic relative biases. Its JSON field is `m_E_vs_Q_per_direction`.

**`R_EQ`** is the number of independent invasion trials **per orientation** for distinguishing the benchmark and quadratic endpoint contrasts. Its JSON field is `R_E_vs_Q_per_orientation`.

**`R_0Q`** is the number of independent invasion trials **per orientation** for distinguishing a symmetric zero-bias comparator from the quadratic endpoint contrast. The comparator uses the same mean probability $\bar\rho_E$; it is not a separately evaluated neutral Moran case. Its JSON field is `R_zero_vs_Q_per_orientation`.

There are two reciprocal directions or orientations. Double a displayed count to obtain the corresponding paired total under the equal-allocation design. Substitution counts, invasion trials, fixation events, and elapsed generations are different quantities. No mutation waiting time or laboratory duration is calculated.

## Built-in validation

`--check` verifies the neutral Moran value, probability domains, the ratio identity $u_n(s)/u_n(-s)=e^{(n-1)s}$, the corresponding hyperbolic-tangent asymmetry, and the near-neutral derivative $(n-1)/(2n)$. It also checks preservation of the reciprocal probability sum by the quadratic endpoint completion, the four manuscript planning rows, and the stated Poisson signal-to-noise count convention.

Expected completion message:

```text
All identity, probability-domain, and table-reproduction checks passed.
```

These are numerical consistency and table-reproduction checks. They are not statistical tests on experimental data or a validation of the wider framework.

## Using additional scenarios

The command-line interface intentionally reproduces the four fixed manuscript scenarios. To calculate another declared scenario without changing those reference checks, import the module from the repository directory:

```python
from dataclasses import asdict
from test3_feasibility import Scenario, calculate, checks

checks()
plan = calculate(Scenario(n=20, s=0.005, y=0.2), alpha=0.05, power=0.80)
print(asdict(plan))
```

Choose $n\geq2$, finite nonzero $s$, and finite positive $y$. `calculate` rejects an inadmissible endpoint completion. Inputs must also respect the scientific domain of the comparison: the near-neutral interpretation requires small $n|s|$, and the proposed low-gradient regime requires small $y$. Passing program validation does not establish these biological assumptions or the feasibility of their calibration.

## Scope and limitations

All model parameters are treated as known. The calculations omit calibration error, shared assay uncertainty, clustering, multiple-testing corrections, censoring, and random-exposure corrections. They do not perform simulation-based power assessment, likelihood fitting, inference from longitudinal trajectories, or empirical comparison with models of reproductive variance, fluctuating selection, spatial structure, or frequency dependence.

Reciprocal introduced-mutant assays avoid waiting for spontaneous mutation and its directional arrival bias, but do not remove the information cost of rare fixation. Consequently, the displayed counts are conditional planning values, not exact minimum sample sizes or evidence of achieved experimental feasibility. A successful fit of the phenomenological response would need its own biological interpretation and would not, by itself, identify the source ontology.

## Attribution and reproducibility

Companion manuscript: **Arellano-Peña, H. (2026). _Biological Invariants Under a Singular-Set-Mediated ECL: A SEQUENTION Account of Canalization, Convergence, and Operational Time._ Version 6, 15 October 2026.**

For reproducible citation of this supplement, include the repository URL, the exact Git commit identifier used, the manuscript version, and any changed input parameters. The script in this repository is the unchanged Version 6 supplement.

SHA-256 of `test3_feasibility.py`:

```text
6f685493c528a849f4d24afaad7e3aad4aa1802d7853db7bb5bf39a1a313aaa7
```
