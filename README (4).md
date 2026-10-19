# Pacejka Magic Formula Tire Model — Calibration and Vehicle-Data Evaluation

**MATLAB implementation of longitudinal and lateral tire-force models under combined-slip conditions, using genetic-algorithm optimization, local refinement, and comparison with tire-test and Formula Student vehicle data.**

This project was developed in a Formula Student vehicle-dynamics context. It explores both the calibration of tire-force models using test-rig data and the evaluation/adaptation of longitudinal-force predictions using vehicle acceleration measurements.

## Project overview

The repository contains four MATLAB scripts covering two complementary parts of tire modeling:

1. **Tire-test data fitting:** calibrate longitudinal ($F_x$) and lateral ($F_y$) Magic Formula models, including combined-slip reduction factors.
2. **Vehicle-data evaluation:** estimate wheel normal loads, wheel slip ratios, and longitudinal tire forces from measured vehicle signals; then compare or fit Magic Formula curves for selected load ranges.

The focus is on **nonlinear mathematical modeling, parameter estimation, optimization, and model-to-data comparison**.

## Repository structure

```text
pacejka-magic-formula-tire-model/
├── README.md
├── Fx_combined.m
├── Fy_combined.m
├── dynamis_fz.m
└── dynamis_reverse_eng.m
```

### MATLAB scripts

| File | Purpose |
| --- | --- |
| [`Fx_combined.m`](Fx_combined.m) | Combines three tire-test runs, groups observations by test conditions, interpolates longitudinal-force data, and fits a **19-parameter** pure/combined-slip $F_x$ model using a genetic algorithm followed by `fmincon`. Includes a sample comparison of measured data, combined-slip predictions, and a reference pure-slip model. |
| [`Fy_combined.m`](Fy_combined.m) | Combines four tire-test runs and fits a **23-parameter pure lateral-force model**; where nonzero longitudinal slip is available, fits an additional **five-parameter** combined-slip factor. Uses weighted fitting and plots measured versus predicted $F_y$. |
| [`dynamis_fz.m`](dynamis_fz.m) | Uses a Formula Student acceleration run to estimate wheel loads, slip ratios, and longitudinal forces; constructs $F_x$–slip maps by normal-load bin; extracts simplified Magic Formula parameters and compares the resulting curves with a reference 15-parameter tire-rig model, including an attempted parameter rescaling. |
| [`dynamis_reverse_eng.m`](dynamis_reverse_eng.m) | Uses the same kind of vehicle-data processing and performs **GA + `fmincon`** fitting of a simplified Magic Formula to estimated rear-right-wheel $F_x$–slip data, with a fixed shape factor $C=1.46$. |

The two `dynamis` scripts are **alternative vehicle-data analysis approaches**, not stages that must be executed in sequence.

## Mathematical model

The core Pacejka-type Magic Formula used in the scripts has the form

$$
F(x)=D\sin\!\left[C\arctan\!\left(Bx-E\left(Bx-\arctan(Bx)\right)\right)\right],
$$

where $x$ denotes the relevant slip variable. The coefficients $B$, $C$, $D$, and $E$ influence the stiffness, shape, peak amplitude, and curvature of the force–slip relationship. In the extended models, these quantities or their underlying coefficients also depend on operating conditions, such as normal load and inflation pressure.

For combined-slip operation, the MATLAB implementations use multiplicative reduction factors:

$$
F_x(\kappa,\alpha)=F_{x0}(\kappa)\,G_{x\alpha}(\alpha),
$$

$$
F_y(\alpha,\kappa)=F_{y0}(\alpha)\,G_{y\kappa}(\kappa).
$$

Here, $\kappa$ is longitudinal slip ratio, $\alpha$ is slip angle, and $F_{x0}$ and $F_{y0}$ are the respective pure-slip predictions. The reduction factors are implemented with cosine/arctangent expressions whose coefficients are fitted from data. **These expressions summarize the model structure; the MATLAB functions contain the actual load-, pressure-, and parameter-dependent formulations.**

## Calibration workflow

### Longitudinal force — `Fx_combined.m`

- Loads `B2356run71.mat`, `B2356run72.mat`, and `B2356run73.mat`.
- Groups observations by **normal load, inflation pressure, inclination angle, and test speed**.
- Aggregates repeated slip-ratio values and constructs smooth $F_x$ data using shape-preserving cubic interpolation (`pchip`).
- Uses a **19-parameter** model: 15 pure-longitudinal coefficients plus four combined-slip coefficients.
- Minimizes force-prediction root-mean-square error (RMSE) with a **genetic algorithm (`ga`)**, then refines the parameter vector using **`fmincon`**.
- Generates a sample overlay of measured $F_x$, combined-slip predictions, and reference pure-slip predictions.

### Lateral force — `Fy_combined.m`

- Loads `B2356run2.mat`, `B2356run4.mat`, `B2356run5.mat`, and `B2356run6.mat`.
- Groups observations by the same test-condition categories, and fits selected **raw** lateral-force data.
- Filters to observations with negative slip angle and positive lateral force; the fitting subset also uses zero-inclination-angle categories.
- Uses a weighted fitting objective, reducing the influence of large deviations and emphasizing small slip angles.
- Optimizes the **23-parameter** pure-slip model using `ga` followed by `fmincon`.
- If the retained dataset includes meaningful longitudinal-slip variation, additionally optimizes a **28-parameter** combined-slip model (23 pure-slip + 5 reduction-factor coefficients).

### Vehicle-data analysis — `dynamis_fz.m` and `dynamis_reverse_eng.m`

Both scripts load `DynamisPRC_Acc (1).mat` and process signals such as vehicle speed, longitudinal acceleration, and individual wheel/motor speeds. They:

1. Compute approximate **static and dynamic normal-load distributions**, including longitudinal load transfer and aerodynamic forces.
2. Estimate wheel slip ratios from wheel angular speeds and longitudinal vehicle speed.
3. Infer total longitudinal force from vehicle acceleration, aerodynamic drag, and rolling resistance, and distribute it among wheels in proportion to estimated normal load.
4. Filter the rear-right-wheel data and form **50 N normal-load bins** (nominally $\pm25$ N).
5. Construct force–slip plots and compare or fit simplified Magic Formula representations.

The methods differ in their final modeling step:

- **`dynamis_fz.m`** derives approximate $B,C,D,E$ coefficients from force–slip characteristics, compares these curves with reference tire-rig predictions, and explores adjusting the rig-model coefficients. It reports curve-comparison RMSE/MAE where available.
- **`dynamis_reverse_eng.m`** holds $C=1.46$ fixed and optimizes $B,D,E$ with a weighted error criterion. Higher weights are assigned near the measured force peak.

**Interpretation:** The forces in these two scripts are **estimated from vehicle dynamics**, not obtained from individual tire-force sensors. Their comparison with Magic Formula curves should therefore be interpreted in light of the vehicle-model assumptions. These scripts also use an equal left/right normal-load split and do not explicitly model lateral load transfer.

## Data requirements and availability

**The raw `.mat` datasets are not included in this repository.** The scripts require the following input files:

| Script | Required MATLAB data files |
| --- | --- |
| `Fx_combined.m` | `B2356run71.mat`, `B2356run72.mat`, `B2356run73.mat` |
| `Fy_combined.m` | `B2356run2.mat`, `B2356run4.mat`, `B2356run5.mat`, `B2356run6.mat` |
| Both `dynamis` scripts | `DynamisPRC_Acc (1).mat` |

The tire-test scripts expect arrays named `FZ`, `V`, `IA`, `P`, `FX`, `FY`, `SA`, and `SL` or `SR` (depending on the script). The vehicle-data scripts expect a `run01` structure with compatible vehicle-speed, acceleration, and wheel-speed fields.

**Because the source datasets are not provided, the scripts cannot be reproduced end-to-end from the repository alone.** The repository documents the implementation and analysis workflow. Any future release of raw test or telemetry data should first be checked for the relevant ownership and sharing permissions.

## Requirements

- **MATLAB** with support for script-local functions and table operations.
- **Global Optimization Toolbox** for `ga`.
- **Optimization Toolbox** for `fmincon` and `optimoptions`.
- **Statistics and Machine Learning Toolbox** for the `iqr` calculation used in `Fy_combined.m` (depending on MATLAB release and installed functions).
- **Parallel Computing Toolbox** if the `UseParallel = true` genetic-algorithm settings are retained and parallel execution is desired. Change these settings to `false` if parallel execution is unavailable.

The code has **not been re-executed in this repository environment**, because the required datasets and MATLAB runtime were not supplied.

## Running the scripts

1. Open MATLAB and set the **Current Folder** to the repository root.
2. Place the required `.mat` files in that folder, or update the `load(...)` paths in the scripts.
3. Execute the relevant workflow individually:

```matlab
run('Fx_combined.m')
% or
run('Fy_combined.m')
% or
run('dynamis_fz.m')
% or
run('dynamis_reverse_eng.m')
```

**Important: do not run the two tire-test fitting scripts back-to-back without saving the results you need.** Both `Fx_combined.m` and `Fy_combined.m` attempt to delete and recreate the folder named `Data Plots and Mats` at startup. Running one can therefore remove outputs generated by the other.

Generated plots and `.mat` files from the test-rig fitting scripts are saved under `Data Plots and Mats/`. The scripts also display some comparison figures interactively. The vehicle-data scripts primarily generate figures and printed results, with an additional `.mat` curve export in `dynamis_fz.m`.

## Suggested figures for the repository

If permission is available to publish the relevant outputs, the following figures provide a concise visual overview:

- **Longitudinal-force calibration:** measured $F_x$ versus slip ratio, overlaid with combined- and pure-slip predictions (`Fx_combined.m`).
- **Lateral-force calibration:** measured $F_y$ versus slip angle, overlaid with the combined-slip prediction (`Fy_combined.m`).
- **Vehicle-data evaluation:** rear-right-wheel $F_x$ versus slip ratio for selected normal-load bins (`dynamis_fz.m` or `dynamis_reverse_eng.m`).
- **Model adaptation:** comparison of the vehicle-derived Magic Formula curve and the original/adjusted tire-rig prediction (`dynamis_fz.m`).

Figures can be added later in a `results/` folder and linked from this README. No numerical accuracy claims are made here because the fitting results and input datasets were not supplied with the code.

## Current limitations

- **Code check before execution:** the supplied `dynamis_fz.m` appears to contain an unmatched `end` just before the final error-summary table (near line 413). Its block structure should be reviewed and corrected before running.
- In `dynamis_fz.m`, the generated file `dyn_MF_curve_Fz850_pm25.mat` has a fixed name even when the script processes several different load bins; later iterations can overwrite earlier exported curves.
- The tire-test comparison plots are based on subsets of data used in the fitting workflow; they should **not** be described as independent, held-out validation without a separate evaluation split.
- The vehicle-data force estimates depend on simplified assumptions for load transfer, aerodynamic forces, rolling resistance, and tire-force allocation.
- `dynamis_reverse_eng.m` fits one selected rear-right-wheel load range (set to approximately **850 N** in the script) rather than automatically validating the model across all vehicle conditions.
- Some run-specific paths, load bins, reference coefficients, and plotting choices are hard-coded and may need adjustment for other datasets.
- The MATLAB scripts have not been packaged into an automated end-to-end pipeline or tested here without the original data.

## Context

This work illustrates how **physics-based tire models** can be combined with **data processing, global/local optimization, and experimental or vehicle-operational measurements** in a Formula Student engineering workflow.
