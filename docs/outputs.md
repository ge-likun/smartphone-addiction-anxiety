# Manuscript outputs

The analysis stages generate the following outputs under `results/`.

| Stage | Manuscript items | Generated content |
|---|---|---|
| prepare | Input for all analyses | Dataset counts, input checksums and software record |
| descriptive | Main Table 1; Supplementary Tables 1–2 | Descriptive/equivalence tables, continuous path results and retest summaries |
| measurement | Questionnaire reliability; Supplementary Tables 3–4 | GAD-7/MPIQ internal consistency, support CFA fit indices, parameters and reliability |
| binary | Binary component of Supplementary Table 5 | Binary MPIQ-group path coefficients and contrasts |
| path_tables | Supplementary Table 5 | Combined continuous and binary path table |
| gamlss | Main Figures 2–3; Supplementary Figures 2–5; Supplementary Tables 6–8 | Distribution-family models, residual diagnostics, association curves and covariate sensitivity |
| gaussian_fit | Model underlying Main Figures 4–5 and Supplementary Tables 9–12 | Joint Gaussian fit and software record |
| gaussian_results | Main Figures 4–5; Supplementary Tables 9–12 | Networks, metric contrasts, distances, posterior draws and sampling diagnostics |
| gaussian_ppc | Supplementary Note 4 | Joint response-distribution, correlation and residual checks |
| ordinal_fit | Models underlying Supplementary Tables 13–15 | Target-specific ordinal fits and sampling diagnostics |
| ordinal_results | Supplementary Tables 13–15; Supplementary Note 4 | Derivative matrices, category-frequency checks, metrics, Gaussian comparisons and draw-pairing sensitivity |

CSV outputs preserve numeric precision; publication tables use formatted intervals and probabilities. Extreme posterior probabilities are marked as Monte Carlo estimates. Figure source data are generated from the fitted models.

Supplementary Table 15 requires both Gaussian and ordinal group summaries. Run `gaussian_results` before `ordinal_results`; its numerical source CSV contains the sample sizes and all within-group support indices for both likelihoods.

Figure 1 is a study-design illustration. Supplementary Figure 1 is an AMOS CFA diagram. These manually prepared images are outside the pipeline; measurement analyses reproduce the CFA specification and fitted parameters. AMOS-specific chi-square and Bollen–Stine bootstrap results require AMOS; the R measurement outputs report lavaan estimates and fit statistics.

Within each model set and parameter, GAMLSS panels share axis limits calculated from the fitted curves and intervals. Curve and histogram values are saved as `GAMLSS_association_curve_data.csv` and `GAMLSS_support_histogram_data.csv`.

Network panels share a data-derived edge-width scale. The default raw-coefficient display threshold is 0.12 and can be configured. All edges contribute to the metrics regardless of this threshold.

Ordinal category checks use exact cumulative-probit inverse-CDF simulation by default. The supplied configuration uses `brms::posterior_predict` for the kinship target, as specified by `ordinal_ppc_engine_overrides`. Engine, seed and replicate count are recorded in the outputs.
