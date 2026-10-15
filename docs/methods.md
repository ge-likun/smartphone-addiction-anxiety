# Analysis definitions

Sex-grouped path models relate MPIQ, support scores and GAD-7, with survey-year adjustment in both equations. Continuous MPIQ and support scores are mean-centred. The binary sensitivity specification replaces continuous MPIQ with the >32 group indicator. Intervals use bias-corrected bootstrap resampling; the replicate count is configured in `local_config.R`.

Measurement analyses fit three correlated support factors, each with three items. The final CFA includes residual covariances for SS5–SS8, SS4–SS7, SS2–SS5, SS2–SS8 and SS1–SS7. ML complete-case, ML FIML and ordinal WLSMV fits are reported. Cronbach's alpha summarizes GAD-7, MPIQ and the support scales; composite reliability and average variance extracted summarize the support factors. Retest stability uses Pearson correlations and Fisher-z intervals.

GAMLSS models use penalized B-splines for all three support scores. Overall models adjust for sex and survey year; sex-stratified models adjust for survey year. The same predictors enter μ, σ and ν; NO has no ν equation. ZABB fixes the GAD-7 maximum at 21 and serves as a bounded-discrete sensitivity specification. Its use does not assume 21 exchangeable Bernoulli trials. Information criteria are compared within compatible likelihood sets: ZAGA/ZAIG and ZANBI/ZABB. NO is retained separately because its density at zero is not a point mass. The order-safe ZANBI implementation is checked against scalar distribution functions. Randomized quantile residuals use recorded seeds and repeated draws.

For ZAGA and ZAIG, μ is the mean of the positive-score component, ν is the probability of a zero score and σ is the distribution-specific scale parameter. Exponentiated link-scale contrasts represent mean ratios for μ, scale ratios for σ and zero-score odds ratios for ν. The cross-sectional fitted patterns describe associations.

Curve ribbons summarize post-estimation quadratic LOESS fits; they are descriptive pointwise summaries, not parameter confidence intervals.

The Gaussian longitudinal network jointly models all follow-up nodes, with all baseline nodes and the full addiction × sex × baseline-node interactions. Residual correlations are shared across groups. Regression slopes have normal(0, 0.5) priors per raw baseline unit. Coefficients are standardized after fitting using each node's SD pooled across both waves and all participants. Rows are baseline sources and columns follow-up targets. All entries, including self-loops, contribute to metrics; plot thresholds do not select metric edges.

| Metric | Definition |
|---|---|
| Sub-network Stability | Largest absolute eigenvalue of the within-domain coefficient matrix |
| Bridge Centrality | Sum of absolute incoming and outgoing cross-domain weights |
| Potential Influence | Signed outgoing sum from a support node to anxiety or MPIQ nodes |
| Cumulative Propagated Response | Signed cumulative domain response over the configured propagation steps after a +1 SD input |

Propagation uses `z[h] = t(W) %*% z[h-1]`, starting with a unit vector at the support source. Positive responses indicate higher modelled symptom values; negative responses indicate lower values. Comparisons of responses are distinct from absolute responses within a group. Spectral radii below one imply asymptotic decay for repeated application of the isolated fixed domain matrix. Network metrics are observational summaries; propagated responses are not causal estimates.

Within-sex contrasts subtract non-addicted from addicted metrics draw by draw. The overall contrast weights both sex-specific contrasts equally. The direct sex comparison subtracts the female contrast from the male contrast. Posterior means, 95% credible intervals, positive/negative probabilities and Pd are calculated; Pd is the larger sign probability. Display symbols distinguish intervals excluding zero, remaining Pd >0.90 findings, and remaining estimates. These display rules are not multiplicity-adjusted tests. Distance summaries include all pairwise Frobenius distances and MAEDs, distance differences and probabilities of largest distance.

Joint Gaussian predictive checks hold baseline values and group assignments fixed and simulate follow-up responses with the fitted residual correlations. Checks cover all nodes, marginal distributions, legal score ranges, follow-up and cross-wave correlations, group-specific residual covariance and fixed-bin residuals. Tail areas and pointwise predictive intervals are descriptive model checks. Continuous predictions are not rounded or truncated to legal categories.

Ordinal sensitivity fits one Bayesian cumulative-probit model per follow-up node. Predictors are baseline-standardized; slope priors are normal(0, 0.5) per baseline SD and thresholds have normal(0, 3) priors. This changes both likelihood and slope-prior scale. Separate targets do not estimate cross-outcome residual correlations. Expected-score derivatives with respect to raw predictors are averaged within groups, then standardized with the Gaussian model's common SD ratios. Propagation holds this average derivative matrix fixed. Category-frequency checks, comparisons of all metric contrasts, and alternative whole-target draw pairings are generated. Dependencies among sources and groups within a target are preserved when pairing independent target draws.
