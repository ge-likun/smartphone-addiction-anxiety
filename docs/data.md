# Input data and scoring

The default `shared` mode reads the study's controlled-access `cross_sectional.csv` and `longitudinal.csv`. Configure their paths under `shared` in `local_config.R`. The accompanying data dictionary defines all variables. Random record identifiers link the longitudinal baseline to the cross-sectional file. GAD-7 items are already coded 0–3; MPIQ items are 1–7 and support items are 0–5.

Preparation checks scores and baseline consistency. It generates sex, addiction, symptom and grade-group labels from the source codes and scores. File row order is preserved for seeded analyses; support retest correlations use all eligible pairs. Network follow-up values and the `retest_` support values represent their respective selection rules. Retest intervals are supplied in months, without actual dates. Original identifiers or timestamps are not required. The stage writes `analysis_data.csv`, `network_data.csv` and `retest_data.csv` to `data/`.

The alternative `raw` mode starts from screened source questionnaire files. Configure source files, worksheet names and columns under `raw` in `local_config.R`. Preparation validates coding, calculates scores and links follow-up records. Source questionnaire and timestamp files are read from `data_raw/` or configured paths. Each survey input supplies a participant identifier, sex code and item responses. Demographic columns are mapped separately. Each follow-up input supplies identifiers, timestamps and the same item responses. Baseline timestamp files supply assessment times for interval summaries; an empty `baseline_times` specification leaves unavailable intervals missing.

| Measure | Raw-input coding in the example configuration | Analysis coding |
|---|---|---|
| Seven GAD-7 items | 1–4 | 0–3 for totals; 1–4 for network models |
| Eight MPIQ items | 1–7 | 1–7; total 8–56 |
| Nine support items | 1–6 | 0–5; three consecutive items per support source |
| Sex | 1 = Male; 2 = Female | Male/Female factors |

In `raw` mode, set input offsets to zero if GAD-7/support items are already zero-based. The example maps the questionnaire missing-value code `-2` to missing item responses before scoring. Configure `item_missing_codes` for the source export. Other invalid non-missing item responses are rejected and missing responses remain missing. Scores are sums of complete item sets; no imputation is performed. If reported totals are present, they are checked against the computed sums. The optional `followup_reported_totals` map can specify different follow-up columns. An empty map disables that check for follow-up totals, not item scoring. Resolve discrepancies between reported totals and item sums before analysis. MPIQ total >32 defines the addicted group. GAD-7 total ≥5 defines the descriptive symptom category. These are analysis grouping rules.

In `raw` mode, cross-sectional rows from the configured survey waves are pooled. Baseline IDs must be unique. Follow-up linkage uses normalized identifiers and the configured minimum identifier length. Integral scientific-notation IDs are expanded when conversion is safe; digit-string IDs retain their leading zeros. Unsafe conversions are reported and left unchanged. Use exports that preserve complete identifiers. Retest analyses select the earliest eligible follow-up. With the default `network_batch_rule = "earliest"`, network preparation retains the earliest submission for each ID within each batch, then selects the latest of those batch records. The alternative `"all"` selects the latest submission across all records. These rules can select different follow-up observations. Ties are resolved by source batch and row order. Assessment intervals are calculated from timestamps when available.

Participant responses require controlled access. Raw identifiers remain local to source-file preparation. In `raw` mode, ordering seeds determine the shuffled analysis-row and retest-linkage order; network rows retain sorted linkage order. Row order matters for seeded bootstrap reproduction in both input modes.

`analysis_data.csv` contains survey year, sex/sex code, demographic variables, item scores (`gad1_item`…`gad7_item`, `mpiq1_item`…`mpiq8_item`, `support1_item`…`support9_item`), totals, support-source scores and retest status. `retest_data.csv` and `network_data.csv` contain paired variables suffixed `_t1` and `_t2`, baseline sex and addiction groups, and linkage identifiers.

Analysis modules apply complete-case exclusions to their required variables and calculate the resulting sample sizes. The specified models require nonempty analysis groups, varying predictors and identifiable parameter equations.
