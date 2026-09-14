required_packages(c("dplyr", "flextable", "officer"))

analysis_file <- Sys.getenv(
  "ANALYSIS_DATA_FILE",
  unset = file.path("..", "analysis_data.csv")
)
output_dir <- Sys.getenv(
  "DESCRIPTIVE_OUTPUT_DIR",
  unset = file.path("..", "results", "01_descriptive_tables")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

continuous_equivalence_bound <- as.numeric(Sys.getenv("EQUIVALENCE_BOUND_D", unset = "0.10"))
categorical_equivalence_bound <- as.numeric(Sys.getenv("EQUIVALENCE_BOUND_V", unset = "0.10"))
if (!is.finite(continuous_equivalence_bound) || continuous_equivalence_bound <= 0) {
  stop("EQUIVALENCE_BOUND_D must be positive.")
}
if (!is.finite(categorical_equivalence_bound) || categorical_equivalence_bound <= 0) {
  stop("EQUIVALENCE_BOUND_V must be positive.")
}

if (!file.exists(analysis_file)) stop("Analysis data not found: ", analysis_file)

dat <- read.csv(analysis_file, check.names = FALSE)
required <- c(
  "survey_year", "sex_code", "age_years", "bmi_kg_m2", "grade_manuscript", "mpiq_total",
  "hukou_code", "only_child_code", "addiction_group", "gad7_symptomatic_code",
  "kinship_support", "offline_friend_support", "online_friend_support", "retest_status"
)
check_columns(dat, required, "analysis_data.csv")

dat$sex <- factor(dat$sex_code, levels = c(1, 2), labels = c("Male", "Female"))
dat$grade <- factor(dat$grade_manuscript,
                    levels = c("Year 1", "Year 2", "Year 3", "Year 4", "Postgraduate"))
dat$hukou <- factor(dat$hukou_code, levels = c(1, 2), labels = c("Urban", "Rural"))
dat$only_child <- factor(dat$only_child_code, levels = c(1, 2), labels = c("Yes", "No"))
dat$addiction_group <- factor(dat$addiction_group,
                              levels = c("Non-addicted", "Addicted"))
expected_addiction_group <- ifelse(dat$mpiq_total > 32, "Addicted", "Non-addicted")
if (any(is.na(dat$addiction_group)) || any(as.character(dat$addiction_group) != expected_addiction_group)) {
  stop("addiction_group does not match the MPIQ total >32 rule.")
}
group_counts <- table(dat$addiction_group)
if (any(group_counts == 0L)) {
  stop("Both addiction groups must contain observations.")
}
dat$gad7_category <- factor(dat$gad7_symptomatic_code, levels = c(0, 1),
                            labels = c("Minimal (score <5)", "Symptomatic (score >=5)"))
dat$survey_wave <- factor(dat$survey_year, levels = c(2023, 2024), labels = c("2023", "2024"))
if (any(table(dat$survey_wave) == 0L)) {
  stop("Both the 2023 and 2024 survey waves must contain observations.")
}

grade_labels <- setNames(levels(dat$grade), levels(dat$grade))
sex_labels <- setNames(levels(dat$sex), levels(dat$sex))
hukou_labels <- setNames(levels(dat$hukou), levels(dat$hukou))
addiction_labels <- c(
  "Non-addicted" = "Non-addicted (score \u226432)",
  "Addicted" = "Addicted (score >32)"
)
gad_labels <- setNames(levels(dat$gad7_category), levels(dat$gad7_category))

demographic_specs <- list(
  list(variable = "age_years", label = "Age, years", type = "continuous"),
  list(variable = "sex", label = "Sex, n (%)", type = "categorical", levels = levels(dat$sex), level_labels = sex_labels),
  list(variable = "grade", label = "Grade, n (%)", type = "categorical", levels = levels(dat$grade), level_labels = grade_labels),
  list(variable = "hukou", label = "Hukou, n (%)", type = "categorical", levels = levels(dat$hukou), level_labels = hukou_labels),
  list(variable = "only_child", label = "Only child, n (%)", type = "categorical", single_level = "Yes"),
  list(variable = "bmi_kg_m2", label = "Body mass index", type = "continuous")
)
support_specs <- list(
  list(variable = "kinship_support", label = "Kinship support", type = "continuous"),
  list(variable = "offline_friend_support", label = "Offline-friend support", type = "continuous"),
  list(variable = "online_friend_support", label = "Online-friend support", type = "continuous")
)
gad_spec <- list(
  list(variable = "gad7_category", label = "GAD-7 category, n (%)", type = "categorical",
       levels = levels(dat$gad7_category), level_labels = gad_labels)
)

main_table <- comparison_table(
  dat, dat$addiction_group, c("Non-addicted", "Addicted"),
  c(demographic_specs, support_specs, gad_spec),
  paste0("Pooled\nn = ", format(nrow(dat), big.mark = ","))
)
names(main_table)[4:5] <- c(
  paste0("Non-addicted\nn = ", format(sum(dat$addiction_group == "Non-addicted"), big.mark = ",")),
  paste0("Addicted\nn = ", format(sum(dat$addiction_group == "Addicted"), big.mark = ","))
)

equivalence_specs <- list(
  Age = list(type = "continuous", variable = "age_years"),
  Sex = list(type = "categorical", variable = "sex"),
  Grade = list(type = "categorical", variable = "grade"),
  Hukou = list(type = "categorical", variable = "hukou"),
  `Only-child status` = list(type = "categorical", variable = "only_child"),
  `Body mass index` = list(type = "continuous", variable = "bmi_kg_m2")
)
equivalence_rows <- lapply(names(equivalence_specs), function(label) {
  spec <- equivalence_specs[[label]]
  value <- dat[[spec$variable]]
  result <- if (spec$type == "continuous") {
    continuous_tost(
      value[dat$addiction_group == "Non-addicted"],
      value[dat$addiction_group == "Addicted"],
      bound = continuous_equivalence_bound
    )
  } else {
    categorical_equivalence(
      value,
      dat$addiction_group,
      bound = categorical_equivalence_bound
    )
  }
  cbind(Characteristic = label, result)
})
equivalence_results <- do.call(rbind, equivalence_rows)
if (!all(equivalence_results$equivalent)) {
  warning("At least one demographic comparison did not meet the equivalence criterion.")
}

wave_specs <- c(
  list(list(variable = "addiction_group", label = "Smartphone-addiction group, n (%)", type = "categorical",
            levels = levels(dat$addiction_group), level_labels = addiction_labels)),
  demographic_specs,
  support_specs,
  gad_spec
)
supplementary_table_1 <- comparison_table(
  dat, dat$survey_wave, c("2023", "2024"), wave_specs,
  paste0("Pooled\nn = ", format(nrow(dat), big.mark = ","))
)
names(supplementary_table_1)[4:5] <- c(
  paste0("2023\nn = ", format(sum(dat$survey_wave == "2023"), big.mark = ",")),
  paste0("2024\nn = ", format(sum(dat$survey_wave == "2024"), big.mark = ","))
)

baseline <- dat[dat$survey_year == 2023L, , drop = FALSE]
baseline$sex <- factor(baseline$sex_code, levels = c(1, 2), labels = c("Male", "Female"))
baseline$grade <- factor(
  baseline$grade_manuscript,
  levels = c("Year 1", "Year 2", "Year 3", "Year 4", "Postgraduate")
)
baseline$hukou <- factor(baseline$hukou_code, levels = c(1, 2), labels = c("Urban", "Rural"))
baseline$only_child <- factor(baseline$only_child_code, levels = c(1, 2), labels = c("Yes", "No"))
baseline$addiction_group <- factor(baseline$addiction_group, levels = c("Non-addicted", "Addicted"))
baseline$gad7_category <- factor(
  baseline$gad7_symptomatic_code, levels = c(0, 1),
  labels = c("Minimal (score <5)", "Symptomatic (score >=5)")
)
baseline$retest_status <- factor(baseline$retest_status, levels = c("Not retested", "Retested"))
retest_counts <- table(baseline$retest_status)
if (nrow(baseline) == 0L || any(retest_counts == 0L)) {
  stop("The 2023 cohort must contain both retested and non-retested participants.")
}

baseline_demographic_specs <- list(
  list(variable = "age_years", label = "Age, years", type = "continuous"),
  list(variable = "sex", label = "Sex, n (%)", type = "categorical", levels = levels(baseline$sex), level_labels = sex_labels),
  list(variable = "grade", label = "Grade, n (%)", type = "categorical",
       levels = c("Year 1", "Year 2", "Year 3", "Year 4"), level_labels = grade_labels),
  list(variable = "hukou", label = "Hukou, n (%)", type = "categorical", levels = levels(baseline$hukou), level_labels = hukou_labels),
  list(variable = "only_child", label = "Only child, n (%)", type = "categorical", single_level = "Yes"),
  list(variable = "bmi_kg_m2", label = "Body mass index", type = "continuous")
)

supplementary_table_2 <- comparison_table(
  baseline, baseline$retest_status, c("Not retested", "Retested"),
  c(
    list(list(variable = "addiction_group", label = "Smartphone-addiction group, n (%)", type = "categorical",
              levels = levels(baseline$addiction_group), level_labels = addiction_labels)),
    baseline_demographic_specs,
    support_specs,
    gad_spec
  ),
  paste0("2023 cohort\nn = ", format(nrow(baseline), big.mark = ","))
)
names(supplementary_table_2)[4:5] <- c(
  paste0("Not retested\nn = ", format(sum(baseline$retest_status == "Not retested"), big.mark = ",")),
  paste0("Retested\nn = ", format(sum(baseline$retest_status == "Retested"), big.mark = ","))
)

write.csv(main_table[, -1], file.path(output_dir, "Main_Table_1_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(equivalence_results, file.path(output_dir, "Main_Table_1_equivalence_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(supplementary_table_1[, -1], file.path(output_dir, "Supplementary_Table_1_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(supplementary_table_2[, -1], file.path(output_dir, "Supplementary_Table_2_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

save_table_docx(
  main_table,
  "Table 1 | Participant characteristics by smartphone-addiction group",
  paste(
    "Values are mean (s.d.) or n (%). P values are from two-sided Wilcoxon rank-sum tests or Pearson chi-square tests.",
    "Continuous-variable summaries and tests use available observations.",
    "Cohen's d is calculated as non-addicted minus addicted; Cramer's V is unsigned.",
    "Equivalence tests were limited to age, sex, grade, hukou, only-child status and body mass index.",
    paste0(
      "The equivalence bounds were |d| < ", format(continuous_equivalence_bound, trim = TRUE),
      " and V < ", format(categorical_equivalence_bound, trim = TRUE), "; ",
      sum(equivalence_results$equivalent), " of ", nrow(equivalence_results),
      " demographic comparisons met these criteria."
    )
  ),
  file.path(output_dir, "Main_Table_1_publication.docx")
)
save_table_docx(
  supplementary_table_1,
  "Supplementary Table 1 | Descriptive characteristics by survey wave",
  paste(
    "Values are mean (s.d.) or n (%). P values are from two-sided Wilcoxon rank-sum tests or Pearson chi-square tests.",
    "Continuous-variable summaries and tests use available observations.",
    "Cohen's d is calculated as 2023 minus 2024; Cramer's V is unsigned."
  ),
  file.path(output_dir, "Supplementary_Table_1_publication.docx")
)
save_table_docx(
  supplementary_table_2,
  "Supplementary Table 2 | Baseline characteristics of participants with and without linked follow-up data",
  paste(
    "The comparison is restricted to participants surveyed in 2023. Values are mean (s.d.) or n (%).",
    "P values are from two-sided Wilcoxon rank-sum tests or Pearson chi-square tests.",
    "Continuous-variable summaries and tests use available observations.",
    "Cohen's d is calculated as not retested minus retested; Cramer's V is unsigned."
  ),
  file.path(output_dir, "Supplementary_Table_2_publication.docx")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Descriptive tables written to ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
