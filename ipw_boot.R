####################################################################################################
# Description: This contains steps of the IPW analysis that are called by ipw_run.R
#
# The most important exported functions are: 
#   ipw_est()  <- creates point estimates 
#   ipw_boot()   <- creates bootstrapped estimates that can be used to construct a 95% CI
#   summarize_ipw_est()  <- creates summary statistics for the point estimates
#   summarize_ipw_boot()  <- creates summary statistics for the bootrapped estimates
#
# CALLED FROM - ipw_run.R
#
# DEPENDENT SCRIPTS - None
#
# INPUT FILES - None
#
# OUTPUT FILES - None
#
# Authors: Theodore Pak, Anna Agan, Jessica G. Young
####################################################################################################

# Load packages
suppressMessages({
    require(plyr)
    require(ggplot2)
    require(dplyr)
    require(splines)
    require(haven)
    require(tidyr)
    require(data.table)
    require(purrr)

    require(progressr)
})

# utility function for appending text to a file (and creating it if it doesn't exist)
cat_append <- function(text_to_append, file_path) {
    if (!is.character(text_to_append)) { text_to_append <- capture.output(text_to_append) }
    cat(text_to_append, file = file_path, append = TRUE, sep = "\n")
}

# utility function that gets weighted estimates of hazards without MSM
non_parametric_cum_haz <- function(weightVector, inputdata, cutTimes, outcomeEvent = TRUE) {
    outputHazards <- rep(NA, length.out = length(cutTimes))
    counter <- 1
    for (i in cutTimes) {
        if (outcomeEvent) {
            indices <- inputdata$t0 == i & inputdata$D == 0
            eventIndicator <- indices & inputdata$Y == 1
        } else {
            indices <- inputdata$t0 == i
            eventIndicator <- indices & inputdata$D == 1
        }
        outputHazards[counter] <- sum(weightVector[eventIndicator]) / sum(weightVector[indices])
        counter <- counter + 1
    }
    return(outputHazards)
}

# utility function that gets weighted estimates of cumulative incidence from nonparametric weighted estimates
# of hazards without MSM
non_parametric_cum_inc <- function(hazard1, hazard2, cutTimes, competing = FALSE) {
    inc <- rep(NA, length.out = length(cutTimes))
    cumulativeSurvival <- c(1, cumprod((1 - hazard1) * (1 - hazard2)))
    counter <- 1
    for (i in seq_along(cutTimes)) {
        if (!competing) {
            inc[i] <- hazard1[i] * (1 - hazard2[i]) * cumulativeSurvival[i]
        } else {
            inc[i] <- hazard2[i] * cumulativeSurvival[i]
        }
    }
    cumInc <- cumsum(inc)
    return(cumInc)
}

# Debug function for counting patients after `ipw_est()` has run and generated weights.
# These counts provide several sanity checks that weighting and input data are as expected
debug_counts <- function(dt, out_debug) {
    cat_append("\n\n#### DEBUG PATIENT EPISODE COUNTS ####\n", out_debug)
    cat_append("Distribution of t0 values per ID\n", out_debug)
    cat_append(summary(dt[, .(uniqt0 = uniqueN(t0)), by=id][, uniqt0]), out_debug)
    
    discharged_before_6h <- length(unique(dt[discharge == 1 & t0 < 11, id]))
    cat_append(paste0("\n1) Discharged before 6h (0 expected): ", discharged_before_6h), out_debug)
    death_before_6h <- length(unique(dt[death == 1 & t0 < 11, id]))
    cat_append(paste0("\n2a) Death before 6h (0 expected): ", death_before_6h), out_debug)
    # Filter out all of the patients in (1) and (2a), since they should have been excluded!
    if (discharged_before_6h + death_before_6h > 0) {
        cat_append("\nWARN: Excluding episodes meeting 1 or 2a from subsequent counts!", out_debug)
    }
    bad_ids <- unique(dt[(discharge == 1 | death == 1) & t0 < 11, id])
    dt.filt <- dt[!(id %in% bad_ids), ]
    
    deaths <- length(unique(dt.filt[death == 1 & t0 == 11, id]))
    cat_append(paste0("\n2b) Deaths at or beyond 6h (non-zero expected): ", deaths), out_debug)
    N <- nrow(dt.filt[t0 == 0, ])
    cat_append(paste0("\n3) Total number of episodes N (non-zero expected): ", N), out_debug)
    N_alt <- length(unique(dt.filt[, id]))
    cat_append(paste0("\n4) Total number of distinct IDs (should equal N): ", N_alt), out_debug)
    ids_0_3 <- setdiff(unique(dt.filt[, id]), unique(dt.filt[weights_0_3 == 0, id]))
    N_0_3 <- length(ids_0_3)
    cat_append(paste0("\n5a) Number of episodes with nonzero weights_0_3 for all t0 (N_0_3): ", 
        N_0_3), out_debug)
    deaths_0_3 <- length(unique(dt.filt[id %in% ids_0_3 & death == 1, id]))
    cat_append(paste0("\n5b) Number of deaths with nonzero weights_0_3 for all t0 (deaths_0_3): ", 
        deaths_0_3), out_debug)
    ids_3_6 <- setdiff(unique(dt.filt[, id]), unique(dt.filt[weights_3_6 == 0, id]))
    N_3_6 <- length(ids_3_6)
    cat_append(paste0("\n5c) Number of episodes with nonzero weights_3_6 for all t0 (N_3_6): ", 
        N_3_6), out_debug)
    deaths_3_6 <- length(unique(dt.filt[id %in% ids_3_6 & death == 1, id]))
    cat_append(paste0("\n5d) Number of deaths with nonzero weights_3_6 for all t0 (deaths_3_6): ", 
        deaths_3_6), out_debug)
    cat_append(paste0("\n6) N - N_0_3 - N_3_6 (0 expected): ", N - N_0_3 - N_3_6), out_debug)
    cat_append(paste0("\n7) deaths - deaths_0_3 - deaths_3_6 (0 expected): ",
        deaths - deaths_0_3 - deaths_3_6), out_debug)
}

# Weights for intervention occuring during or between intervals `markl` to `marku`
#    e.g., markl=0 and markl=1 means within the first hour, for half-hour intervals
int_weights <- function(dt, markl, marku, id_col = "id") {
    # Make a limited copy of `dt`, with only the columns we need to calculate the weights
    cols <- c("t0", "cumsumA", "p0denom", "p1denom", id_col)
    dt.temp <- dt[, ..cols]
    dt.temp[, wgt_temp := 9999]
    dt.temp[, intcheckl := (t0 - markl) + 1]
    dt.temp[, intchecku := (t0 - marku) + 1]

    # If we are in any half hour prior to when a person should start under the intervention then weight 
    # contribution is 1 over chance of not starting by that time given covariates for people who have not 
    # yet started by that time.   
    # O.w. weight contribution is 0.
    if (markl > 0) {
        dt.temp[t0 < markl, wgt_temp := 0]
        dt.temp[t0 < markl & cumsumA == 0, wgt_temp := 1 / p0denom]
    }

    # During the period where treatment is allowed but doesn't have to start (grace period), the weight 
    # contribution is 1.
    # Whatever the person is doing in that period is allowed under intervention
    dt.temp[t0 >= markl & t0 < marku, wgt_temp := 1]

    # In the half hour constituting the end of the intervention period, for an indiviual who just starts 
    # at that time, the weight contribution is 1 over the chance of just starting at that time given covaiates
    # for an individual who started earlier during the grace period, the weight contribution is 1
    # O.w. weight contribution is zero
    dt.temp[t0 == marku, wgt_temp := 0]
    dt.temp[t0 == marku & cumsumA == 1, wgt_temp := 1 / p1denom]
    dt.temp[t0 == marku & cumsumA > 1 & cumsumA <= intcheckl, wgt_temp := 1]

    # In the half hour constituting any time after the intervention period for an individual who started 
    # during the grace period or at the end of the intervention period weight contribution is 1.
    # O.w. weight contribution is zero
    dt.temp[t0 > marku, wgt_temp := 0]
    dt.temp[t0 > marku & cumsumA >= intchecku & cumsumA <= intcheckl, wgt_temp := 1]

    # Calculate final inverse probability weights for complete data
    # It's important that you specify the grouping column for this step in `id_col`
    dt.stabwts <- dt.temp[, .(stabwts = cumprod(wgt_temp)), by = id_col]

    return(dt.stabwts$stabwts)
}

# Fit a logistic regression model that is used to generate the IPW weights
ipw_glm_denom_prob <- function(dt, Kint, adj = "max") {
    form <- if (adj == "max") {
        (A ~ t0_cat
            # ----- Time-invariant covariates -----
            + hospitalType + agecat + raceEthnicity + sex + elix_AIDS + elix_ALCOHOL + elix_ANEMDEF 
            + elix_AUTOIMMUNE + elix_BLDLOSS + elix_CANCER_LEUK + elix_CANCER_LYMPH+ elix_CANCER_METS 
            + elix_CANCER_SOLID + elix_COAG + elix_DEMENTIA + elix_DEPRESS + elix_DIAB_CX +  elix_DIAB_UNCX 
            + elix_DRUG_ABUSE + elix_HF + elix_HTN_CX + elix_HTN_UNCX + elix_LIVER_MLD + elix_LIVER_SEV 
            + elix_LUNG_CHRONIC + elix_NEURO_MOVT + elix_NEURO_OTH + elix_NEURO_SEIZ + elix_OBESE 
            + elix_PARALYSIS + elix_PERIVASC + elix_PSYCHOSES + elix_PULMCIRC + elix_RENLFL_MOD + elix_RENLFL_SEV
            + elix_THYROID_HYPO + elix_THYROID_OTH + elix_ULCER_PEPTIC + elix_VALVE + elix_WGHTLOSS + elix_CBVD 
            + elixcat + dischargePrior90Days + insType + facility_admission + firstCodeStatusGrp 
            + ems + intubatedPreHosp
            ### AHRQ infection categories were added to the max adjusted model (5/8/23)
            + ahrq_infcat_pulmonary + ahrq_infcat_cns + ahrq_infcat_intra_abdominal + ahrq_infcat_genitourinary
            + ahrq_infcat_bone_joint + ahrq_infcat_septicemia_bacteremia + ahrq_infcat_skin_soft_tissue_infection
            + ahrq_infcat_other + ahrq_infcat_obstetric_gynecologic
            ### Year, beta-lactam allergy, and spoken language of English (6/30/23)
            + edArrivalYR + allergy_bl + primLangEnglish
            # ----- Time-varying covariates -----
            # Vitals including O2
            + tempMax_cat + HR_cat + SBP_cat + RR_cat + GCS_cat + o2f
            # Labs
            + creat_cat  + plt_cat + lac_cat + hct_cat + wbc_cat + tbili_cat + ag_cat + ast_cat
            + sdm_cat + alb_cat + glu_cat +
            # IVF boluses (added 6/30/23)
            + fluid_bolus_cumsum_cat)
    } else if (adj == "max_HCA") {
        (A ~ t0_cat
            # ----- Time-invariant covariates -----
            ### HCA does not have firstCodeStatusGrp
            + hospitalType + agecat + raceEthnicity + sex + elix_AIDS + elix_ALCOHOL + elix_ANEMDEF 
            + elix_AUTOIMMUNE + elix_BLDLOSS + elix_CANCER_LEUK + elix_CANCER_LYMPH+ elix_CANCER_METS 
            + elix_CANCER_SOLID + elix_COAG + elix_DEMENTIA + elix_DEPRESS + elix_DIAB_CX +  elix_DIAB_UNCX 
            + elix_DRUG_ABUSE + elix_HF + elix_HTN_CX + elix_HTN_UNCX + elix_LIVER_MLD + elix_LIVER_SEV 
            + elix_LUNG_CHRONIC + elix_NEURO_MOVT + elix_NEURO_OTH + elix_NEURO_SEIZ + elix_OBESE 
            + elix_PARALYSIS + elix_PERIVASC + elix_PSYCHOSES + elix_PULMCIRC + elix_RENLFL_MOD + elix_RENLFL_SEV
            + elix_THYROID_HYPO + elix_THYROID_OTH + elix_ULCER_PEPTIC + elix_VALVE + elix_WGHTLOSS + elix_CBVD 
            + elixcat + dischargePrior90Days + insType + facility_admission + ems + intubatedPreHosp
            ### AHRQ infection categories are named slightly differently in HCA
            + ahrq_pulmonary + ahrq_cns + ahrq_intra_abdominal + ahrq_genitourinary
            + ahrq_bone_joint + ahrq_septicemia + ahrq_skin_tissue
            + ahrq_other + ahrq_obs_gyno
            ### HCA does not have spoken language
            + edArrivalYR + allergy_bl
            # ----- Time-varying covariates -----
            # Vitals including O2
            + tempMax_cat + HR_cat + SBP_cat + RR_cat + GCS_cat + o2f
            # Labs
            + creat_cat  + plt_cat + lac_cat + hct_cat + wbc_cat + tbili_cat + ag_cat + ast_cat
            + sdm_cat + alb_cat + glu_cat
            # NOTE: Cannot do IVF boluses in the HCA dataset, yet
            )
    } else if (adj == "baseline") {
        (A ~ t0_cat
            # ----- Time-invariant covariates only -----
            + hospitalType + agecat + raceEthnicity + sex + elix_AIDS + elix_ALCOHOL + elix_ANEMDEF 
            + elix_AUTOIMMUNE + elix_BLDLOSS + elix_CANCER_LEUK + elix_CANCER_LYMPH+ elix_CANCER_METS 
            + elix_CANCER_SOLID + elix_COAG + elix_DEMENTIA + elix_DEPRESS + elix_DIAB_CX +  elix_DIAB_UNCX 
            + elix_DRUG_ABUSE+ elix_HF + elix_HTN_CX + elix_HTN_UNCX + elix_LIVER_MLD + elix_LIVER_SEV 
            + elix_LUNG_CHRONIC + elix_NEURO_MOVT + elix_NEURO_OTH + elix_NEURO_SEIZ + elix_OBESE 
            + elix_PARALYSIS + elix_PERIVASC + elix_PSYCHOSES + elix_PULMCIRC + elix_RENLFL_MOD + elix_RENLFL_SEV
            + elix_THYROID_HYPO + elix_THYROID_OTH + elix_ULCER_PEPTIC + elix_VALVE + elix_WGHTLOSS + elix_CBVD 
            + elixcat + dischargePrior90Days + insType + facility_admission + firstCodeStatusGrp 
            + ems + intubatedPreHosp
            ### AHRQ infection categories were added to the max adjusted model (5/8/23)
            + ahrq_infcat_pulmonary + ahrq_infcat_cns + ahrq_infcat_intra_abdominal + ahrq_infcat_genitourinary
            + ahrq_infcat_bone_joint + ahrq_infcat_septicemia_bacteremia + ahrq_infcat_skin_soft_tissue_infection
            + ahrq_infcat_other + ahrq_infcat_obstetric_gynecologic
            ### Year, beta-lactam allergy, and spoken language of English (6/30/23)
            + edArrivalYR + allergy_bl + primLangEnglish)
    } else if (adj == "baseline_HCA") {
        (A ~ t0_cat
            # ----- Time-invariant covariates only -----
            ### HCA does not have firstCodeStatusGrp
            + hospitalType + agecat + raceEthnicity + sex + elix_AIDS + elix_ALCOHOL + elix_ANEMDEF 
            + elix_AUTOIMMUNE + elix_BLDLOSS + elix_CANCER_LEUK + elix_CANCER_LYMPH+ elix_CANCER_METS 
            + elix_CANCER_SOLID + elix_COAG + elix_DEMENTIA + elix_DEPRESS + elix_DIAB_CX +  elix_DIAB_UNCX 
            + elix_DRUG_ABUSE + elix_HF + elix_HTN_CX + elix_HTN_UNCX + elix_LIVER_MLD + elix_LIVER_SEV 
            + elix_LUNG_CHRONIC + elix_NEURO_MOVT + elix_NEURO_OTH + elix_NEURO_SEIZ + elix_OBESE 
            + elix_PARALYSIS + elix_PERIVASC + elix_PSYCHOSES + elix_PULMCIRC + elix_RENLFL_MOD + elix_RENLFL_SEV
            + elix_THYROID_HYPO + elix_THYROID_OTH + elix_ULCER_PEPTIC + elix_VALVE + elix_WGHTLOSS + elix_CBVD 
            + elixcat + dischargePrior90Days + insType + facility_admission + ems + intubatedPreHosp
            ### AHRQ infection categories are named slightly differently in HCA
            + ahrq_pulmonary + ahrq_cns + ahrq_intra_abdominal + ahrq_genitourinary
            + ahrq_bone_joint + ahrq_septicemia + ahrq_skin_tissue
            + ahrq_other + ahrq_obs_gyno
            ### HCA does not have spoken language
            + edArrivalYR + allergy_bl)
    } else if (adj == "none") {
        A ~ t0_cat
    } else {
        stop(paste0("Invalid value for adj: ", adj))
    }
    glm(form,
        data = dt[Alag1 == 0 & t0 < Kint, ],
        family = binomial(link = "logit"))
}

# Create point estimates for the weighted estimates of intervention hazards of death and discharge
# NOTE: this destructively modifies the data.table `dt`!
#   - `Kint` is the number of time intervals
#   - `cutTimes` is a vector with all valid time intervals, usually c(0:(Kint-1))
#   - `compareTimes` is a two-column matrix with start and end intervals for each time window to compare
#   - `adj` is the name of the adjustment model (used for `ipw_glm_denom_prob()`)
#   - `intSize` is the size of each time interval in hours (by default, 0.5 hours)
ipw_est <- function(dt, Kint, cutTimes, compareTimes, adj = "max", intSize = 0.5, out_debug = NULL) {
    denom_prob <- ipw_glm_denom_prob(dt, Kint, adj)
    if (!is.null(out_debug)) { cat_append(summary(denom_prob), out_debug) }
    dt.est <- data.table(t0 = cutTimes)

    dt[, p0denom := 9999]
    dt[, p1denom := 9999]
    dt[
        A == 0 & Alag1 == 0 & t0 < Kint, 
        p0denom := (1 - (predict(denom_prob, dt[A == 0 & Alag1 == 0 & t0 < Kint, ], type = "response")))
    ]
    dt[
        A == 1 & Alag1 == 0 & t0 < Kint, 
        p1denom := predict(denom_prob, dt[A == 1 & Alag1 == 0 & t0 < Kint, ], type = "response")
    ]
    
    # For each time window that we will compare ...
    for (i in seq_len(nrow(compareTimes))) {
        compare_lower <- compareTimes[i, 1]
        compare_upper <- compareTimes[i, 2]
        hours_lower <- compare_lower * intSize
        hours_upper <- (compare_upper + 1) * intSize
        hours <- paste0(hours_lower, "_", hours_upper)
        weights_col <- paste0("weights_", hours)
        if (!is.null(out_debug)) { 
            cat_append(paste0("\n#### WEIGHTS ", hours_lower, "-", hours_upper, "h ####\n"), out_debug)
        }
        
        # Calculate weights
        dt[, (weights_col) := int_weights(dt, compare_lower, compare_upper)]
        if (!is.null(out_debug)) { cat_append(summary(dt[[weights_col]]), out_debug) }
        # Clamp weight columns at the 99th percentile
        cut99 <- quantile(dt[[weights_col]], probs = c(.99))
        dt[dt[[weights_col]] > cut99, (weights_col) := cut99]
        if (!is.null(out_debug)) { cat_append(summary(dt[[weights_col]]), out_debug) }
        # Calculate cumulative hazards of death and discharge
        haz_death <- non_parametric_cum_haz(dt[[weights_col]], dt, cutTimes, TRUE)
        haz_discharge <- non_parametric_cum_haz(dt[[weights_col]], dt, cutTimes, FALSE)
        # Calculate weighted estimates of cumulative incidence of death and discharge
        dt.est[, (paste0("death_cum_inc_", hours)) := 
            non_parametric_cum_inc(haz_death, haz_discharge, t0, FALSE)]
        dt.est[, (paste0("discharge_cum_inc_", hours)) := 
            non_parametric_cum_inc(haz_death, haz_discharge, t0, TRUE)]
    }
    
    # Save debug information on patient counts
    if (!is.null(out_debug)) { debug_counts(dt, out_debug) }
    
    # Unweighted analysis
    dt[, weights_nc := 1]
    haz_death_nc <- non_parametric_cum_haz(dt[, weights_nc], dt, cutTimes, TRUE)
    haz_discharge_nc <- non_parametric_cum_haz(dt[, weights_nc], dt, cutTimes, FALSE)
    dt.est[, `:=`(
        death_cum_inc_nc = non_parametric_cum_inc(haz_death_nc, haz_discharge_nc, t0, FALSE),
        discharge_cum_inc_nc = non_parametric_cum_inc(haz_death_nc, haz_discharge_nc, t0, TRUE)
    )]
    
    # We only care about the absolute risk at the max timepoint
    dt.est <- dt.est[t0 == tail(cutTimes, n = 1), ]
    dt.est
}

# Resample the dataset (with replacements) and rerun a point estimate for the IPW model
ipw_boot_resamp <- function(dt, Kint, cutTimes, compareTimes, adj = "max", id_col = "id") {
    setkeyv(dt, id_col)
    ids <- unique(dt[[id_col]])
    ids_star <- sample(ids, replace = TRUE)

    dt.samp_ids <- data.table(samp_id = seq_along(length(ids)))
    dt.samp_ids[, (id_col) := ids_star]
    setkeyv(dt.samp_ids, id_col)

    dt.resamp <- dt[dt.samp_ids, allow.cartesian = TRUE]
    dt.resamp[, old_id := dt.resamp[[id_col]]]
    dt.resamp[, (id_col) := samp_id]
    dt.est <- ipw_est(dt.resamp, Kint, cutTimes, compareTimes, adj)
    return(dt.est)
}

# Run a bootstrap of the weighted estimates of intervention hazards of death and discharge
# This returns a data.table with the results for each bootstrap replicate
# The arguments are similar to `ipw_est()`
ipw_boot <- function(dt, Kint, cutTimes, compareTimes, n_reps = 1000, adj = "max", id_col = "id", 
        .parallel = FALSE) {
    reps <- 1:n_reps
    if (.parallel) { require(doFuture) }
    prog <- progressor(along = reps)
    dts <- llply(reps, function(x) { 
        dt.est <- ipw_boot_resamp(dt, Kint, cutTimes, compareTimes, adj, id_col) 
        prog(sprintf("finished resample # %g", x))
        dt.est[, boot := x]
        setcolorder(dt.est, c("boot", "t0"))
        dt.est
    }, .parallel = .parallel)
    rbindlist(dts)
}

##########
# The summarize_* functions convert the above outputs, which are cumulative hazards of death and discharge,
# into absolute risks, relative risks, and absolute risk differences vs the 1st interval (easier to interpret)
##########

### Summarizes an IPW point estimate for a single sample -- see ipw_est()
### Additional columns for relative risk and absolute risk difference are also calculated
### The possible baseline intervals are any of the intervals starting at time 0
summarize_ipw_est <- function(dt.est, compareTimes, intSize = 0.5) {
    dt.est.more <- copy(dt.est)
    orig_measures <- grep("_cum_inc_", colnames(dt.est), value=TRUE)
    compare_intervals <- cbind(compareTimes[, 1] * intSize, (compareTimes[, 2] + 1) * intSize)
    baseline_intervals <- compare_intervals[compare_intervals[, 1] == 0, , drop=FALSE]
    relative_intervals <- apply(compare_intervals[compare_intervals[, 1] != 0, , drop=FALSE], 1,
        function(row) { paste0(row, collapse="_") })
    relative_interval_regexp <- paste0("(", paste0(relative_intervals, collapse = "|"), ")$")
    relative_measures <- grep(relative_interval_regexp, orig_measures, value=TRUE)
    for (rel_meas in relative_measures) {
        hour_lo <- as.numeric(strsplit(rel_meas, split="_")[[1]][4])
        allowed_baselines <- baseline_intervals[baseline_intervals[, 2] <= hour_lo, , drop=FALSE]
        for (i in seq_len(allowed_baselines)) {
            rel_to_interval <- paste0(allowed_baselines[i, 1], "_", allowed_baselines[i, 2])
            rel_to <- sub(relative_interval_regexp, rel_to_interval, rel_meas)
            dt.est.more[, (paste0(rel_meas, ".rr.", rel_to_interval)) := .SD[[rel_meas]] / .SD[[rel_to]]]
            dt.est.more[, (paste0(rel_meas, ".arr.", rel_to_interval)) := .SD[[rel_meas]] - .SD[[rel_to]]]
        }
    }
    
    dt.melted <- melt(dt.est.more, id.vars = c("t0"), variable.name = "measure", value.name = "val")
    dt.summary <- dt.melted[measure %in% orig_measures, .(abs_risk = val), by = "measure"]
    for (i in rev(seq_len(nrow(baseline_intervals)))) {
        rel_to_interval <- paste0(baseline_intervals[i, 1], "_", baseline_intervals[i, 2])
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".rel_risk")) := 
                dt.melted[measure == paste0(.BY$measure, ".rr.", rel_to_interval), val]]
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".abs_risk_diff")) := 
                dt.melted[measure == paste0(.BY$measure, ".arr.", rel_to_interval), val]]
    }
    dt.summary
}

### Summarizes the IPW bootstrapped output -- see ipw_boot()
### Because this output is a distribution of estimates, medians and 95% CIs are provided
### Of note, we merge this with the summarize_ipw_est() output above for ease of interpretation
summarize_ipw_boot <- function(dt.boot, dt.est.summary, compareTimes, intSize = 0.5) {
    dt.boot.more <- copy(dt.boot)
    orig_measures <- grep("_cum_inc_", colnames(dt.boot), value=TRUE)
    compare_intervals <- cbind(compareTimes[, 1] * intSize, (compareTimes[, 2] + 1) * intSize)
    baseline_intervals <- compare_intervals[compare_intervals[, 1] == 0, , drop=FALSE]
    relative_intervals <- apply(compare_intervals[compare_intervals[, 1] != 0, , drop=FALSE], 1,
        function(row) { paste0(row, collapse="_") })
    relative_interval_regexp <- paste0("(", paste0(relative_intervals, collapse = "|"), ")$")
    relative_measures <- grep(relative_interval_regexp, orig_measures, value=TRUE)
    for (rel_meas in relative_measures) {
        hour_lo <- as.numeric(strsplit(rel_meas, split="_")[[1]][4])
        allowed_baselines <- baseline_intervals[baseline_intervals[, 2] <= hour_lo, , drop=FALSE]
        for (i in seq_len(nrow(allowed_baselines))) {
            rel_to_interval <- paste0(allowed_baselines[i, 1], "_", allowed_baselines[i, 2])
            rel_to <- sub(relative_interval_regexp, rel_to_interval, rel_meas)
            dt.boot.more[, (paste0(rel_meas, ".rr.", rel_to_interval)) := .SD[[rel_meas]] / .SD[[rel_to]]]
            dt.boot.more[, (paste0(rel_meas, ".arr.", rel_to_interval)) := .SD[[rel_meas]] - .SD[[rel_to]]]
        }
    }

    dt.melted <- melt(dt.boot.more, id.vars = c("boot", "t0"), variable.name = "measure", value.name = "val")
    dt.summary <- dt.melted[
        measure %in% orig_measures,
        .(
            abs_risk.med = median(val),
            abs_risk.ci95_lo = quantile(val, 0.025),
            abs_risk.ci95_up = quantile(val, 0.975)
        ),
        by = "measure"
    ]
    for (i in rev(seq_len(nrow(baseline_intervals)))) {
        rel_to_interval <- paste0(baseline_intervals[i, 1], "_", baseline_intervals[i, 2])
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".abs_risk_diff.med")) := 
                median(dt.melted[measure == paste0(.BY$measure, ".arr.", rel_to_interval), val])]
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".abs_risk_diff.ci95_lo")) :=
                quantile(dt.melted[measure == paste0(.BY$measure, ".arr.", rel_to_interval), val], 0.025)]
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".abs_risk_diff.ci95_up")) :=
                quantile(dt.melted[measure == paste0(.BY$measure, ".arr.", rel_to_interval), val], 0.975)]
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".rel_risk.med")) :=
                median(dt.melted[measure == paste0(.BY$measure, ".rr.", rel_to_interval), val])]
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".rel_risk.ci95_lo")) :=
                quantile(dt.melted[measure == paste0(.BY$measure, ".rr.", rel_to_interval), val], 0.025)]
        dt.summary[measure %in% relative_measures, by = "measure",
            (paste0("vs_", rel_to_interval, ".rel_risk.ci95_up")) := 
                quantile(dt.melted[measure == paste0(.BY$measure, ".rr.", rel_to_interval), val], 0.975)]
    }
        
    # Join in the point estimates summary and add some extra initial columns for readability
    dt.summary <- dt.est.summary[dt.summary, on = "measure"]
    dt.summary[, measure := as.character(measure)]
    dt.summary[, outcome := unlist(lapply(strsplit(measure, split="_"), function(x) x[1]))]
    dt.summary[
        !grepl("_nc$", measure), 
        `:=`(
            hour_lo = as.numeric(unlist(lapply(strsplit(measure, split="_"), function(x) x[4]))),
            hour_up = as.numeric(unlist(lapply(strsplit(measure, split="_"), function(x) x[5])))
        ),
    ]
    setcolorder(dt.summary, sort(colnames(dt.summary)))
    setcolorder(dt.summary, c("measure", "outcome", "hour_lo", "hour_up"))
    dt.summary[order(outcome, hour_lo, -hour_up)]
}