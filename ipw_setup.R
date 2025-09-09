####################################################################################################
# Description: This contains setup functions for the input person-time datasets for ipw_run.R
# The most important exported function is `ipw_setup()`, which runs all of the other functions
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
})

# Ensures that column names are named consistently between datasets from the different sources
ipw_normalize_hca_colnames <- function(dt) {
    arvyr_cols <- grep("^arvYR_[0-9]+$", colnames(dt), value=TRUE)
    if ("ExclOutsideHosp" %in% colnames(dt)) {
        setnames(dt, old = "ExclOutsideHosp", new = "exclOutsideHosp")
    }
    if ("excl_def_CC6hrs" %in% colnames(dt)) {
        setnames(dt, old = "excl_def_CC6hrs", new = "excl_def_CC6Hrs")
    }
    if ("ivabxInit" %in% colnames(dt)) {
        setnames(dt, old = "ivabxInit", new = "ivAbxInit")
    }
    if ("defToIVabx6hrs" %in% colnames(dt)) {
        setnames(dt, old = "defToIVabx6hrs", new = "defToIVabx6Hrs")
    }
    if (length(arvyr_cols) > 0) {
        if (!("edArrivalYR" %in% colnames(dt))) {
            message("creating edArrivalYR")
            dt[, edArrivalYR := as.integer(NA)]
            for (col in arvyr_cols) {
                as_year <- as.integer(sub("^arvYR_", "", col))
                dt[dt[[col]] == 1, edArrivalYR := as_year]
            }
        }
        dt[, (arvyr_cols) := NULL]
    }
    return(dt)
}

# Renames/copies certain columns to the short names expected by IPW functions
ipw_add_standard_colnames <- function(dt) {
    dt[, cumsumA := ivAbxInitCumSum]
    dt[, A := ivAbxInit]
    dt[, Alag1 := lag1ivAbxInit]
    dt[t == 0, Alag1 := 0]
    dt[, t0 := t]
    dt[, t0_cat := as.factor(t0)]
    dt[, D := discharge]
    dt[, Y := death]
    dt[, id := hospID]
    dt[, edArrivalYR := as.factor(edArrivalYR)]
    setkey(dt, "hospID")
}

# Ensure appropriate exclusions are applied to the the input dataset
ipw_apply_exclusions <- function(dt) {
    # apply main exclusions
    dt.new <- dt[exclOutsideHosp==0 & excl_def_CC6Hrs==0 & exclAdmitPsych==0 & 
            exclAdmitOB==0 & excl_def_vitals==0 & exclKeyLabs==0 & excl_def_POIVabx==0 & 
            exclfirstIVabxDischarge==0 & exclCovDay1or2==0, ]

    # create subset based on infectionCriteria (Suspected Infection/Sepsis/Septic Shock)
    # furthermore, restrict to those receiving IV antibiotics within 6h of the ST0 timepoint
    dt.new <- dt.new[infectionCriteria==1 & defToIVabx6Hrs==1, ]
    
    # Added 2025-03-28: ensure patients who died or were discharged before 6h are excluded
    # using the `death`, `discharge`, and `t` columns, not just the above flags, as some
    # exceptions have been found
    bad_ids <- unique(dt.new[(discharge == 1 | death == 1) & t < 11, hospID])
    dt.new <- dt.new[!(hospID %in% bad_ids), ]
    
    return(dt.new)
}

# Create categorical versions of age, Elixhauser, and allergy variables
ipw_categorize_age_elix_and_allergy_bl <- function(dt) {
    # Categorical age
    if (!("agecat" %in% colnames(dt))) {
        dt[, agecat := as.factor(case_when(
            encAge < quantile(encAge, c(.25, .5, .75))[1] ~ 1,
            encAge < quantile(encAge, c(.25, .5, .75))[2] ~ 2,
            encAge < quantile(encAge, c(.25, .5, .75))[3] ~ 3,
            TRUE ~ 4))]
    }

    if (!("agesubjcat" %in% colnames(dt))) {
        dt[, agesubjcat := as.factor(case_when(
            encAge %in% 0:40 ~ 0,
            encAge %in% 41:50 ~ 1,
            encAge %in% 51:60 ~ 2,
            encAge %in% 61:70 ~ 3,
            TRUE ~ 4))]
    }

    # Categorical Elixhauser score
    if (!("elixcat" %in% colnames(dt))) {
        dt[, elixcat := as.factor(case_when(
            elix_index_mortality < quantile(elix_index_mortality, c(.25, .5, .75))[1] ~ 1,
            elix_index_mortality < quantile(elix_index_mortality, c(.25, .5, .75))[2] ~ 2,
            elix_index_mortality < quantile(elix_index_mortality, c(.25, .5, .75))[3] ~ 3,
            TRUE ~ 4))]
    }
    
    # Categorical allergy to beta-lactams
    if (("allergy_bl" %in% colnames(dt)) && !is.factor(dt[["allergy_bl"]])) {
        dt[, allergy_bl := addNA(factor(allergy_bl, levels=1:4,
            labels=c("1: PCN allergy", "2: Non-PCN allergy", "3: Other abx allergy", 
                "4: More than one type of allergy")))]
    }

    return(dt)
}

# Create categorical versions of the vitals variables
ipw_categorize_vitals <- function(dt) {
    #Max temp
    if (!("tempMax_cat" %in% colnames(dt)) || !is.factor(dt[["tempMax_cat"]])) {
        dt[, tempMax_cat := cut(tempMax, c(79.99, 96.79, 100.39, Inf), 
            c("1: <96.8", "2: 96.8-100.3", "3: 100.4+"), include.lowest=TRUE)]
    }

    #Heart rate
    if (!("HR_cat" %in% colnames(dt)) || !is.factor(dt[["HR_cat"]])) {
        dt[, HR_cat := cut(HR, c(0, 39, 90, 120, Inf), 
            c("1: <40", "2: 40-90", "3: 91-120", "4 :121+"), include.lowest=TRUE)]
    }

    #Systolic blood pressure
    if (!("SBP_cat" %in% colnames(dt)) || !is.factor(dt[["SBP_cat"]])) {
        dt[, SBP_cat := cut(SBP, c(0, 39, 90, 120, 160, Inf), 
            c("1: <40", "2: 40-90", "3: 91-120", "4: 121-160", "5: 161+"), include.lowest=TRUE)]
    }

    #Respiratory rate
    if (!("RR_cat" %in% colnames(dt)) || !is.factor(dt[["RR_cat"]])) {
        dt[, RR_cat := cut(RR, c(0, 10, 20, 30, Inf), 
            c("1: <=10", "2: 11-20", "3: 21-30", "4: 31+"), include.lowest=TRUE)]
    }

    #GCS - Note that in HCA data, there is no uncut `GCS` column
    GCS_labels <- c("5: <6", "4: 6-9", "3: 10-12", "2: 13-14", "1: 15")
    if (!("GCS_cat" %in% colnames(dt))) {
        dt[, GCS_cat := addNA(cut(GCS, c(3, 5, 9, 12, 14, Inf), 
            GCS_labels, include.lowest=TRUE))]
    } else if (!is.factor(dt[["GCS_cat"]])) {
        dt[, GCS_cat := addNA(factor(GCS_cat, levels = 5:1, labels = GCS_labels))]
    }
}

# Create a categorical version of the O2 device requirement variable combined with spO2
ipw_categorize_o2 <- function(dt) {
    if (!("o2" %in% colnames(dt))) {
        dt[, o2 := case_when(
            o2Device == "NONE" & spO2 >= 95 ~ 0,
            o2Device == "NONE" & spO2 < 95 ~ 1,
            o2Device == "NASAL CANNULA" ~ 2,
            o2Device == "SIMPLE MASK" ~ 3,
            o2Device == "OXYMIZER" ~ 4,
            o2Device == "ADVANCED MASK" ~ 5,
            o2Device == "HIGH FLOW" ~ 6,
            o2Device == "BIPAP" ~ 7,
            o2Device == "VENTILATOR" ~ 8,
            o2Device == "ECMO" ~ 9)]
    }

    if (!("o2f" %in% colnames(dt)) || !is.factor(dt[["o2f"]])) {
        dt[, o2f := factor(if ("o2f" %in% colnames(dt)) o2f else o2, 
            levels = 0:9, 
            labels = c("No o2 device with sp02 >=95", "No o2 device with sp02 <95",
                "Nasal cannula", "Simple mask", "O2 conserving device", 
                "Non-rebreather", "High flow", "BiPAP", "Ventilator", "ECMO"))]
    }
}

# Create categorical versions of the labs variables
ipw_categorize_labs <- function(dt) {
    #Creatinine
    if (!("creat_cat" %in% colnames(dt)) || !is.factor(dt[["creat_cat"]])) {
        dt[, creat_cat := addNA(cut(creat, c(0, 1.199, 1.999, 3.499, 4.999, Inf), 
            c("1: <1.2", "2: 1.2-<2.0", "3: 2.0-<3.5", "4: 3.5-<5.0", "5: 5.0+"), include.lowest=TRUE))]
    }

    #Platelet
    if (!("plt_cat" %in% colnames(dt)) || !is.factor(dt[["plt_cat"]])) {
        dt[, plt_cat := addNA(cut(plt, c(0, 19, 49, 99, 149, Inf), 
            c("5: <20", "4: <50", "3: 50-99", "2: 100-149", "1: 150+"), include.lowest=TRUE))]
    }

    #Lactate
    if (!("lac_cat" %in% colnames(dt)) || !is.factor(dt[["lac_cat"]])) {
        dt[, lac_cat := addNA(cut(lac, c(0, 1.999, 2.999, 3.999, Inf), 
            c("1: <2", "2: 2.0-<3.0", "3: 3.0-<4.0", "4: 4.0+"), include.lowest=TRUE))]
    }

    #Hematocrit
    if (!("hct_cat" %in% colnames(dt)) || !is.factor(dt[["hct_cat"]])) {
        dt[, hct_cat := addNA(cut(hct, c(0, 20.9, 27.9, 35.9, Inf), 
            c("1: <21.0", "2: 21.0-27.9", "3: 28.0-35.9", "4: 36.0+"), include.lowest=TRUE))]
    }

    #WBC
    if (!("wbc_cat" %in% colnames(dt)) || !is.factor(dt[["wbc_cat"]])) {
        dt[, wbc_cat := addNA(cut(wbc, c(0, 0.999, 3.999, 11.999, 19.999, Inf), 
            c("1: <1", "2: 1.0<-4.0", "3: 4.0<-12.0", "4: 12.0-<20", "5: 20+"), include.lowest=TRUE))]
    }

    #Bilirubin
    if (!("tbili_cat" %in% colnames(dt)) || !is.factor(dt[["tbili_cat"]])) {
        dt[, tbili_cat := addNA(cut(tbili, c(0, 1.199, 1.999, 5.999, 11.999, Inf), 
            c("1: <1.2", "2: 1.2-<2.0", "3: 2.0-<6.0", "4: 6.0-<12.0", "5: 12.0+"), include.lowest=TRUE))]
    }

    #Anion gap
    if (!("ag_cat" %in% colnames(dt)) || !is.factor(dt[["ag_cat"]])) {
        dt[, ag_cat := addNA(cut(ag, c(0, 11.9, 16.0, 23.0, Inf), 
            c("1: <12", "2: 12-16", "3: 17-23", "4: 24+"), include.lowest=TRUE))]
    }

    #Ast
    if (!("ast_cat" %in% colnames(dt)) || !is.factor(dt[["ast_cat"]])) {
        dt[, ast_cat := addNA(cut(ast, c(0, 49, 99, 199, 499, Inf), 
            c("1: <50", "2: 50-99", "3: 100-199", "4: 200-499", "5: 500+"), include.lowest=TRUE))]
    }

    #Sodium
    if (!("sdm_cat" %in% colnames(dt)) || !is.factor(dt[["sdm_cat"]])) {
        dt[, sdm_cat := addNA(cut(sdm, c(0, 129, 150, Inf), 
            c("1: <130", "2: 130-150", "3: 151+"), include.lowest=TRUE))]
    }

    #Albumin - Note that in HCA data, there is no uncut `alb` column
    alb_levels <- c("1: 0-1.9", "2: 2.0-2.9", "3: 3-3.5", "4: 3.6+")
    if (!("alb_cat" %in% colnames(dt))) {
        dt[, alb_cat := addNA(cut(alb, c(0, 1.9, 2.9, 3.5, Inf), 
            alb_levels, include.lowest=TRUE))]
    } else if (!is.factor(dt[["alb_cat"]])) {
        dt[, alb_cat := addNA(factor(alb_cat, levels = 1:4, labels = alb_levels))]
    }

    #Glucose - Note that in HCA data, there is no uncut `glu` column
    glu_levels <- c("1: <60", "2: 60-119", "3: 120-199", "4: 200-299", "5:300+")
    if (!("glu_cat" %in% colnames(dt))) {
        dt[, glu_cat := addNA(cut(glu, c(0, 59.0, 119.0, 199.0, 299.0, Inf), 
            glu_levels, include.lowest=TRUE))]
    } else if (!is.factor(dt[["glu_cat"]])) {
        dt[, glu_cat := addNA(factor(glu_cat, levels = 1:5, labels = glu_levels))]
    }
}

# Create a categorical version of fluid bolus cumulative sum normalized to patient weight
ipw_normalize_fluids_to_weight <- function(dt) {
    if ("fluid_bolus_cumsum_cat" %in% colnames(dt)) { return(dt) }
    if (!("fluid_bolus_cumsum_dose" %in% colnames(dt))) { return(dt) }
    if (!("weightKg_imp" %in% colnames(dt))) { return(dt) }

    dt[, fluid_bolus_cumsum := replace_na(fluid_bolus_cumsum_dose, 0)]
    dt[order(t), fluid_bolus_cumsum := cumsum(fluid_bolus_cumsum), by = hospID]

    dt[, 
        fluid_bolus_cumsum_cat := addNA(
            cut(fluid_bolus_cumsum / weightKg_imp, c(0, 10, 20, 30, 40, Inf), 
            c("1: 0-9", "2: 10-19", "3: 20-29", "4: 30-39", "5:40+"), include.lowest=TRUE)
        )
    ]
}

# This is the main setup function that runs all of the above setup functions in the correct sequence
ipw_setup <- function(dt) {
    dt <- ipw_normalize_hca_colnames(dt)
    dt <- ipw_apply_exclusions(dt)
    ipw_add_standard_colnames(dt)
    
    ipw_categorize_age_elix_and_allergy_bl(dt)
    ipw_categorize_vitals(dt)
    ipw_categorize_o2(dt)
    ipw_categorize_labs(dt)
    ipw_normalize_fluids_to_weight(dt)

    return(dt)
}