#' Studies of pain relief medications for osteoarthritis
#'
#' A dataset containing results on the WOMAC pain scale (0-10) over time for studies investigating 29
#' treatments for pain relief in patients with osteoarthritis. Standard deviations have been imputed for
#' 269 observations.
#'
#' `osteopain` is a data frame in long format (one row per observation, arm and study),
#' with the variables `studyID`, `time`, `y`, `se`, `treatment`, `arm` and `treatname`.
#'
#' @format A data frame with 417 rows and 7 variables:
#' * `studyID` Study identifiers
#' * `time` Numeric data indicating follow-up times
#' * `y` Numeric data indicating the mean response for a given observation
#' * `se` Numeric data indicating the standard error for a given observation
#' * `treatment` Treatment identifiers as factors. Labels are shortened treatment names.
#' * `arm` Arm identifiers coded for each study
#' * `treatname` Character data giving the full names of each treatment
#'
#' @source Pfizer Ltd.
"osteopain"





#' Studies of combined treatments for reducing serum uric acid in patients with gout
#'
#' A dataset from a systematic review of interventions for lowering Serum Uric Acid (SUA) concentration in
#' patients with gout **(not published previously)**. The outcome is continuous, and aggregate data responses
#' correspond to the mean change from baseline in SUA in mg/dL.
#' Treatments with similar doses have been pooled together to improve
#' network connectivity and facilitate evidence synthesis, resulting in 19 treatments of 7 agents included
#' in the network. Standard deviations have been imputed for 181 observations.
#'
#' `goutSUA_CFBcomb` is a data frame in long format (one row per observation, arm and study),
#' with the variables `studyID`, `time`, `y`, `se`, `treatment`, `treatname` and `class`.
#'
#' @format A data frame with 224 rows and 7 variables:
#' * `studyID` Study identifiers
#' * `time` Numeric data indicating follow-up times
#' * `y` Numeric data indicating the mean response for a given observation
#' * `se` Numeric data indicating the standard error for a given observation
#' * `treatment` Treatment identifiers as factors. Labels are shortened treatment names.
#' * `treatname` Character data giving the full names of each treatment in the format agent_dose
#' * `class` Shortened agent names stored as factors.
#'
#' @source Pfizer Ltd.
"goutSUA_CFBcomb"





#' Studies of treatments for reducing serum uric acid in patients with gout
#'
#' A dataset from a systematic review of interventions for lowering Serum Uric Acid (SUA) concentration in
#' patients with gout **(not published previously)**. The outcome is continuous, and aggregate data responses
#' correspond to the mean change from baseline in SUA in mg/dL.
#' Overall there are 41 treatments of 8 agents in the network. Standard
#' deviations have been imputed for 181 observations.
#'
#' `goutSUA_CFB` is a data frame in long format (one row per observation, arm and study),
#' with the variables `studyID`, `time`, `y`, `se`, `treatment`, `treatname` and `class`.
#'
#' @format A data frame with 224 rows and 7 variables:
#' * `studyID` Study identifiers
#' * `time` Numeric data indicating follow-up times
#' * `y` Numeric data indicating the mean response for a given observation
#' * `se` Numeric data indicating the standard error for a given observation
#' * `treatment` Treatment identifiers as factors. Labels are shortened treatment names.
#' * `treatname` Character data giving the full names of each treatment in the format agent_dose
#' * `class` Shortened agent names stored as factors.
#'
#' @source Pfizer Ltd.
"goutSUA_CFB"





#' Studies of treatments for reducing body weight in patients with obesity
#'
#' A dataset from a systematic review of pharmacological treatments for reducing body weight in patients with
#' obesity. The outcome is continuous, and aggregate data responses are given as mean change from baseline in
#' body weight (KG). Overall there are 35 RCTs investigating
#' 26 treatments of 16 agents (/combinations of agents) in the network. Standard
#' deviations have been imputed for 421 observations.
#'
#' `obesityBW_CFB` is a data frame in long format (one row per observation, arm and study),
#' with the variables `studyID`, `time`, `y`, `se`, `N`, `treatment`, `treatname`, `agent` and `class`.
#'
#' @format A data frame with 710 rows and 7 variables:
#' * `studyID` Study identifiers
#' * `time` Numeric data indicating follow-up times
#' * `y` Numeric data indicating the mean response for a given observation
#' * `se` Numeric data indicating the standard error for a given observation
#' * `N` Numeric data indicating the number of participants used to calculate means for each observation
#' * `treatment` Treatment identifiers as factors. Labels are shortened treatment names.
#' * `treatname` Character data giving the full names of each treatment in the format agent_dose
#' * `agent` Agent (drug) names stored as characters
#' * `class` The drug class of the agent (a broader category than `agent`) stored as characters
#'
#' @source Pfizer Ltd.
"obesityBW_CFB"





#' Studies of alogliptin for lowering blood glucose concentration in patients with type II diabetes
#'
#' A dataset from a systematic review of Randomised-Controlled Trials (RCTs) comparing different doses of
#' alogliptin with placebo \insertCite{langford2016}{MBNMAtime}. The systematic review was simply performed and was intended to
#' provide data to illustrate a statistical methodology rather than for clinical inference. Alogliptin is
#' a treatment aimed at reducing blood glucose concentration in type II diabetes. The outcome is continuous,
#' and aggregate data responses correspond to the mean change in HbA1c from baseline to follow-up.
#' The dataset includes 14 Randomised-Controlled Trials (RCTs), comparing 5
#' different doses of alogliptin with placebo, leading to 6 different treatments (combination of dose and agent)
#' within the network.
#'
#' `alog_pcfb` is a data frame in long format (one row per observation, arm and study),
#' with the variables `studyID`, `clinicaltrialGov_ID`, `agent`, `dose`, `treatment`, `time`, `y`, `se`, and `N`.
#'
#' @format A data frame in long format (one row per arm and study), with 46 rows and 9 variables:
#' * `studyID` Study identifiers
#' * `clinicaltrialGov_ID` The clinicaltrial.gov ID code
#' * `agent` Character data indicating the agent to which participants were randomised
#' * `dose` Numeric data indicating the standardised dose received
#' * `treatment` Character data indicating the treatment (combination of agent and dose) to which participants were randomised
#' * `time` Numeric data indicating the time at which the observation was measured (given in weeks)
#' * `y` Numeric data indicating the mean change from baseline in blood glucose concentration (mg/dL) in a study arm
#' * `se` Numeric data indicating the standard error for the mean change from baseline in blood glucose concentration (mg/dL) in a study arm
#' * `N` Numeric data indicating the number in each arm at each follow-up time
#'
#' @references
#' \insertAllCited{}
#'
"alog_pcfb"






#' Studies comparing Tiotropium, Aclidinium and Placebo for maintenance treatment of moderate to severe chronic obstructive pulmonary disease
#'
#' A dataset from a systematic review of Randomised-Controlled Trials (RCTs) for maintenance treatment of moderate to severe chronic
#' obstructive pulmonary disease (COPD) \insertCite{karabis2013}{MBNMAtime}. Data are extracted from \insertCite{tallarita2019}{MBNMAtime}.
#' SEs were imputed for three studies, and number of patients randomised were imputed for one study (LAS 39) in which they were missing,
#' using the median standard deviation calculated from other studies in the
#' dataset. The outcome is trough Forced Expiratory Volume in 1 second (FEV1), measured in litres and reported in each study arm as mean
#' change from baseline to follow-up. The dataset includes 13 Randomised-Controlled Trials (RCTs), comparing 2 treatments (Tiotropium and
#' Aclidinium) and placebo.
#'
#' `copd` is a data frame in long format (one row per observation, arm and study),
#' with the variables `studyID`, `time`, `y`, `se`, `treatment`, and `n`.
#'
#' @format A data frame in long format (one row per arm and study), with 80 rows and 6 variables:
#' * `studyID` Study identifiers
#' * `time` Numeric data indicating the time at which the observation was measured (given in weeks)
#' * `y` Numeric data indicating the mean change from baseline in FEV1 (litres) in a study arm
#' * `se` Numeric data indicating the standard error for the mean change from baseline in FEV1 in a study arm
#' * `treatment` Factor data indicating the treatment to which participants were randomised
#' * `n` Numeric data indicating the number of participants randomised to each arm
#'
#' @references
#' \insertAllCited{}
#'
"copd"
