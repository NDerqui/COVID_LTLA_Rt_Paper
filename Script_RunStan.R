#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Impact of vaccination on SARS-CoV-2 transmission in the UK: a modelling study #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


## Script to run model via rstan
# 1: Get data, clean and prepare stan data list
# 2: Run model
# 3: Retrieve pars of interest and save



# SET UP ------------------------------------------------------------------


rm(list = ls())
memory.limit(10000)

# Packages

library(dplyr)
library(rstudioapi)
library(StanHeaders)
library(rstan)
options(mc.cores = parallel::detectCores())

# Data

data_vax <- read.csv("Data/rtm_incoming_vaccination_20211116-173218-f36a1245_vacc_coverage_ltla.csv")
data_rt <- read.csv("Data/UK_hotspot_Rt_estimates.csv")
data_var <- read.csv("Data/vam_by_ltla.csv")



# OVERVIEW DATA ------------------------------------------------------------


#### Vaccination Data ####

names(data_vax)
dim(data_vax)

# Data on 306 LTLA vaccination, per 20 age groups, per 342 dates

l_code <- length(unique(data_vax$ltla_code))

l_agemin <- length(unique(data_vax$age_band_min))

l_date <- length(unique(data_vax$vaccination_date))

dim(data_vax)
l_code*l_date*l_agemin

# As described before, when age is NA, this is the sum of all age_groups

sum(data_vax$First[is.na(data_vax$age_band_min)],
    data_vax$Second[is.na(data_vax$age_band_min)],
    data_vax$Third[is.na(data_vax$age_band_min)])
sum(data_vax$First[!is.na(data_vax$age_band_min)],
    data_vax$Second[!is.na(data_vax$age_band_min)],
    data_vax$Third[!is.na(data_vax$age_band_min)])

# Set as date

data_vax$vaccination_date <- as.Date(data_vax$vaccination_date)


#### Rt Data ####

names(data_rt)
dim(data_rt)

# Data on 391 LTLS Rt, per 338 dates

l_name <- length(unique(data_rt$area))

l_date_2 <- length(unique(data_rt$date))

dim(data_rt)
l_name*l_date_2

# Set as date

data_rt$date <- as.Date(data_rt$date)


#### Var Data ####

names(data_var)
dim(data_var)

# Data on 314 LTLA and 871 dates

l_code2 <- length(unique(data_var$ltlacode))

l_date3 <- length(unique(data_var$date))

dim(data_var)
l_code2*l_date3


#### Clean ####

rm(l_agemin, l_code, l_code2, l_name, l_date, l_date_2, l_date3)



# SET PARAMETERS ----------------------------------------------------------


# Set date parameters

Date_Start <- as.Date("01/01/2021", format = "%d/%m/%Y")
Date_End <- as.Date("31/12/2021", format = "%d/%m/%Y")

# Covariates

covar_var <- c("Var_PreAlpha", "Var_Alpha", "Var_Delta")
covar_vax <- c("First_Prop", "Second_Prop", "Third_Prop")



# DATA FOR Stan PREP FUNC --------------------------------------------------


get_data <- function(data_vax, data_rt, data_var,
                     covar_var, covar_vax,
                     Date_Start, Date_End,
                     DoVariants, DoVaxVariants, DoAge, DoVaxAge) {
  
  # Substract data
  
  dvax <- data_vax
  drt <- data_rt
  dvar <- data_var
  dage <- data_vax
  
  
  ## Making Vax/Var/Rt Data comparable ##
  
  # No/Yes age groups
  
  dvax <- filter(dvax, is.na(age_band_min))
  
  dage <- filter(dage, !is.na(age_band_min))
  
  # Dates of interest
  
  dvax <- rename(dvax, date = vaccination_date)
  dvax <- filter(dvax, date >= Date_Start)
  dvax <- filter(dvax, date <= Date_End)
  
  drt <- filter(drt, date >= Date_Start)
  drt <- filter(drt, date <= Date_End)
  
  dvar <- filter(dvar, date >= Date_Start)
  dvar <- filter(dvar, date <= Date_End)
  
  dage <- rename(dage, date = vaccination_date)
  dage <- filter(dage, date >= Date_Start)
  dage <- filter(dage, date <= Date_End)
  
  # LTLA names
  
  drt <- rename(drt, ltla_name = area)
  dvar <- rename(dvar, ltla_code = ltlacode)
  
  # Three age groups
  
  dage <- dage %>%
    mutate(total = First/First_Prop) %>%
    mutate(group = case_when(age_band_min >= 15 & age_band_min <= 45 ~ "15-49",
                             age_band_min >= 50 & age_band_min <= 65 ~ "50-69",
                             age_band_min >= 70 & age_band_min <= 90 ~ "70plus"))
  
  #To calculate new group proportions, only one weekly entry per age band
  dage$week <- round(as.numeric(floor((dage$date - Date_Start)/7)), digits = 0)
  dage <- dage %>%
    mutate(combi = paste0(week, ltla_name, age_band_min))
  dage <- filter(dage, !duplicated(dage$combi))
  
  
  ## Merge ##
  
  # Select the vars of interest
  
  dvax <- select(dvax, "ltla_code", "ltla_name", "date",
                 "First_Prop", "Second_Prop", "Third_Prop")
  drt <- select(drt, "ltla_name", "date", "Rt")
  dvar <- select(dvar, "ltla_code", "date", "n_all_wildtype_variant",
                 "n_all_alpha_variant", "n_all_delta_variant")
  dage <- select(dage, "ltla_code", "ltla_name", "date", "group",
                 "First", "Second", "Third", "total")
  
  # Merge
  
  data_merge <- merge(dvax, drt, by = c("ltla_name", "date"), all = TRUE)
  data_merge <- merge(data_merge, dvar, by = c("ltla_code", "date"), all = TRUE)
  
  data_merge_age <- merge(dage, drt, by = c("ltla_name", "date"), all = TRUE)
  data_merge_age <- merge(data_merge_age, dvar, by = c("ltla_code", "date"), all = TRUE)
  
  
  ## Final Cleaning ##
  
  # No NA
  
  data_merge <- data_merge[complete.cases(data_merge),]
  data_merge_age <- data_merge_age[complete.cases(data_merge_age),]
  
  # Weekly dates: var for no. of weeks & only one obs per week
  
  #Var for the number of weeks
  data_merge$week <- round(as.numeric(floor((data_merge$date - min(data_merge$date))/7)), digits = 0)
  data_merge_age$week <- round(as.numeric(floor((data_merge_age$date - min(data_merge_age$date))/7)), digits = 0)
  
  #Combi for unique combination of LTLA*week
  data_merge <- data_merge %>%
    mutate(combi = paste0(week, "/", ltla_name))
  data_merge_age <- data_merge_age %>%
    mutate(combi = paste0(week, "/", group, "/", ltla_name))
  
  # Age group proportion: prop of vaccinated and age prop in each region
  
  data_merge_age <- data_merge_age %>%
    group_by(combi) %>%
    mutate(FirstDose = sum(First),
           SecondDose = sum(Second),
           ThirdDose = sum(Third)) %>%
    mutate(First_Prop = FirstDose / sum(total),
           Second_Prop = SecondDose / sum(total),
           Third_Prop = ThirdDose / sum(total)) %>%
    ungroup() %>%
    arrange(ltla_name, week, group)
  
  # Calculate overall proportion in age group, across all LTLAs - chose 1 timepoint
  
  age_prop <- data_merge_age %>%
    filter(week == 1) %>%
    group_by(group) %>%
    mutate(total_group = sum(total)) %>%
    ungroup() %>%
    mutate(total_total = sum(total)) %>%
    group_by(group) %>%
    mutate(age_prop = total_group/total_total) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    arrange(group) %>%
    select(age_prop)
  
  #Remove the duplicates
  
  data_merge <- filter(data_merge, !duplicated(data_merge$combi))
  data_merge_age <- filter(data_merge_age, !duplicated(data_merge_age$combi))
  
  # Create vars for the proportion of alpha vs delta
  
  data_merge <- data_merge %>%
    mutate(total = (n_all_wildtype_variant + n_all_alpha_variant + n_all_delta_variant)) %>%
    mutate(Var_PreAlpha = case_when(total == 0 ~ 0,
                                    total != 0 ~ n_all_wildtype_variant/total)) %>%
    mutate(Var_Alpha = case_when(total == 0 ~ 1,
                                 total != 0 ~ n_all_alpha_variant/total)) %>%
    mutate(Var_Delta = case_when(total == 0 ~ 0,
                                 total != 0 ~ n_all_delta_variant/total))
  
  data_merge_age <- data_merge_age %>%
    mutate(total = (n_all_wildtype_variant + n_all_alpha_variant + n_all_delta_variant)) %>%
    mutate(Var_PreAlpha = case_when(total == 0 ~ 0,
                                    total != 0 ~ n_all_wildtype_variant/total)) %>%
    mutate(Var_Alpha = case_when(total == 0 ~ 1,
                                 total != 0 ~ n_all_alpha_variant/total)) %>%
    mutate(Var_Delta = case_when(total == 0 ~ 0,
                                 total != 0 ~ n_all_delta_variant/total))
  
  # Select only the ltla with max no of weeks
  
  data_merge <- data_merge %>%
    group_by(ltla_name) %>%
    mutate(max_week = length(unique(week))) %>%
    ungroup() %>%
    filter(max_week == length(unique(week)))
  
  data_merge_age <- data_merge_age %>%
    group_by(ltla_name) %>%
    mutate(max_week = length(unique(week))) %>%
    ungroup() %>%
    filter(max_week == length(unique(week)))
  
  # Select cols
  
  data_model <- select(data_merge, "ltla_name", "date", "week", "Rt",
                       "First_Prop", "Second_Prop", "Third_Prop",
                       "Var_PreAlpha", "Var_Alpha", "Var_Delta")
  data_model <- arrange(data_model, ltla_name, week)
  data_model_age <- select(data_merge_age, "ltla_name", "date", "week", "Rt","group",
                           "First_Prop", "Second_Prop", "Third_Prop",
                           "Var_PreAlpha", "Var_Alpha", "Var_Delta")
  data_model_age <- arrange(data_model_age, group, ltla_name, week)
  
  ## Now at this point, the data_model is dim timepoints * LTLAs;
  ## but the data_model_age is dim timepoints * LTLAs * age_geroups!
  
  ## This is useful to extract the correct information in each case,
  ## but it is NEVER the dimensions required to be passed to the model!!!
  
  
  ## Stan list data ##
  
  if (DoAge == 0) {
    
    ## Options without age
    
    dim(data_model)
    
    # No of LTLA
    
    NumLTLAs <- length(unique(data_model$ltla_name))
    NumLTLAs
    
    # No of LTLAs in the data
    
    NamesLTLAs <- unique(data_model$ltla_name)
    
    data_model$LTLAs <- NA
    for(i in 1:NumLTLAs){
      data_model$LTLAs[data_model$ltla_name == NamesLTLAs[i]] = i
    }
    
    LTLAs <- data_model$LTLAs
    LTLAs
    
    # No of weeks
    
    NumTimepoints <- length(unique(data_model$week))
    NumTimepoints
    
    Timepoints <- data_model$week
    Timepoints
    
    # No of weeks per LTLA
    
    NumWeeksByLTLA <- rep(NA, NumLTLAs)
    for(i in 1: NumLTLAs){   
      
      d_sub = data_model[data_model$ltla_name == NamesLTLAs[i], ]
      
      NumWeeksByLTLA[i] = length(unique(d_sub$week))
    }
    NumWeeksByLTLA
    
    # No of age groups
    
    NumGroup <- 1
    NumGroup
    NumVaxGroup <- 1
    NumVaxGroup
    
    # Names of age groups
    
    Groups <- rep(1, times = nrow(data_model))
    Groups
    
    # No of total obs 
    
    NumDatapoints <- NumLTLAs * NumTimepoints
    NumDatapoints
    
    nrow(data_model)
    
    # Variants
    
    if (DoVariants == 1) {
      NumVar <- length(covar_var)
      VarProp <- data_model[,covar_var]
    } else {
      NumVar <- 1
      VarProp <- as.data.frame(matrix(1, nrow = NumDatapoints))
    }
    
    # Data

    RtVals <- data_model$Rt
    VaxProp <- array(data = as.matrix(data_model[,covar_vax]),
                     dim = c(NumDatapoints, length(covar_vax), NumGroup)) # Should have a 3rd dimension with age groups
    VarProp <- array(data = as.matrix(VarProp),
                     dim = c(NumDatapoints, length(covar_var)))
    AgeProp <- array(1, dim = 1)
  
  } else {
    
    ## Options with age
    
    dim(data_model_age)
    
    # Rt, VarProp, etc should only be TimeRegion, thus...
    
    data_model_age_unique <- data_model_age %>%
      group_by(ltla_name, week) %>%
      filter(row_number() == 1) %>%
      ungroup()
    
    # No of LTLA
    
    NumLTLAs <- length(unique(data_model_age_unique$ltla_name))
    NumLTLAs
    
    # LTLAs in the data
    
    NamesLTLAs <- unique(data_model_age_unique$ltla_name)
    
    data_model_age_unique$LTLAs <- NA
    for(i in 1:NumLTLAs){
      data_model_age_unique$LTLAs[data_model_age_unique$ltla_name == NamesLTLAs[i]] = i
    }
    
    LTLAs <- data_model_age_unique$LTLAs
    LTLAs
    
    # No of weeks
    
    NumTimepoints <- length(unique(data_model_age_unique$week))
    NumTimepoints
    
    Timepoints <- data_model_age_unique$week
    Timepoints
    
    # No of weeks per LTLA
    
    NumWeeksByLTLA <- rep(NA, NumLTLAs)
    for(i in 1: NumLTLAs){   
      
      d_sub = data_model_age[data_model_age$ltla_name == NamesLTLAs[i], ]
      
      NumWeeksByLTLA[i] = length(unique(d_sub$week))
    }
    NumWeeksByLTLA
    
    # No of age groups
    
    NumGroup <- length(unique(data_model_age$group))
    NumGroup
    
    if (DoVaxAge == 1) {
      NumVaxGroup <- NumGroup} else {
      NumVaxGroup <- 1
    }
    NumVaxGroup
    
    # Names of age groups
    
    NamesGroups <- unique(data_model_age_unique$group)
    
    data_model_age_unique$Groups <- NA
    for(i in 1:NumGroup){
      data_model_age_unique$Groups[data_model_age_unique$group == NamesGroups[i]] = i
    }
    
    Groups <- data_model_age_unique$Groups
    Groups
    
    # No of total obs 
    
    NumDatapoints <- nrow(data_model_age_unique)
    NumDatapoints
    
    # Variants
    
    if (DoVariants == 1) {
      NumVar <- length(covar_var)
      VarProp <- data_model_age_unique[,covar_var]
    } else {
      NumVar <- 1
      VarProp <- as.data.frame(matrix(1, nrow = NumDatapoints))
    }
    
    # And the age should be added as a third dimension
    
    VaxProp <- data_model_age %>%
      select(ltla_name, week, group, First_Prop, Second_Prop, Third_Prop) %>%
      mutate(combi = paste0(ltla_name, "_", week)) %>%
      select(ltla_name, week, combi, group, First_Prop, Second_Prop, Third_Prop)
    
    x <- as.matrix(VaxProp[VaxProp$group == "15-49",covar_vax])
    y <- as.matrix(VaxProp[VaxProp$group == "50-69",covar_vax])
    z <- as.matrix(VaxProp[VaxProp$group == "70plus",covar_vax])
    
    VaxProp <- array(data = c(x, y, z), dim = c(NumDatapoints, length(covar_vax), NumGroup))
    rm(x,y,z)
   
    # Data
    
    RtVals <- data_model_age_unique$Rt
    VaxProp <- VaxProp
    VarProp <- array(data = as.matrix(VarProp),
                     dim = c(NumDatapoints, length(covar_var)))
    AgeProp <- age_prop$age_prop
  
  }
  
  # Num of time pars
  
  NumTrendPar <- NumTimepoints
  
  # Variants (common ground)
  
  if (DoVaxVariants == 1) {
    NumVaxVar <- NumVar
  } else {
    NumVaxVar <- 1
  }
  
  # Check point
  
  VaxProp[VaxProp > 1] <- 1
  
  # Stan list
  
  data_stan <- list(
    DoVariants = DoVariants,
    DoVaxVariants = DoVaxVariants,
    DoAge = DoAge,
    DoVaxAge = DoVaxAge,
    
    NumDatapoints = NumDatapoints,
    NumLTLAs = NumLTLAs,
    NumDoses = length(covar_vax),
    NumVar = NumVar,
    NumVaxVar = NumVaxVar,
    NumGroup = NumGroup,
    NumVaxGroup = NumVaxGroup,
    NumTimepoints = NumTimepoints,
    NumTrendPar = NumTrendPar,
    
    Timepoints = Timepoints,
    LTLAs = LTLAs,
    Groups = Groups,
    
    RtVals = RtVals,
    VaxProp = VaxProp,
    VarProp = VarProp,
    AgeProp = AgeProp,
    
    NamesLTLAs = NamesLTLAs
  )  
  
  return(data_stan)
}

# Run function to get data

data_stan <- get_data(data_vax = data_vax, data_rt = data_rt, data_var = data_var,
                      covar_var = covar_var, covar_vax = covar_vax,
                      Date_Start = Date_Start, Date_End = Date_End,
                      DoVariants = 1, DoVaxVariants = 0, DoAge = 0, DoVaxAge = 0)



# STAN MODEL --------------------------------------------------------------


#### Model ####

library(here)

ModelChar <- "StanModel"

cat(paste0("Begin Model compilation ", Sys.time(), "\n"))
StanModel <- stan_model(paste0(ModelChar, ".stan"))
cat(paste0("Model compilation done ", Sys.time(), "\n"))

# Create and write meta data

ModelMetaData 				= c()
ModelMetaData$iter 			= 10000 #Increase
ModelMetaData$warmup 		= 2500 #Increase
ModelMetaData$thin 			= 1
ModelMetaData$chains 		= 10 #Increase
ModelMetaData$adapt_delta 	= 0.9
ModelMetaData$max_treedepth = 15
ModelMetaData$ModelChar 	= ModelChar

ModelMetaData_dummy = as.data.frame(unlist(ModelMetaData))
colnames(ModelMetaData_dummy) = NULL


#### Run ####

cat(paste0("Begin Sampling ", Sys.time(), "\n"))
fit = sampling(StanModel, data = data_stan, 
                 iter 	= ModelMetaData$iter, 
                 warmup 	= ModelMetaData$warmup, 
                 thin 	= ModelMetaData$thin, 
                 chains 	= ModelMetaData$chains, 
                 pars 	= c("VaxEffect",
                           "LogPredictions", "RegionalTrends",
                           "NationalTrend", "RegionalScale",
                           "VarAdvantage"), 
                 control = list(adapt_delta = ModelMetaData$adapt_delta, max_treedepth = ModelMetaData$max_treedepth))
cat(paste0("Sampling done ", Sys.time(), "\n"))



# RESULTS -----------------------------------------------------------------


#### Get pars of interest

## LOO-CV

loo_cv <- loo(fit)

## Fit results

model_matrix <- as.matrix(fit)


# Vaccine Effect

VE <- colMeans(model_matrix[, grep("VaxEffect", colnames(model_matrix))])

VE_Quan <- colQuantiles(model_matrix[, grep("VaxEffect", colnames(model_matrix))], probs=c(0.025,0.975))

VE_data <- round(data.frame(VE, VE_Quan), digits = 4)

# Var Advantage

if (data_stan[[5]] == 1) {
  
  VarAdvantage <- colMeans(model_matrix[, grep('^VarAdvantage\\[', colnames(model_matrix))])
  
  VarAdvantage_Quan <- colQuantiles(model_matrix[, grep('^VarAdvantage\\[', colnames(model_matrix))], probs=c(0.025,0.975))
  
  VarAdvantage_data <- round(data.frame(VarAdvantage, VarAdvantage_Quan), digits = 4)
}

# Rt predictions

Rt <- colMeans(model_matrix[, grep("LogPredictions", colnames(model_matrix))])

Rt_Quan <- colQuantiles(model_matrix[, grep("LogPredictions", colnames(model_matrix))], probs=c(0.025,0.975))

Rt_data <- round(data.frame(Rt, Rt_Quan), digits = 4)

# Regional Trends

RegionalTrends <- colMeans(model_matrix[, grep("RegionalTrends", colnames(model_matrix))])

RegionalTrends_Quan <- colQuantiles(model_matrix[, grep("RegionalTrends", colnames(model_matrix))], probs=c(0.025,0.975))

RegionalTrends_data <- round(data.frame(RegionalTrends, RegionalTrends_Quan), digits = 4)

# National Trend

NationalTrend <- colMeans(model_matrix[, grep('^NationalTrend\\[', colnames(model_matrix))])

NationalTrend_Quan <- colQuantiles(model_matrix[, grep('^NationalTrend\\[', colnames(model_matrix))], probs=c(0.025,0.975))

NationalTrend_data <- round(data.frame(NationalTrend, NationalTrend_Quan), digits = 4)

# Scaling

Scaling <- colMeans(model_matrix[, grep('^gamma\\[', colnames(model_matrix))])

Scaling_Quan <- colQuantiles(model_matrix[, grep('^gamma\\[', colnames(model_matrix))], probs=c(0.025,0.975))

Scaling_data <- round(data.frame(Scaling, Scaling_Quan), digits = 4)

## Return object

list_result <- list(loo = loo_cv,
                    VaccineEffect = VE_data,
                    if(data_stan[[5]] == 1) {VarAdvantage = VarAdvantage_data} else {VarAdvantage = NA},
                    Rt_Predictions = Rt_data,
                    RegionalRends = RegionalTrends_data,
                    NationalTrend = NationalTrend_data,
                    Scaling = Scaling_data)

list_result


#### Save data ####

model_name <- "Main_Model"

saveRDS(list_result, paste0("Results/Final_for_paper_LOCAL_", model_name, ".Rds"))
