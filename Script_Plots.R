#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Impact of vaccination on SARS-CoV-2 transmission in the UK: a modelling study #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


## Script plot model results
# 1: Retrieve saved list
# 2: Plot!



# SET UP ------------------------------------------------------------------

# Packages

library(tidyverse)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(rcartocolor)
library(here)



# DATA --------------------------------------------------------------------


#### 'Observed' data ####

# Data imputed into the model
# Saved in session from previous script

Rt_data <- data_stan[[17]]

LTLA <- data_stan[[15]]

VaxProp <- as.data.frame(data_stan[[18]])

VarProp <- data_stan[[19]]

NamesLTLAs <- data.frame(NamesLTLAs = data_stan[[21]])
NamesLTLAs$LTLA <- 1:221

date <- max(c(min(data_rt$date), min(data_var$date))) + data_stan[[14]]*7


#### Model data ####

## Main Model: with Variants, but not Variant's VE

list_result <- readRDS("Results/Final_for_paper_LOCAL_Main_Model.Rds")


#### Substract parameters ####

DoVariants <- 1

DoAge <- 0

# Parameters

VaxEffect_data <- list_result[[2]]
VaxEffect <- list_result[[2]][1]

if (DoVariants == 1) {
  VarAdvantage_data <- list_result[[3]]
  VarAdvantage_data <- rbind(c(1, 1, 1), VarAdvantage_data)
  VarAdvantage <- list_result[[3]][1]
  VarAdvantage <- rbind(1, VarAdvantage)
} else {
  VarAdvantage_data <- data.frame(Null = NA)
}

Rt_Predictions_data <- list_result[[4]]
Rt_Predictions <- list_result[[4]][1]
colnames(Rt_Predictions_data) <- c("Rt", "Rt_low", "Rt_upp")

RegionalTrends_data <- list_result[[5]]
RegionalTrends <- list_result[[5]][1]
colnames(RegionalTrends_data) <- c("RegionalTrends", "RegionalTrends_low", "RegionalTrends_upp")
RegionalTrends_low <- list_result[[5]][2]
RegionalTrends_upp <- list_result[[5]][3]

NationalTrend_data <- list_result[[6]]
NationalTrend <- list_result[[6]][1]

Scaling_data <- list_result[[7]]
Scaling <- list_result[[7]][1]

Rt_NoVax <- (as.matrix(VarProp) %*% as.matrix(VarAdvantage))*RegionalTrends
Rt_NoVax_low <- (as.matrix(VarProp) %*% as.matrix(VarAdvantage))*RegionalTrends_low
Rt_NoVax_upp <- (as.matrix(VarProp) %*% as.matrix(VarAdvantage))*RegionalTrends_upp
Rt_NoVax_data <- data.frame(Rt_NoVax, Rt_NoVax_low, Rt_NoVax_upp)
colnames(Rt_NoVax_data) <- c("Rt_NoVax", "Rt_NoVax_low", "Rt_NoVax_upp")

estimates <- data.frame(LTLA, date, Rt_data, Rt_Predictions_data,
                        Rt_NoVax_data, RegionalTrends_data, NationalTrend_data)
estimates <- merge(estimates, NamesLTLAs)


#### VaxEffect ####

# Change rownames for easier interpretation

if (nrow(VaxEffect_data) == 3) {
  row.names(VaxEffect_data) <- c("Dose 1", "Dose 2", "Dose 3")
  
} else {
  
  if(nrow(VaxEffect_data) == 6) {
    row.names(VaxEffect_data) <- c("Alpha1", "Alpha2", "Alpha3",
                                   "Delta1", "Delta2", "Delta3")
  } else {
    
    if (DoAge == 1) {
      row.names(VaxEffect_data) <- c("15-49_D1", "15-49_D2", "15-49_D3",
                                     "50-69_D1", "50-69_D2", "50-69_D3",
                                     "70+_D1", "70+_D2", "70+_D3")
    } else {
      row.names(VaxEffect_data) <- c("PreAl_1", "PreAl_2", "PreAl_3",
                                     "Alpha1", "Alpha2", "Alpha3",
                                     "Delta1", "Delta2", "Delta3")
    }
  }
}


#### Save ####

list_table <- list("VaxEffect" = VaxEffect_data,
                   "VarAdvantage" = VarAdvantage_data,
                   "Scaling" = Scaling_data,
                   "NationalTrends" = NationalTrend_data,
                   "RegionalTrends" = RegionalTrends_data,
                   "Rt_NoVax_[VarAd_x_RegTrend]" = Rt_NoVax,
                   "RtPredictions" = Rt_Predictions_data,
                   "RtData" = Rt_data)



# MINOR -------------------------------------------------------------------


## How much does vaccination reduce Rt overall

mean(estimates$Rt)
mean(estimates$Rt[estimates$date < as.Date("30/04/2021", format = "%d/%m/%Y")])
mean(estimates$Rt[estimates$date < as.Date("15/07/2021", format = "%d/%m/%Y")])

mean(estimates$Rt_NoVax)
mean(estimates$Rt_NoVax[estimates$date < as.Date("30/04/2021", format = "%d/%m/%Y")])
mean(estimates$Rt_NoVax[estimates$date < as.Date("15/07/2021", format = "%d/%m/%Y")])

1 - mean(estimates$Rt)/mean(estimates$Rt_NoVax)
1 - mean(estimates$Rt[estimates$date < as.Date("30/04/2021", format = "%d/%m/%Y")])/mean(estimates$Rt_NoVax[estimates$date < as.Date("30/04/2021", format = "%d/%m/%Y")])
1 - mean(estimates$Rt[estimates$date < as.Date("15/07/2021", format = "%d/%m/%Y")])/mean(estimates$Rt_NoVax[estimates$date < as.Date("15/07/2021", format = "%d/%m/%Y")])



# PLOTS -------------------------------------------------------------------


#### Figure 1 ####

# Observed Rt

fig_1a <- ggplot(data = estimates,
                 mapping = aes(x = date, y = Rt_data, group = date)) +
          geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
          geom_boxplot(color = carto_pal(name = "Safe")[10],
                       fill = carto_pal(name = "Safe")[10], alpha = 0.2) +
          scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
          theme_bw() +
          labs(y = "Observed Reproduction Number (Rt)") +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
                axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.05)),
                legend.text = element_text(size = rel(1.1)))

# Vaccination data

fig_1b <- ggplot(data = estimates,
                 mapping = aes(x = date, group = date)) +
          geom_boxplot(data = VaxProp, mapping = aes(y = V1, fill = "V1", color = "V1"), alpha = 0.2) +
          geom_boxplot(data = VaxProp, mapping = aes(y = V2, fill = "V2", color = "V2"), alpha = 0.2) +
          geom_boxplot(data = VaxProp, mapping = aes(y = V3, fill = "V3", color = "V3"), alpha = 0.2) +
          scale_fill_manual(breaks = c("V1", "V2", "V3"),
                            values = c(carto_pal(name = "Safe")[4],
                                       carto_pal(name = "Safe")[8],
                                       carto_pal(name = "Safe")[3]),
                            labels = c("Dose 1", "Dose 2", "Dose 3"),
                            name = "Dose") +
          scale_color_manual(breaks = c("V1", "V2", "V3"),
                             values = c(carto_pal(name = "Safe")[4],
                                        carto_pal(name = "Safe")[8],
                                        carto_pal(name = "Safe")[3]),
                            labels = c("Dose 1", "Dose 2", "Dose 3"),
                            name = "Dose") +
          scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
          scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
          theme_bw() +
          labs(y = "Proportion of vaccinated population") +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
                axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.05)),
                legend.text = element_text(size = rel(1.1)), 
                legend.position = "bottom",
                legend.title = element_blank())

# Variants data

VarProp_plot <- data.frame(date = estimates$date, VarProp) %>%
  pivot_longer(-date, names_to = "variant", values_to = "prop")

fig_1c <- ggplot(data = VarProp_plot,
                 mapping = aes(x = date, y = prop, fill = variant)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(breaks = c("X1", "X2", "X3"),
                    values = c(carto_pal(name = "Safe")[4],
                               carto_pal(name = "Safe")[8],
                               carto_pal(name = "Safe")[3]),
                    labels = c("Wild-type", "Alpha", "Delta")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_bw() +
  labs(y = "Proportion of circulating variants") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.05)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "bottom",
        legend.title = element_blank())

### All

# fig_1 <- ggarrange(fig_1a, fig_1b, fig_1c, ncol = 1, labels = c("A", "B", "C"))
fig_1 <- ggarrange(fig_1b, fig_1c, ncol = 1, labels = c("A", "B"))

png(filename = "Results/Figure1.png",
    height = 10, width = 8, res = 1200, units = "in")
fig_1
dev.off()


#### Figure 2 ####

# Observed vs predicted Rt

fig_2 <- ggplot(data = estimates,
                  mapping = aes(x = date, group = date)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
  geom_boxplot(mapping = aes(y = Rt_data, color = "RtData", fill = "RtData"), alpha = 0.2) +
  geom_boxplot(mapping = aes(y = Rt, color = "RtLogP", fill = "RtLogP"), alpha = 0.2) +
  scale_color_manual(breaks = c("RtLogP", "RtData"),
                     values = c(carto_pal(name = "Safe")[4], carto_pal(name = "Safe")[10]),
                     labels = c("Model-Predicted Rt", "Observed Rt"), name = "Rt") +
  scale_fill_manual(breaks = c("RtLogP", "RtData"),
                    values = c(carto_pal(name = "Safe")[4], carto_pal(name = "Safe")[10]),
                    labels = c("Model-Predicted Rt", "Observed Rt"), name = "Rt") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  scale_y_continuous(limits = c(0.0, 3)) +
  theme_bw() +
  labs(y = "Reproduction Number (Rt)") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank())

### All

png(filename = "Results/Figure2.png",
    height = 5, width = 8, res = 1200, units = "in")
fig_2
dev.off()


#### Figure 3 ####

## Parameter in a couple of LTLAs

estimates_sum <- estimates %>%
  filter(LTLA %in% c(1,47,93,104,122,147,208,179,221))

fig_3 <- ggplot(data = estimates_sum) +
  geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
  geom_line (mapping = aes(x = date, y = Rt_NoVax, color = "WithVarAd"), linewidth = 1.03) +
  geom_ribbon(mapping = aes(x = date, ymin = Rt_NoVax_low, ymax = Rt_NoVax_upp, fill = "WithVarAd"), alpha = 0.3) +
  geom_line (mapping = aes(x = date, y = Rt, color = "RtLogP"), linewidth = 1.03) +
  geom_ribbon(mapping = aes(x = date, ymin = Rt_low, ymax = Rt_upp, fill = "RtLogP"), alpha = 0.3) +
  geom_line (mapping = aes(x = date, y = Rt_data, color = "RtData"), linewidth = 1.03) +
  scale_fill_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                    values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                               carto_pal(name = "Safe")[10]),
                    labels = c("Rt in the absence of vaccination",
                               "Model-Predicted Rt", "Observed Rt"),
                    name = "Type") +
  scale_color_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                     values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                              carto_pal(name = "Safe")[10]),
                     labels = c("Rt in the absence of vaccination",
                                "Model-Predicted Rt", "Observed Rt"),
                     name = "Type") +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y") +
  scale_y_continuous(limits = c(0.0, 3)) +
  labs(y = "Reproduction Number (Rt)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.05), angle = 10), axis.text.y = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)), strip.text = element_text(size = rel(1.2)),
        legend.position = "bottom", legend.title = element_blank()) +
  facet_wrap(NamesLTLAs ~ .) + guides(fill = "none")

png(filename = "Results/Figure3.png",
    height = 8, width = 12, res = 1200, units = "in")
fig_3
dev.off()


#### Figure 4 ####

# VE

VaxEffect_plot <- cbind(dose = c("Dose 1", "Dose 2", "Dose 3"),
                        VaxEffect_data[1:3,])

fig_4a <- ggplot(data = VaxEffect_plot,
                 mapping = aes(x = dose, y = VE, color = dose)) +
  geom_point(shape = 15, size = 3) +
  geom_errorbar(mapping = aes(ymin = X2.5., ymax = X97.5.), width = 0.1) +
  scale_color_manual(values = c(carto_pal(name = "Safe")[4],
                                carto_pal(name = "Safe")[8],
                                carto_pal(name = "Safe")[3]),
                     labels = c("Dose 1", "Dose 2", "Dose 3")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  theme_bw() +
  labs(y = "Vaccine Effectiveness") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

# VarAd

VarAdvantage_plot <- cbind(variant = c(" Pre-Alpha", "Alpha", "Delta"),
                           VarAdvantage_data)

fig_4b <- ggplot(data = VarAdvantage_plot,
                 mapping = aes(x = variant, y = VarAdvantage, color = variant)) +
  geom_point(shape = 15, size = 3) +
  geom_errorbar(mapping = aes(ymin = X2.5., ymax = X97.5.), width = 0.1) +
  scale_x_discrete(breaks = c(" Pre-Alpha", "Alpha", "Delta"),
                   labels = c("Wild-type", "Alpha", "Delta")) +
  scale_color_manual(breaks = c(" Pre-Alpha", "Alpha", "Delta"),
                     values = c(carto_pal(name = "Safe")[4],
                                carto_pal(name = "Safe")[8],
                                carto_pal(name = "Safe")[3]),
                     labels = c("Wild-type", "Alpha", "Delta")) +
  theme_bw() +
  labs(y = "Variant Advantage") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

## All

fig_4 <- ggarrange(fig_4a, fig_4b, ncol = 2, labels = c("A", "B"))

png(filename = "Results/Figure4.png",
    height = 4, width = 6, res = 1200, units = "in")
fig_4
dev.off()


#### Figure 5 ####

# Regional Trend with VarAd

fig_5 <- ggplot(data = estimates,
                 mapping = aes(x = date, group = LTLA)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
  geom_line(mapping = aes(y = Rt_NoVax, color = "WithVarAd", alpha = "Two")) +
  geom_line(mapping = aes(y = Rt_data, color = "Rt_data", alpha = "One")) +
  scale_color_manual(breaks = c("Rt_data", "WithVarAd"),
                     values = c(carto_pal(name = "Safe")[10], carto_pal(name = "Safe")[7]),
                     labels = c("Observed Rt in each LTLA", "Fitted Rt in the absence of vaccination in each LTLA"),
                     name = "Type") +
  scale_alpha_manual(breaks = c("One", "Two"),
                     values = c(0.2, 0.3),
                     labels = c("Observed Rt in each LTLA", "Fitted Rt in the absence of vaccination in each LTLA"),
                     name = "Type") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  scale_y_continuous(limits = c(0.0, 2.7)) +
  theme_bw() +
  labs(y = "Reproduction Number (Rt)") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.1), angle = 10), axis.text.y = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

png(filename = "Results/Figure5.png",
    height = 5, width = 8, res = 1200, units = "in")
fig_5
dev.off()


#### Supplementary 1 ####

# Model vs observed Rt

fig_s1a <- ggplot(data = estimates) +
  geom_point(mapping = aes(x = Rt_data, y = Rt)) +
  geom_abline(x = 0, y = 1, color = carto_pal(name = "Safe")[9]) +
  theme_bw() +
  labs(x = "Observed Reproduction Number (Rt)",
       y = "Model-Predicted Reproduction Number (Rt)") 

png(filename = "Results/FigureS1.png",
    height = 5, width = 5, res = 1200, units = "in")
fig_s1a
dev.off()


#### Supplementary 2 ####

# Fig 3 for the rest of LTLAs

estimates_rest <- estimates %>%
  filter(!(LTLA %in% c(1,47,93,104,122,147,208,179,221)))

cut_offs <- list(c(1:20), c(21:40), c(41:60), c(61:80), c(81:100),
                 c(101:120), c(121:140), c(141:160), c(161:180), c(181:200),
                 c(201:212))

cut_titles <- c("FigureS2a", "FigureS2b", "FigureS2c", "FigureS2d", "FigureS2e",
                "FigureS2f", "FigureS2g", "FigureS2h", "FigureS2i", "FigureS2j",
                "FigureS2k")

for (cut in 1:length(cut_offs)) {
  
  p <- ggplot(data = filter(estimates_rest, LTLA %in% cut_offs[[cut]])) +
    geom_hline(yintercept = 1, linetype = "dashed", color = carto_pal(name = "Safe")[9]) +
    geom_line (mapping = aes(x = date, y = Rt_NoVax, color = "WithVarAd"), linewidth = 1.03) +
    geom_ribbon(mapping = aes(x = date, ymin = Rt_NoVax_low, ymax = Rt_NoVax_upp, fill = "WithVarAd"), alpha = 0.3) +
    geom_line (mapping = aes(x = date, y = Rt, color = "RtLogP"), linewidth = 1.03) +
    geom_ribbon(mapping = aes(x = date, ymin = Rt_low, ymax = Rt_upp, fill = "RtLogP"), alpha = 0.3) +
    geom_line (mapping = aes(x = date, y = Rt_data, color = "RtData"), linewidth = 1.03) +
    scale_fill_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                      values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                                 carto_pal(name = "Safe")[10]),
                      labels = c("Rt in the absence of vaccination",
                                 "Model-Predicted Rt", "Observed Rt"),
                      name = "Type") +
    scale_color_manual(breaks = c("WithVarAd", "RtLogP", "RtData"),
                       values = c(carto_pal(name = "Safe")[7], carto_pal(name = "Safe")[4],
                                  carto_pal(name = "Safe")[10]),
                       labels = c("Rt in the absence of vaccination",
                                  "Model-Predicted Rt", "Observed Rt"),
                       name = "Type") +
    scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y") +
    scale_y_continuous(limits = c(0.0, 3)) +
    labs(y = "Reproduction Number (Rt)") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1.1)),
          axis.text.x = element_text(size = rel(0.85), angle = 10), axis.text.y = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(1.1)), strip.text = element_text(size = rel(1.2)),
          legend.position = "bottom", legend.title = element_blank()) +
    facet_wrap(NamesLTLAs ~ .) + guides(fill = "none")
  
  png(filename = paste0("Results/", cut_titles[cut], ".png"),
      height = 12, width = 16, res = 1200, units = "in")
  print(p)
  dev.off()
  
}
rm(p, cut)


#### Supplementary 3 ####

# Comparing vaccine effect

all_ve <- data.frame()

## Loop to retrieve all

names <- c("Main_Model", "VarVE", "Age", "VarAgeVE")
models <- c("  Main model", " Variant-specific-VE", "Age-model", "Age-variant-specific-VE")

for (variation in 1:length(names)) {
  
  data <- readRDS(paste0("Results/Final_for_paper_LOCAL_", names[variation], ".Rds"))
  
  # Extract VE
  
  add <- data[[2]]
  add$model <- paste0(models[variation])
  
  if (nrow(add) == 3) {
    
    add$dose <- c("Dose 1", "Dose 2", "Dose 3")
    add$variant <- "Non-variant specific"
    add$age_group <- "Non-age-group specific"
    
  } else {
    
    if (nrow(add) == 9) {
      
      add$dose <- rep(c("Dose 1", "Dose 2", "Dose 3"), times = 3)
      add$variant <- c(rep(" Wild-type", 3), rep("Alpha", 3), rep("Delta", 3))
      add$age_group <- "Non-age-group specific"
      
    } else {
      
      add$dose <- rep(c("Dose 1", "Dose 2", "Dose 3"), times = 9)
      add$variant <- rep(c(rep(" Wild-type", 3), rep("Alpha", 3), rep("Delta", 3)), 3)
      add$age_group <- c(rep("15-49", 9), rep("50-69", 9), rep("70-plus", 9))
      
    }
  }
  
  colnames(add) <- c("ve", "low", "upp", "model", "dose", "variant", "age_group")
  all_ve <- rbind(all_ve, add)
}
rm(add, variation, data)

## Plot VE

png(paste0("Results/FigureS3.png"),
    width = 8, height = 5, units = 'in', res = 1200)

ggplot(data = all_ve,mapping = aes(x = dose, y = ve, color = variant, shape = age_group)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(mapping = aes(ymin = low, ymax = upp),
                position = position_dodge(0.5), width = 0.2) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(1:2, 4, 12)], name = "Variant") +
  scale_shape_manual(values = c(0, 2, 5, 19), name = "Age Group") +
  labs(y = "Vaccine Effectiveness", color = "Model") +
  theme_bw() +
  theme(title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  facet_wrap(. ~ model) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

dev.off()


#### Table S1 ####

# Comparing loo-cv

all_loo <- data.frame()

## Loop to retrieve all

names <- c("Main_Model", "VarVE", "Age", "VarAgeVE")
models <- c("  Main model", " Variant-specific-VE", "Age-model", "Age-variant-specific-VE")

for (variation in 1:length(names)) {
  
  data <- readRDS(paste0("Results/Final_for_paper_LOCAL_", names[variation], ".Rds"))
  
  # Extract and add
  
  add <- as.data.frame(data[[1]][5:10])
  add$model <- paste0(models[variation])
  
  all_loo <- rbind(all_loo, add)
}
rm(add, variation, data)

write.csv(all_loo, "Results/TableS1.csv", row.names = FALSE)