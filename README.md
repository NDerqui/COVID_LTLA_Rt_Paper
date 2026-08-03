# Modelling the impact of vaccination on SARS-CoV-2 transmission in England

This repository contains the code for data cleaning, analysis and figures generation used for the manuscript "Modelling the impact of vaccination on SARS-CoV-2 transmission in England", currently under review but whose early version is available online as a [preprint](https://www.medrxiv.org/content/10.1101/2024.09.06.24313210v1).

## Abstract

**Background** Efficacy and effectiveness of vaccines against SARS-CoV-2 infection, severe disease and death have been widely assessed. However, the impact of vaccination against SARS-CoV-2 transmission is far less well-characterized, and has major implications for public health, because it informs the indirect effects of vaccination in addition to its direct effects. Analysing the effects of SARS-CoV-2 vaccination on transmission is challenging, because they must be considered in tandem with the time-varying reproduction number (Rt), while also accounting for regional variability, for example due to the presence of more transmissible variants.

**Methods** We fitted a Bayesian hierarchical model to previously obtained estimates of Rt to estimate the effectiveness of vaccination with one, two and three doses on SARS-CoV-2 transmission in the UK during 2021. Vaccine effectiveness is defined as the proportional reduction in the time-varying reproduction number Rt. The model accounts for transmission at national and Lower Tier Local Authority (LTLA)-level, and uses vaccination data provided by the UK Health Security Agency (UKHSA), detailing the LTLA-specific proportions of people who have received doses one, two and three. The model also incorporates data on the proportion of wild-type, Alpha and Delta SARS-CoV-2 variants over time in each LTLA, obtained from UKHSA and the COVID-19 Genomics UK (COG-UK) Consortium.

**Results** We find that vaccination had moderate-to-large effectiveness against transmission for dose 1 (39.30%, 95% CrI 26.64% - 48.07%), and for dose 3 (48.69%, 95% CrI 27.97% - 71.30%), but negligible effects on dose 2, likely attributable to the coincident importation and dominance of the Delta variant in the UK. Nationally, our model fitted the previously estimated values of time-series of Rt values well, largely reproducing the reproduction number averaged across LTLAs for each timepoint. This lends support to our hypothesis that the extent of vaccination (or lack thereof) was a major determinant of transmission intensity. Our model fits further reproduced well the reproduction numbers at regional level, although outliers were less well captured, implying some degree of variation that is not explained by our model.

**Conclusions** To our knowledge, our analysis is the first evidence of the effectiveness of SARS-CoV-2 vaccination against its transmission at population level. We find that vaccination is an effective tool for the control of SARS-CoV-2 transmission, in addition to its well-documented effects on disease burden and mortality. Our results allow future assessment of the impact of vaccination accounting for several circulating variants and sociodemographic factors.

## Citation

If you would like to cite this work, please use:

```
@article{
   author = {Derqui, Nieves and Mishra, Swapnil and Hinsley, Wes R. and Bhatt, Samir and Laydon, Daniel J.},
   title = {Impact of vaccination on SARS-CoV-2 transmission in the UK: a modelling study},
   journal = {medRxiv},
   DOI = {10.1101/2024.09.06.24313210},
   url = {https://www.medrxiv.org/content/10.1101/2024.09.06.24313210v1},
   year = {2024},
   type = {Journal Article}
}
```
