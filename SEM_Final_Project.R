### DATA PREPARATIONS ### 
psycho <- read.spss("Conradi et al (2015) Data.sav", to.data.frame=TRUE, use.value.labels=FALSE)
colnames(psycho)[44:96] <- c(paste0("Q", 1:50), "Q23r", "Q35r", "Q49r")
colnames(psycho)[97:106] <- paste0(c("DC", "GR", "LY", "MA", "REM", "UNEMO", "CAL", 
                                     "THRILL", "IMPULS", "IRRESP"), "_summed")
length(which(is.na(psycho)))

#Select variables
psycho <- psycho[,c(44:106, 125, 127)]

# split the data into training and test
set.seed(666) # because otherwise this will differ

psycho_training <- data.frame(splitSample(psycho)[[1]])
psycho_testing <- data.frame(splitSample(psycho)[[2]])
View(psycho)
length(psycho$Geslacht[psycho == "Vrouw"])
# calculate the attachment style for each person 
# this has to be done separately for training and testing sets because otherwise our training dataset 
# is not unbiased, so first for training:
med_avoi = median(psycho_training$AVOIm)
med_anxi = median(psycho_training$ANXIm)

for (i in 1:nrow(psycho_training)) {
  
  if (psycho_training[i, "AVOIm"] >= med_avoi &
      psycho_training[i, "ANXIm"] >= med_anxi) {
    psycho_training[i, "attachment"] = "high-high"
  } else if (psycho_training[i, "AVOIm"] >= med_avoi &
             psycho_training[i, "ANXIm"] < med_anxi) {
    psycho_training[i, "attachment"] = "high-low"
  } else if (psycho_training[i, "AVOIm"] < med_avoi &
             psycho_training[i, "ANXIm"] >= med_anxi) {
    psycho_training[i, "attachment"] = "low-high"
  } else if (psycho_training[i, "AVOIm"] < med_avoi &
             psycho_training[i, "ANXIm"] < med_anxi) {
    psycho_training[i, "attachment"] = "low-low"
  }
  
}

# now for testing:
med_avoi = median(psycho_testing$AVOIm)
med_anxi = median(psycho_testing$ANXIm)

for (i in 1:nrow(psycho_testing)) {
  
  if (psycho_testing[i, "AVOIm"] >= med_avoi &
      psycho_testing[i, "ANXIm"] >= med_anxi) {
    psycho_testing[i, "attachment"] = "high-high"
  } else if (psycho_testing[i, "AVOIm"] >= med_avoi &
             psycho_testing[i, "ANXIm"] < med_anxi) {
    psycho_testing[i, "attachment"] = "high-low"
  } else if (psycho_testing[i, "AVOIm"] < med_avoi &
             psycho_testing[i, "ANXIm"] >= med_anxi) {
    psycho_testing[i, "attachment"] = "low-high"
  } else if (psycho_testing[i, "AVOIm"] < med_avoi &
             psycho_testing[i, "ANXIm"] < med_anxi) {
    psycho_testing[i, "attachment"] = "low-low"
  }
  
}


colnames(psycho_training) <- c(paste0("Q", 1:50), "Q23r", "Q35r", "Q49r", 
                               "DC_summed", "GR_summed", "LY_summed", "MA_summed", 
                               "REM_summed", "UNEMO_summed", "CAL_summed", 
                               "THRILL_summed", "IMPULS_summed", "IRRESP_summed",
                               "AVOIm", "ANXIm", "attachment")
colnames(psycho_testing) <- c(paste0("Q", 1:50), "Q23r", "Q35r", "Q49r", 
                              "DC_summed", "GR_summed", "LY_summed", "MA_summed", 
                              "REM_summed", "UNEMO_summed", "CAL_summed", 
                              "THRILL_summed", "IMPULS_summed", "IRRESP_summed",
                              "AVOIm", "ANXIm", "attachment")

### 3-FACTOR MODEL WITH SUM SCORES ### 

mod_three_sum <- '
  # Factor loadings:
  GM =~ DC_summed + GR_summed + LY_summed + MA_summed
  CU =~ UNEMO_summed + CAL_summed + REM_summed
  II =~ THRILL_summed + IMPULS_summed + IRRESP_summed
'
configural_three_sum <- cfa(mod_three_sum, psycho_training)
fitMeasures(configural_three_sum, c("chisq", "df", "pvalue"))
fitMeasures(configural_three_sum, c("RMSEA",
                                    "CFI", "RNI", "TLI",
                                    "GFI", "AGFI"))
# plot this using the test data
configural_three_sum_test <- cfa(mod_three_sum, psycho_testing)
semPaths(configural_three_sum_test,
         what = "col", whatLabels = "est",
         layout = "tree2", curvePivot = TRUE, groups = "latents", pastel = TRUE,
         title = FALSE)

modificationindices(configural_three_sum) %>% arrange(-mi)

### try to modify three_sum because the RMSEA is so big (Sacha's recommendation)

modified_again_three_sum <- '
  # Factor loadings:
  GM =~ DC_summed + GR_summed + LY_summed + MA_summed + REM_summed # added REM
  CU =~ UNEMO_summed + CAL_summed + REM_summed + GR_summed + IMPULS_summed # added GR and IMPULS
  II =~ THRILL_summed + IMPULS_summed + IRRESP_summed + MA_summed + CAL_summed  # added MA and CAL
  
  #residual variances
  UNEMO_summed ~~ CAL_summed
  UNEMO_summed ~~ REM_summed
'
configural_modified_again <- cfa(modified_again_three_sum, psycho_testing)
fitMeasures(configural_modified_again, c("chisq", "df", "pvalue"))
fitMeasures(configural_modified_again, c("RMSEA",
                                         "CFI", "RNI", "TLI",
                                         "GFI", "AGFI"))
#p value is signficant, but it might be overpowered --> fit indices indicate good fit

# HT of close fit
DFm <- 25
Tm <- 58.391
n <- nrow(psycho_testing)
lambda_c <- 0.05^2 * n * DFm
pchisq(Tm, DFm, lambda_c, lower.tail = FALSE)
## [1] 0.4782273
### could not reject the hypothesis that the model fits well

# HT of non-close fit
lambda_c <- 0.08^2 * 500 * DFm
pchisq(Tm, DFm, lambda_c, lower.tail = TRUE)
##  [1] 0.002923551
### could indeed reject hypothesis that model fits poorly

## use this for means and variances 

# fit configural model with attachment 
conf <- cfa(modified_again_three_sum, psycho_testing, group = "attachment")
fitMeasures(conf, c("chisq", "df", "pvalue"))
fitMeasures(conf, c("RMSEA","CFI", "RNI", "TLI","GFI", "AGFI"))
#p value is significant, but may be overpowered --> fit indices indicate good fit

# weak invariance testing
weak <- cfa(modified_again_three_sum, psycho_testing, group = "attachment", 
            group.equal = "loadings")
anova(conf, weak)

# this is not significant and it has a lower AIC & BIC, so we can go on with the strong model 

strong <- cfa(modified_again_three_sum, psycho_testing, group = "attachment",
              group.equal = c("loadings", "intercepts"))
anova(weak, strong)

#Strong invariance does not hold --> check for partial strong invariance

# lagrange multiplier test:
lavTestScore(strong)$uni %>%
  arrange(-X2)
parameterEstimates(strong)[7,]

strong2 <- cfa(modified_again_three_sum, psycho_testing, group = "attachment",
             group.equal = c("loadings", "intercepts"),
             group.partial = "CU =~ CAL_summed")
anova(weak, strong2)
#It holds, so go on with strict invariance

strict <- cfa(modified_again_three_sum, psycho_testing, group = "attachment",
              group.equal = c("loadings", "intercepts", "residuals", "residual.covariances"),
              group.partial = "CU =~ CAL_summed")

anova(strong2, strict)
#It does not hold - include more group.partial parameters
# lagrange multiplier test:
lavTestScore(strict)$uni %>%
  arrange(-X2)
parameterEstimates(strict)[c(17), ]

strict2 <- cfa(modified_again_three_sum, psycho_testing, group = "attachment",
              group.equal = c("loadings", "intercepts", "residuals", "residual.covariances"),
              group.partial = c("CU =~ CAL_summed", "REM_summed ~~ UNEMO_summed"))

anova(strong2, strict2)


semPlot <- semPaths(strict2, rotation = 2,
         what = "col", whatLabels = "est",
         layout = "tree2", curvePivot = TRUE, groups = "latents", pastel = TRUE,
         title = FALSE)

#save Paths
pdf(file = "SEMPlot_Project.pdf")
semPaths(strict2, 
         what = "col", whatLabels = "est",
         layout = "tree2", curvePivot = TRUE, groups = "latents", pastel = TRUE,
         title = FALSE)
dev.off()

# strict2 shows that partial invariance holds, as the test is non-significant

# means and variances 
# test for homogeneity 
# do the variances differ?
eqvars <- cfa(modified_again_three_sum, psycho_testing, group = "attachment",
              group.equal = c("loadings", "intercepts", "residuals",
                              "residual.covariances", "lv.variances",
                              "lv.covariances"),
              group.partial = c("CU =~ CAL_summed", "REM_summed ~~ UNEMO_summed"))
anova(strict2, eqvars)

# it is  significant, meaning that meaning that groups differ in variances 

# do the means differ?
eqmeans <- cfa(modified_again_three_sum, psycho_testing, group = "attachment",
               group.equal = c("loadings", "intercepts", "residuals",
                               "residual.covariances", "lv.variances",
                               "lv.covariances", "means"),
               group.partial = c("CU =~ CAL_summed", "REM_summed ~~ UNEMO_summed"))
anova(strict2, eqmeans)

# it is significant, so inspect which parameters could be freed

# lagrange multiplier test:
lavTestScore(eqmeans)$uni %>%
  arrange(-X2)
parameterEstimates(eqmeans)[c(23), ]

eqmeans2 <- cfa(modified_again_three_sum, psycho_testing, group = "attachment",
                group.equal = c("loadings", "intercepts", "residuals",
                                "residual.covariances", "lv.variances",
                                "lv.covariances", "means"),
                group.partial = c( "CU =~ CAL_summed", "REM_summed ~~ UNEMO_summed",
                                   "UNEMO_summed ~~ UNEMO_summed"))
anova(strict2, eqmeans2)

# it is  significant, meaning that groups differ in means 
# We can inspect the mean differences:

parameterestimates(eqmeans2) %>%
  filter(op == "~1") %>%
  arrange(-est)

#The means for Thrillseeking and Impulsivity and Dishonest Charm seem to differ 
#the most between the groups (arguably, also all other observed variable means differ across groups)

# does the model still fit? 
fitMeasures(eqmeans,
            c("chisq", "df", "pvalue", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi"))
#yes, it does, it fits very well indeed... 
