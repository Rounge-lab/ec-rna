#### R script related to differential expression analysis in the manuscript "Differential levels of circulating RNA prior to endometrial cancer diagnosis"
#### Date: 18 October 2023

### Aim: Using Limma Bioconductor package for differential expression (DE) analysis of RNA data
### Models: presented in Table 2 in the manuscript (Models 1, 2 and 3)
### Random effect in the analysis: sample batch groups (3 groups)
### Note: the script presented here is shown for Model 1 where the effect of case-control status and binary BMI levels (low BMI and high BMI) is evaluated.
### Note: for Models 2 and 3, the procedure is very similar given that smoking and physical activity were also recoded as binary variables.
### Note: below "cacostat" refers to case-control status, "bmibin" is a binary BMI variable ('low BMI' and 'high BMI'), "age_sample_collection" refers to age at sample collection as a scaled and centred variable, and "bd_grp" refers to sample batch groups.

## R/Bioconductor packages needed for analyses and visualization
library(edgeR) # will load Limma as well
library(ggplot2)
library(tidyverse)

## load meta file and raw RNA count file
meta <- readRDS("meta.rds")
rna <- readRDS("rna.rds")

## make sure that row names and column names match between the meta and RNA dataframes (i.e. samples as row in meta data frame should exactly match samples as column in RNA data frame) 
all(rownames(meta) %in% colnames(rna)) # should report TRUE
all(rownames(meta) == colnames(rna)) # should report TRUE

## introducing a single factor for case-control status for low BMI and high BMI
Group <- factor(paste(meta$cacostat, meta$bmibin, sep = "_"))

## creating the model matrix
design <- model.matrix(~0 + age_sample-collection + Group, meta)

## creating the DGEList object for edgeR
y <- DGEList(counts = rna, samples = meta)

## pre-filtering of RNA (removes RNA that have very low counts in most samples)
keep <- filterByExpr(y, y$samples$cacostat)
y <- y[keep,]
dim(y) # check how many RNA remained

## calculation of normalization factors to normalize for RNA composition 
y <- calcNormFactors(y)
v <- voom(y, design, plot = TRUE)  # visually check the voom plot

## apply duplicateCorrelation for the random effect (bd_grp)
corfit <- duplicateCorrelation(v, design, block = meta$bd_grp)
corfit$consensus.correlation   # check the value of consensus.correlation to evaluate if having a mixed model is advantageous 

## apply voom again (with the block and correlation parameters this time)
v <- voom(y, design, block = meta$bd_grp, correlation = corfit$consensus)
# note: above, we once used voom to estimate v, so we can used it to get corfit. Next time, we calculated v again, with the corfit consensus.

## apply lmFit to fit the linear model
fit <- lmFit(v, design, block = meta$bd_grp, correlation = corfit$consensus)

## construct the contrast matrix based on coefficients form the linear model fit
# in the manuscript, these are referred to, respectively, as "case-control in low BMI", "case-control in high BMI", "case-control for low BMI and high BMI", and "difference in case-control between low BMI and high BMI (interaction)")
contrasts <- makeContrasts(
  CasevsControlforLowbmi = GroupCase_lowbmi - GroupControl_lowbmi,
  CasevsControlforHighbmi = GroupCase_highbmi - GroupControl_highbmi,
  CasevsControl = (GroupCase_highbmi + GroupCase_lowbmi)/2 - (GroupControl_highbmi + GroupControl_lowbmi)/2,
  CasevsControldiff = (GroupCase_highbmi - GroupCase_lowbmi) - (GroupControl_highbmi - GroupControl_lowbmi),
  levels = design)

## compute estimated coefficients and standard errors for the defined contrasts
fit2 <- contrasts.fit(fit, contrasts)

## compute t-statistics and p-values from coefficients and standard errors
fit2 <- eBayes(fit2)

## extract tables of RNA from the linear model fit
topTableCaseConLowbmi <- topTable(fit2, coef = "CasevsControlforLowbmi", sort = "P", n = Inf)
topTableCaseConHighbmi <- topTable(fit2, coef = "CasevsControlforHighbmi", sort = "P", n = Inf)
topTableCaseCon <- topTable(fit2, coef = "CasevsControl", sort = "P", n = Inf)
topTableCaseCondiff <- topTable(fit2, coef = "CasevsControldiff", sort = "P", n = Inf)

## vizualising the associations using volcano plot with a different color for RNA that were statistically differentially expressed
# here shown for case-control for low BMI and high BMI
sum(topTableCaseCon$adj.P.Val < 0.05, na.rm = TRUE) # reports how many RNA were statistically differentially expressed
volcanoTab.modelCaseCon  <- topTableCaseCon %>%
                        rownames_to_column("mirna") %>%
                        mutate(`-log10(adjpvalue)` = -log10(adj.P.Val)) %>%
                        mutate(diff_expressed = adj.P.Val < 0.05)
ggplot(volcanoTab.modelCaseCon, aes(x = logFC, y = `-log10(adjpvalue)`,color = diff_expressed)) +
  geom_point() +
  labs(title = "case-control for low BMI and high BMI")




