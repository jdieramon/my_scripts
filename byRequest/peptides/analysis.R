# Dependencies -----------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("genefilter")


library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(MSnSet.utils)
library(genefilter)
library(MSstats)



# Load functions
source('functions.R')

# Load Data 
prots <- read.csv2("data/T1_JDie_09042025.csv", header=TRUE)



# Data Pre-processing ----------------------------------------------------------

# Extract quant data columns for the analysis

# n. of peptides 
# n. PSM
# iBAQ data

            
target <- c("Protein.IDs", "MS.MS.count", 
            "Peptides.Cold.Liver_LDLRKO1.01", "Peptides.Cold.Liver_LDLRKO1.02", 
            "Peptides.Cold.Liver_LDLRKO2.01", "Peptides.Cold.Liver_LDLRKO2.02", 
            "Peptides.Cold.Liver_LDLRKO3.01", "Peptides.Cold.Liver_LDLRKO3.02", 
            "Peptides.RT.Liver_LDLRKO1.01", "Peptides.RT.Liver_LDLRKO1.02",
            "Peptides.RT.Liver_LDLRKO2.01", "Peptides.RT.Liver_LDLRKO2.02", 
            "Peptides.RT.Liver_LDLRKO3.01", "Peptides.RT.Liver_LDLRKO3.02", 
            "iBAQ.Cold.Liver_LDLRKO1.01", "iBAQ.Cold.Liver_LDLRKO1.02", 
            "iBAQ.Cold.Liver_LDLRKO2.01", "iBAQ.Cold.Liver_LDLRKO2.02", 
            "iBAQ.Cold.Liver_LDLRKO3.01", "iBAQ.Cold.Liver_LDLRKO3.02", 
            "iBAQ.RT.Liver_LDLRKO1.01", "iBAQ.RT.Liver_LDLRKO1.02", 
            "iBAQ.RT.Liver_LDLRKO2.01", "iBAQ.RT.Liver_LDLRKO2.02", 
            "iBAQ.RT.Liver_LDLRKO3.01", "iBAQ.RT.Liver_LDLRKO3.02")    

prots <- prots %>% 
  select(all_of(target))

names_tidy <- str_remove(str_remove(str_remove(names(prots), "iBAQ."), "LDLRK"), ".Liver_")
names(prots) = names_tidy

# Expand rows 
df <- rows_expansion(prots)


# Remove contaminats CON_ & 'decoy' REV_
df$Protein.IDs[which(str_detect(df$Protein.IDs, "CON_"))]
df$Protein.IDs[which(str_detect(df$Protein.IDs, "REV_"))]

df <-  df %>% 
  filter(!str_detect(Protein.IDs, "CON_"),
         !str_detect(Protein.IDs, "REV_"))

# If N.peptide = 1, convert intensities = 0 
df <- df %>% 
  mutate(ColdO1.01 = if_else(Peptides.ColdO1.01 == 1, 0, ColdO1.01), 
         ColdO1.02 = if_else(Peptides.ColdO1.02 == 1, 0, ColdO1.02), 
         ColdO2.01 = if_else(Peptides.ColdO2.01 == 1, 0, ColdO2.01), 
         ColdO2.02 = if_else(Peptides.ColdO2.02 == 1, 0, ColdO2.02), 
         ColdO3.01 = if_else(Peptides.ColdO3.01 == 1, 0, ColdO3.01), 
         ColdO3.02 = if_else(Peptides.ColdO3.02 == 1, 0, ColdO3.02), 
         RTO1.01 = if_else(Peptides.RTO1.01 == 1, 0, RTO1.01), 
         RTO1.02 = if_else(Peptides.RTO1.02 == 1, 0, RTO1.02),
         RTO2.01 = if_else(Peptides.RTO2.01 == 1, 0, RTO2.01), 
         RTO2.02 = if_else(Peptides.RTO2.02 == 1, 0, RTO2.02),
         RTO3.01 = if_else(Peptides.RTO3.01 == 1, 0, RTO3.01), 
         RTO3.02 = if_else(Peptides.RTO3.02 == 1, 0, RTO3.02)) 



# Simplify the dataset 
df <- df %>% 
  select(Protein.IDs, MS.MS.count, ColdO1.01:RTO3.02)


# Exclusivity filter ----------------------------------------------------------
treat_col1 = "ColdO1.01"
treat_col2 = "ColdO3.02"
control_col1 = "RTO1.01"
control_col2 = "RTO3.02"


## Exclusive 'treatment' condition
ids_exc.treat <- exclusive_treat(df, treat_col1, treat_col2, control_col1, control_col2)

# sanity check 

df %>% 
  filter(Protein.IDs %in% ids_exc.treat) %>% View()

df %>% 
  filter(Protein.IDs %in% ids_exc.treat) %>% 
  select(3:14) %>% 
  transmute(nsamples_cold = rowSums(across(1:6, ~ as.integer(. > 0))), 
            nsamples_treat = rowSums(across(7:12, ~ as.integer(. > 0)))) %>% 
  View()

# Create the exclusive 'treatment' condition object and save it to a file
df_treat <- df %>% 
  filter(Protein.IDs %in% ids_exc.treat)

#write.csv(df_treat, file = "res/df_treat", row.names = F)
#df_treat <- read.csv("res/df_treat")



## Exclusive 'control' condition
ids_exc.control <- exclusive_control(df, treat_col1, treat_col2, control_col1, control_col2)

# sanity check 
df %>% 
  filter(Protein.IDs %in% ids_exc.control) %>% View()

df %>% 
  filter(Protein.IDs %in% ids_exc.control) %>% 
  select(3:14) %>% 
  transmute(nsamples_cold = rowSums(across(1:6, ~ as.integer(. > 0))), 
            nsamples_treat = rowSums(across(7:12, ~ as.integer(. > 0)))) %>% 
  View()


# Create the exclusive 'control' condition object and save it to a file
df_control <- df %>% 
  filter(Protein.IDs %in% ids_exc.control)

#write.csv(df_control, file = "res/df_control", row.names = F)
#df_control <- read.csv("res/df_control")


## Both conditions -> multiple hypothesis testing
ids_both <- both_conditions(df, treat_col1, treat_col2, control_col1, control_col2)

proteins <- df %>% 
  filter(Protein.IDs %in% ids_both)


## Discards (Neither exclusive nor for hypothesis testing)
ids_discard <- df %>% 
  filter(!Protein.IDs %in% c(ids_exc.treat, ids_exc.control, ids_both)) %>%
  pull(Protein.IDs)

df %>% 
  filter(Protein.IDs %in% ids_discard) %>% View()

# Create the 'discard'  object and save it to a file
df_discard <- df %>% 
  filter(Protein.IDs %in% ids_discard)

#write.csv(df_discard, file = "res/df_discard", row.names = F)




# Build MSnSet object ---------------------------------------------------------

# Build pdata dataset
pdata <- data.frame(sample = rep(rep(1:3, each = 2),2),
                    rep = rep(rep(1:2, 3), 2),
                    treat = rep(c("cold", "RT"), each = 6))


rownames(pdata) <- colnames(df)[3:14] #select columns with quant (iBAQ) data 


# Build feature dataset
fdata <- data.frame(GenBank = proteins$Protein.IDs)
rownames(fdata) <- proteins$Protein.IDs


# Build expression dataset
edata <- as.matrix(proteins[,3:14]) #select columns with quant (iBAQ) data 
rownames(edata) <- proteins$Protein.IDs
edata[1:10, 1:6]
identical(rownames(pdata), colnames(edata))


# Check the factorial :1 x 2
table(pdata[, c("treat")])

# Build MSnSet object
m <- MSnSet(edata, fdata, pdata)
sampleNames(m)


# # Save MSnSet (and 'proteins' dataset)
# save(m, proteins, file = "data/msnset.RData", compress = TRUE)

# Remove all objects except for functions 
rm(list = setdiff(ls(), lsf.str()))


# Load MSnSet object and proteins dataset 
load(file = "data/msnset.RData")

dim(exprs(m))
dim(fData(m))
dim(pData(m))
# 
edata = exprs(m)
pdata <- pData(m)
fdata <- fData(m)
# 
edata[1:10, 1:6]




# # EDA -------- (raw intensisties) ---------------------------------------------
# 
# check for NA values
sum(is.na(edata))

# N. peptides (not NA) per sample
#m$peptides <-colSums(!is.na(edata))

# check for 0 values
sum(edata == 0)

#N. peptides non-detected over samples
sum(rowSums(edata) == 0)
which(rowSums(edata) == 0)


# N. peptides (not 0) per sample
m$peptides <- colSums(edata != 0)

# Summary over samples
summary(m$peptides)

# Boxplot over samples
boxplot(m$peptides, horizontal = TRUE,
        xlab = "Number of Peptides Detected")

# Boxplot with an ExpressionSet object
# boxplot of the 1st peptide in the expression (abundance) matrix
boxplot(edata[1,] ~ pdata[, "treat"], main = fdata[1,"GenBank"])

# boxplot of the 1000th peptide in the expression (abundance) matrix
boxplot(edata[1000,] ~ pdata[, "treat"], main = fdata[1000,"GenBank"])



# boxplots
boxplot(log2(edata[,1]+1), main = colnames(edata)[1]) #transform: log2+1
boxplot(log2(edata[,7]+1), main = colnames(edata)[7]) #transform: log2+1

boxplot(log2(edata+1), range = 2, cex.axis = 0.65) #plot more than 1 value at once

# density plots
plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=2)
lines(density(log2(edata[,7]+1)),col=3)

#qqplot
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3, pch=19)
abline(c(0,1),lw=2, col=2) #45º line

# PCA
plot_pca(m)
plot_pca(m, phenotype = "treat", legend_title = "treatment") # add a 50% Normal confidence ellipse for each group

#identify the most influential features (peptides) in the PCA space
plot_pca(m, phenotype = "treat",
         biplot = TRUE,
         label_args = list(color = "black"),
         legend_title = "treatment") # add a 50% Normal confidence ellipse for each group

# Heatmap
ehmap <- edata[rowMeans(edata) > 1e6, ] # filter for features with high rowMeans
dim(ehmap)
heatmap(ehmap, cexCol = .7)




# Differential Analysis --------------------------------------------------------
load(file = "data/msnset.RData")

dim(exprs(m))
dim(fData(m))
dim(pData(m))
# 
edata = exprs(m)
pdata <- pData(m)
fdata <- fData(m)
# 
edata[1:10, 1:6]


# Log-transform data **************************************************
#edata <- log(edata+1)
edata <- log2(edata + 1)
edata[1:10, 1:6]



## Linear regression  ----------------------------------------------------

# Contraste : cold vs RT 
col.str = "treat"
levels(as.factor(pData(m)[, col.str])) # baseline es el primer elemento : cold


#coef.str, variable of interest
#model.str, full model, including variable of interest and any covariates

# lm_res <- limma_gen(m, model.str = "~ treat", coef.str = "treat") # valores logFC muy altos!
#lm_res %>% arrange(adj.P.Val) %>% head() #top 6 rows sorted by adjusted P-val

#logFC, slopes of the regression lines 
#AveExpr, average of the values for that feature (contains the y-intercepts)


## Model summary logFC > |1|
#model_summary(lm_res, m, col.str = "treat")


# Features having a significant linear relationship with 'treatment' : top10
# lm_res %>% 
#   filter(adj.P.Val < .05) %>% 
#   arrange(adj.P.Val) %>% 
#   head(10)

# top 3 features sorted by adjusted p-value 
#plot3top(proteins, c("P12787", "O08997", "P20108"))



# Features having a significant linear relationship with 'treatment'sorted by logFC
# lm_res %>% 
#   filter(adj.P.Val < .05, abs(logFC) > 1 ) %>% 
#   arrange(desc(logFC)) %>% head()
# 
# plot3top(proteins, c("P40936", "Q60866", "A2AUR3"))






## t-test : one-comparison ----------------------------------------------------
# differences between the “cold” and “RT” groups
# t_res1 <- limma_a_b(eset = m, model.str = "~ treat", 
#                     coef.str = "treat")
# 
# table(t_res1$adj.P.Val < 0.05) 
# 
# t_res1  %>% arrange(adj.P.Val) %>% head() #top 6 rows sorted by adjusted P-val
# #logFC, difference in means between  "cold" (reference) and "RT" groups
# #AveExpr, overall mean
# 
# levels(as.factor(m$treat))
# #to change the baseline (reference) level and use "RT" as first level (reference): 
# # m$treat <- relevel(as.factor(m$treat), ref = "RT")
# 
# ## Model summary logFC > |1|
# model_summary(t_res1, m, col.str = "treat")
# 
# 
# # Features having a significant difference between the "cold" and "RT" groups: top15
# t_res1 %>% 
#   filter(adj.P.Val < .05, abs(logFC) > 1 ) %>% 
#   arrange(adj.P.Val) %>% 
#   head(15)
# 
# # top 3 features sorted by adjusted p-value
# plot3top(proteins, c("P12787", "O08997", "P20108"))
# 
# 
# # Features having a significant difference between the "cold" and "RT" groups sorted by logFC
# t_res1 %>% 
#   filter(adj.P.Val < .05, abs(logFC) > 1 ) %>% 
#   arrange(desc(logFC)) %>% head()
# 
# plot3top(proteins, c("P40936", "Q60866", "A2AUR3"))


#rowwtests performs Student's t-test (equal variances assumed)

table(pdata[, c("treat")])
g <- rep(c(1,0), each = 6) # define columns per 'group'
ttest.res <- rowttests(edata,factor(g))

#Benjamini & Hochberg correction
fdr <- p.adjust(ttest.res$p.value, method="fdr") 

ttest.res$adj.P.Val <- fdr
head(ttest.res)

# estimate logFC : contraste Cold vs RT !! 

# escoger uno de los 2 casos para el calculo de la media 

# # caso 1: media sobre 6 muestras 
# logFC <- as_tibble(edata) %>% 
#   transmute(mean_Cold = rowMeans(across(1:6), na.rm = TRUE), 
#             mean_RT = rowMeans(across(7:12), na.rm = TRUE),
#             logFC = mean_Cold - mean_RT) %>% 
#   pull(logFC)

# caso2: media sobre valores existentes (0 son tratados como NA)
logFC <- as_tibble(edata) %>%
  rowwise() %>%
  mutate(
    mean_Cold = mean(c_across(1:6)[c_across(1:6) > 0], na.rm = TRUE),
    mean_RT   = mean(c_across(7:12)[c_across(7:12) > 0], na.rm = TRUE),
    logFC = mean_Cold - mean_RT) %>%
  ungroup() %>%
  pull(logFC)


# Combine the FC values to the model 
ttest.res$logFC <- logFC
ttest.res <- ttest.res %>% select(logFC, p.value, adj.P.Val, everything())

head(ttest.res)


## Model summary logFC > |1|
model_summary(ttest.res, m, col.str = "treat")


# Features having a significant difference between the means of the groups: top10
ttest.res %>% 
  filter(adj.P.Val < .05) %>% 
  arrange(adj.P.Val) %>% 
  head(10)

# top 3 features sorted by adjusted p-value 
plot3top(proteins, c("Q9CQ60", "D3Z4X1", "Q8CBG6"))


# Features having a significant difference between the means of the groups: sorted by logFC
ttest.res %>% 
  filter(adj.P.Val < .05) %>% 
  arrange(desc(logFC)) %>% 
  head(10)

plot3top(proteins, c("Q9D881", "P19536", "F7C106"))
proteins %>% filter(Protein.IDs == "Q9D881")
edata['Q9D881',]



# ## Comparison : significant t-test & lm models
# lm_res_sig <- lm_res %>% 
#   filter(adj.P.Val < .05) %>% 
#   arrange(adj.P.Val) 
# 
# 
# dim(t_res1_sig)
# dim(lm_res_sig)
# 
# # Features sig. in t-test not sig. in lm
# sum(!rownames(t_res1_sig) %in% rownames(lm_res_sig))
# # Features sig. in lm not sig. in t-test
# sum(!rownames(lm_res_sig) %in% rownames(t_res1_sig))





## t-test : (Welsch's test) one-comparison -------------------------------------
# Welch t-test does not make the assumption of equal variance 

treat_col1 = "ColdO1.01"
treat_col2 = "ColdO3.02"
control_col1 = "RTO1.01"
control_col2 = "RTO3.02"

dat = edata

treat_cols <- which(colnames(dat) == treat_col1):which(colnames(dat) == treat_col2)
control_cols <- which(colnames(dat) == control_col1):which(colnames(dat) == control_col2)

# Function to perform Welch's t-test for a single row
welch_test <- function(row) {
  group_treat <- row[treat_cols]
  group_control <- row[control_cols]
  test <- t.test(group_treat, group_control, var.equal = FALSE)  # Welch's t-test
  return(c(p.value = test$p.value, t.statistic = test$statistic))
}


# Apply the test to each row of the matrix
wtest.res <- t(apply(edata, 1, welch_test))

# Convert to data frame for easier handling
wtest.res <- as.data.frame(wtest.res)

# Adjust P-values for multiple testing
wtest.res$adj.P.Val <- p.adjust(wtest.res$p.value, method = "fdr")  # Benjamini-Hochberg

# Combine the FC values to the model
wtest.res$logFC <- logFC

# Result 
wtest.res <- wtest.res %>% select(logFC, p.value, adj.P.Val, everything())


## Model summary logFC > |1|
model_summary(wtest.res, m, col.str = "treat")


# Features having a significant difference between the means of the groups: top10
wtest.res %>% 
  filter(adj.P.Val < .05) %>% 
  arrange(adj.P.Val) %>% 
  head(10)

# top 3 features sorted by adjusted p-value 
plot3top(proteins, c("Q921F2", "Q6VYI4", "Q6VYI5"))


# Features having a significant difference between the means of the groups: sorted by logFC
wtest.res %>% 
  filter(adj.P.Val < .05) %>% 
  arrange(desc(logFC)) %>% 
  head(10)

plot3top(proteins, c("Q9D881", "P19536", "F7C106"))
proteins %>% filter(Protein.IDs == "Q9D881")
edata['Q9D881',]



## ANOVA : significant difference between the means of three or more groups. 




## Differential expression analysis with limma --------------------------------

load(file = "data/msnset.RData")

dim(exprs(m))
dim(fData(m))
dim(pData(m))
# 
edata = exprs(m)
pdata <- pData(m)
fdata <- fData(m)
# 
edata[1:10, 1:6]

# log2 transform 
#edata <- log2(edata + 1) #l.309 ya esta log-transformed
edata[1:10, 1:6]


### Pre-procesing for limma ----------------------------------------------------

# Distribution of peptides abundance for all of your samples 
limma::plotDensities(edata, legend = F)
limma::plotDensities(edata, group = pdata[, "treat"], legend = "topright")

# almost all densities lie ~ 0 (not here in this example)
# only a subset of the detected peptides are relevant to the case study 

# because most of the data lies near zero with a very long right tail, the data set 
# likely contains measurements for many peptides that are not relevant for the study 
# and should be removed. 


# step1 : log transform #######OJO! ya estaban log-transformados
# to view the entire distribution
#edata <- log(edata+1)
limma::plotDensities(edata, legend = F)
# distributions are not the same across the samples 


# step2 : Quantile normalize
# to transform each sample to the same empirical distribution
edata <- limma::normalizeBetweenArrays(edata)
limma::plotDensities(edata, legend = F)

# step3 : filtering and remove irrelevant genes/proteins
#choose a cutoff that excludes the peak of genes with low expression levels 
#abline(v = 5)
#keep <- rowMeans(edata) > 5
#edata <- edata[keep, ]
#limma::plotDensities(edata, legend = F)


# (lastly) Accounting for technical batch effects 
# Identify the largest sources of variation in the dataset
limma::plotMDS(edata, labels = pData(m)[,"treat"], gene.selection = "common", main ="Treatment")
limma::plotMDS(edata, labels = pData(m)[,"rep"], gene.selection = "common", main = "Repetition")
limma::plotMDS(edata, labels = pData(m)[,"sample"], gene.selection = "common", main = "Sample")

#The major sources of variation correlate with the variables of interest, so no need to correct for technical batch effects. 



## now perform the differential expression test  !

# Build MSnSet object 
m <- MSnSet(edata, fdata, pdata)
sampleNames(m)

# Save MSnSet (and 'proteins' dataset)
#save(m, df, proteins, file = "data/msnset_processed.RData", compress = TRUE)

# Remove all objects except for functions 
#rm(list = setdiff(ls(), lsf.str()))


# Load dataset processed for limma analysis
load("data/msnset_processed.RData")

fdata <- fData(m)
edata <- exprs(m)
pdata <- pData(m)


#step1 : build the design matrix
group <- as.factor(pdata$treat)
#to change the baseline (reference) level and use "RT" as first level (reference):  
group <- relevel(group, ref = "RT")

design = model.matrix(~0 + group) # 0 means no intercept for the linear model
colnames(design) = gsub("group","", colnames(design))
# count the number of samples modeled by each coefficient
colSums(design)


#step2 : construct the contrast matrix with `makeContrast`
x <- c("cold-RT")
contrast <- limma::makeContrasts(contrasts = x, levels = design)
contrast  # view the contrast matrix : cold vs RT (same as ttest or wtest)

#step3 : fit the model coefficients
fit1 <- limma::lmFit(na.omit(edata), design) # pensar si log o log2

# edata esta pre-procesado ! log transf - step1 pre-processing


# fit the contrasts
fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)

# calculate the t-statistics
limma_res <- limma::eBayes(fit2)

# summarize results 
summary(limma::decideTests(limma_res))
# venn diagram
limma::vennDiagram(limma_res)

limma::volcanoplot(fit = limma_res, coef = x, highlight = 5, names = rownames(limma_res$p.value))


# extract limma results using `topTable` function
# coef = 1 allows you to extract the specific contrast
# n=Inf output all rows

limma_res$contrasts
limma.res <- limma::topTable(limma_res, coef = 1, n= Inf)
head(limma.res)
dim(limma.res)

table(limma.res$ adj.P.Val  < 0.05)





## DEqMS analysis --------------------------------------------------------------

load(file = "data/msnset.RData")

dim(exprs(m))
dim(fData(m))
dim(pData(m))
# 
edata = exprs(m)
pdata <- pData(m)
fdata <- fData(m)
# 


# intensidades: DEqMS requiere que los datos estén en escala log2, porque está construido
# sobre el framework de limma, que asume log2 intensities como input.
# samples require to have medians centered. If not, do median centering. 

# log2 transform 
edata[1:10, 1:6]
edata <- log2(edata + 1)
edata[1:10, 1:6]


boxplot(na.omit(edata), las=2, main = "dataset") # non-median centered
dat.log = DEqMS::equalMedianNormalization(edata)
boxplot(na.omit(dat.log), las=2, main = "dataset") # median centered


# Peptides count : This method models the variance as a function of the number of peptides identified
# per protein, and for that it needs to work with linear intensities,
# not log-transformed.


# Load Data 
psm.count  <- read.csv2("data/T1_JDie_09042025.csv", header=TRUE)

# Extract columns with PSM data for the analysis
target <- c( "Protein.IDs", "Razor...unique.peptides.Cold.Liver_LDLRKO1.01", 
             "Razor...unique.peptides.Cold.Liver_LDLRKO1.02", 
             "Razor...unique.peptides.Cold.Liver_LDLRKO2.01", 
             "Razor...unique.peptides.Cold.Liver_LDLRKO2.02", 
             "Razor...unique.peptides.Cold.Liver_LDLRKO3.01", 
             "Razor...unique.peptides.Cold.Liver_LDLRKO3.02", 
             "Razor...unique.peptides.RT.Liver_LDLRKO1.01", 
             "Razor...unique.peptides.RT.Liver_LDLRKO1.02", 
             "Razor...unique.peptides.RT.Liver_LDLRKO2.01", 
             "Razor...unique.peptides.RT.Liver_LDLRKO2.02", 
             "Razor...unique.peptides.RT.Liver_LDLRKO3.01", 
             "Razor...unique.peptides.RT.Liver_LDLRKO3.02")


psm.count <- psm.count  %>% 
  select(all_of(target))

psm.count %>% tibble()

# Expand rows 
df_long <- rows_expansion(psm.count)

df_long %>% dim()
df_long %>% head()


# The object 'proteinas' contains the filtered data: the column Protein.IDS represent 
# the peptides that are not exclusive to any condition and free of contaminants.

#load("data/msnset_processed.RData")
#proteins

# variance estimate based on minimum number of psms per protein used for quantification
df_long <- df_long %>% 
  filter(Protein.IDs %in% proteins$Protein.IDs) %>% 
  rowwise() %>%
  mutate(min_pept = min(c_across(starts_with("Razor")), na.rm = TRUE)) %>%
  ungroup() %>% 
  select(Protein.IDs, min_pept)

# read the input protein table
psm.count.table = data.frame(count = df_long$min_pept, 
                             row.names =  df_long$Protein.IDs)

psm.count.table %>% head()

summary(psm.count.table$count) ##el minimo PSM algunas es 0
rownames(psm.count.table)[which(psm.count.table$count == 0)]
min.psm <- rownames(psm.count.table)[which(psm.count.table$count == 0)]


# STEP1 : Make design table : how samples are arranged in different goups/classes
group <- as.factor(pdata$treat)
#to change the baseline (reference) level and use "RT" as first level (reference):  
group <- relevel(group, ref = "RT")

# Build the design matrix
design = model.matrix(~0 + group) # 0 means no intercept for the linear model
colnames(design) = gsub("group","", colnames(design))
# count the number of samples modeled by each coefficient
colSums(design)



#STEP2 : Make contrasts 
x <- c("cold-RT")
contrast <- limma::makeContrasts(contrasts = x, levels = design)
contrast  # view the contrast matrix : cold vs RT (same as ttest or wtest)

#STEP3 : fit the model coefficients
# only for proteins with min psm > 0 :
all_proteins <- rownames(dat.log)
target_proteins <- setdiff(all_proteins, min.psm) # elements in 'all_proteins' but not in 'min.psm'
fit1 <- limma::lmFit(na.omit(dat.log[target_proteins, ]), design)
# fit the contrasts
fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)
# calculate the t-statistics
fit3 <- limma::eBayes(fit2)

# SEPT 4 : DEqMS analysis 
# assign an extra variable `count` to `limma` object, telling how many PSMs are 
# quantifed for each protein
fit3$count <- psm.count.table[rownames(na.omit(dat.log[target_proteins, ])), "count"]
deqMS_res = DEqMS::spectraCounteBayes(fit3)


# visualize the fit curve : variance dependence on quantified PSM
# show only proteins quantified by <= 30 PSMs
DEqMS::VarianceBoxplot(deqMS_res, 30, main="Alfonso dataset", xlab="PSM count")

# DEqMS models the relationship between variance and the number of PSMs per protein. 
# Although no fixed threshold is specified in the vignette, the dependence on PSM 
# counts is clearly emphasized, and proteins with few PSMs are treated with greater
# statistical rigor. This approach implies and implicit filter : having a reasonable 
# minimum number of PSMs is a requirement to obtain reliable variance estimates
# (~ min > 10 in this case)

DEqMS::VarianceScatterplot(deqMS_res)
# varianza estimada de la expresion proteica en funcion del n. de PSMs por proteina (log scale)
# proteinas con PSM bajos tienen varianza mas alta 
# proteinas con altos PSMs tienen varrianzas mas bajas y estables (= mas confianza en su cuantificacion)
# figura ayuda a decidir si vale la pena filtrar proteinas con muy pocos PSMs al observar que su 
# varianza es excesivamente alta 


# extract the results as a data frame and save it 
head(deqMS_res$coefficients) # each col is a specific contrast
DEqMS.res = DEqMS::outputResult(deqMS_res,coef_col = 1)
head(DEqMS.res)


# Save results (models) --------------------------------------------------------
save(ttest.res, wtest.res, limma.res, DEqMS.res, file = "res/models.RData")


# Volcano plot
# Usar logFC del modelo limma; sigue siendo válido porque DEqMS solo ajusta la varianza, no el efecto.
ggplot(DEqMS.res, aes(x = logFC, y = -log10(sca.adj.pval))) +
  geom_point(aes(color = sca.adj.pval < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  labs(
    title = "Volcano Plot - DEqMS",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Significant"
  ) +
  theme_minimal()


# Volcano plot con clasificación por significancia + FC
ggplot(DEqMS.res, aes(x = logFC, y = -log10(sca.adj.pval))) +
  geom_point(aes(
    color = (sca.adj.pval < 0.05 & abs(logFC) > 1)
  ), alpha = 0.7) +
  scale_color_manual(
    values = c(`TRUE` = "red", `FALSE` = "grey90"),
    guide = guide_legend(title = NULL),  # Quitar título de la leyenda
    labels = c("Non-significant", "Significant")  # Etiquetas en leyenda
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot - DEqMS",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)"
  ) +
  theme_minimal()



## MSstats analysis --------------------------------------------------------------

# MSstats (when used with MaxQuant data) requires three main files
 
# proteinGroups.txt – Used to associate peptides with proteins.
#                   – Output file from MaxQuant.
 
# evidence.txt – Contains raw peptide-level information.
#               – Output file from MaxQuant.
 
# annotation.csv – Must be created manually.

#annotation_msstats
annotation <- data.frame(
  Raw.file = c("Cold-Liver_LDLRKO1-01", "Cold-Liver_LDLRKO1-02",
               "Cold-Liver_LDLRKO2-01", "Cold-Liver_LDLRKO2-02",
               "Cold-Liver_LDLRKO3-01", "Cold-Liver_LDLRKO3-02",
               "RT-Liver_LDLRKO1-01",   "RT-Liver_LDLRKO1-02",
               "RT-Liver_LDLRKO2-01",   "RT-Liver_LDLRKO2-02",
               "RT-Liver_LDLRKO3-01",   "RT-Liver_LDLRKO3-02"),
  Condition = rep(c("Cold", "RT"), each = 6),
  BioReplicate = rep(1:3, each = 2, times = 2), 
  TechReplicate = rep(1:2, times = 6))

#evidence_df
evidence_df <- read.delim("data/evidence.txt", sep = "\t", header = TRUE)


#proteinGroups : archivo sin modificar que exporta MaxQuant
proteinGroups  <- read.csv2("data/T1_JDie_09042025.csv", header=TRUE)

# at least, you'll need this columns : 
target <- c("Protein.IDs", "Protein.names", "Gene.names", "Score", "Intensity", "MS.MS.count",
            "Reverse", "Potential.contaminant", "id", "MS.MS.IDs", "Best.MS.MS",
            "Oxidation..M..site.IDs","Taxonomy.IDs" )

proteinGroups <- proteinGroups %>% 
  select(all_of(target))

# Convert MaxQuant data to MSstats format
converted <- MaxQtoMSstatsFormat(
  proteinGroups = proteinGroups,
  annotation = annotation,
  evidence = evidence_df)

# MSstats performs log2 transformation automatically
# When you use functions like MaxQtoMSstatsFormat() or dataProcess(), MSstats applies
# log2 transformation to intensity values if they are not already transformed.

# Process the data with MSstats 
msstats_processed_data <- dataProcess(converted, censoredInt = '0')

# uso '0' en lugar de 'NA' porque los valores faltantes son ceros explicitos 
#y evito el warning 
#Warning message:
#  In survreg.fit(X, Y, weights, offset, init = init, controlvals = control,  :
#                   Ran out of iterations and did not converge


# Define the contrast matrix (Cold vs RT)
comparison_msstats_cold_vs_rt <- matrix(c(1, -1), 
                                        nrow = 1, 
                                        byrow = TRUE, 
                                        dimnames = list("Cold_vs_RT", c("Cold", "RT")))

# Perform the group comparison
msstats_results_data <- groupComparison(contrast.matrix = comparison_msstats_cold_vs_rt, 
                                        data = msstats_processed_data)


# Extract the complete results 
all_msstats_results <- msstats_results_data$ComparisonResult

# expand rows with multiple ids 
all_msstats_results[1:3, 1:4]
msstats.res <- rows_expansion(all_msstats_results)
msstats.res <- msstats.res %>% filter(log2FC != Inf)
# lo filtro para los peptidos con PSM min 10 ?? v. dat.log ----------------------



# Save results (models) --------------------------------------------------------
save(ttest.res, wtest.res, limma.res, DEqMS.res, msstats.res, file = "res/models.RData")




## -----------------------------------------------------------------------------





# load models 
load("res/models.RData")


write.csv(ttest_res, file = "~/Desktop/ttest.csv")
write.csv(ttest_res, file = "~/Desktop/wtest.csv")
write.csv(limma.res, file = "~/Desktop/limma.csv")
write.csv(DEqMS.res, file = "~/Desktop/deqMS.csv")
write.csv(DEqMS.res, file = "~/Desktop/msstats.csv")





# Results visualzation --------------------------------------------------------
# Distribution of P-values from the collection of hypothesis tests : linear regression
# hist(lm_res$P.Value, breaks = seq(0,1,0.05),
#      main = "Histogram of P-values from t-test Results", 
#      xlab = "P-value")
# 


# Distribution of P-values from the collection of hypothesis tests : t test 
hist(ttest.res$adj.P.Val, breaks = seq(0,1,0.05),
     main = "Histogram of P-values (adj.) from t-test Results", 
     xlab = "P-value")
#There is a peak around 0 that indicates the null hypothesis is false for some of the tests. 


# Normally, the adjusted p-values would be used, thought if they are high in the results, the unadjusted p-values are commonly used.
plot_volcano(df = ttest.res, logFC = "logFC", 
             pvals = "adj.P.Val", sig_threshold = 0.05)

# plot_volcano(df = t_res1, logFC = "logFC", 
#              pvals = "P.Value", sig_threshold = 0.05)


# label top features
ttest_res$GenBank <- rownames(ttest_res)
plot_volcano(df = ttest_res, logFC = "logFC", 
             pvals = "adj.P.Val", sig_threshold = 0.05, 
             label = "GenBank", 
             num_features = 5)



# Compare p-values from DEqMS to other tests ----------------------------------
head(ttest.res)
dim(ttest.res)

head(wtest.res)
dim(wtest.res)

head(limma.res)
dim(limma.res)

head(DEqMS.res)
dim(DEqMS.res)

head(msstats.res)
dim(msstats.res)


# Number of significant proteins 
table(ttest.res$adj.P.Val < 0.05)
table(wtest.res$adj.P.Val < 0.05)
table(limma.res$ adj.P.Val  < 0.05)
table(DEqMS.res$sca.P.Value < 0.05)
table(msstats.res$adj.pvalue < 0.05)


# visualize the distribution of p-values by different analysis
plot(sort(-log10(limma.res$adj.P.Val),decreasing = TRUE), 
     type="l",lty=2,lwd=2, ylab="-log10(p-value)",ylim = c(0,10),
     xlab="Proteins ranked by p-values",
     col="purple")

lines(sort(-log10(ttest.res$adj.P.Val),decreasing = TRUE), 
      lty=2,lwd=2,col="orange")

lines(sort(-log10(wtest.res$adj.P.Val),decreasing = TRUE), 
      lty=2,lwd=2,col="green")

lines(sort(-log10(DEqMS.res$sca.adj.pval),decreasing = TRUE), 
      lty=1,lwd=2,col="red")

lines(sort(-log10(msstats.res$adj.pvalue),decreasing = TRUE), 
      lty=2,lwd=2,col="steelblue")



legend("topright",legend = c("Limma","DEqMS","t.test", "w.test", "MSstats"),
       col = c("purple","red","orange", "green", "steelblue"),lty=c(2,1,2),lwd=2)


abline(h = - log10(0.05))


# visualize the distribution of adjusted p-values by different analysis
# el tutorial en bioconductor usa P vals pero lo hago con adj P vals y sale parecido

# plot(sort(-log10(limma.res$adj.P.Val),decreasing = TRUE), 
#      type="l",lty=2,lwd=2, ylab="-log10(p-value)",ylim = c(0,10),
#      xlab="Proteins ranked by adjusted p-values",
#      col="purple")
# 
# lines(sort(-log10(DEqMS.res$sca.adj.pval),decreasing = TRUE), 
#       lty=1,lwd=2,col="red")
# 
# lines(sort(-log10(t_res1$adj.P.Val),decreasing = TRUE), 
#       lty=2,lwd=2,col="orange")
# 
# 
# legend("topright",legend = c("Limma","DEqMS","t.test"),
#        col = c("purple","red","orange"),lty=c(2,1,2),lwd=2)
# 
# 
# abline(h = - log10(0.05))





# Plotting top500 proteins ranked by P-values
plot(sort(-log10(limma.res$adj.P.Val),decreasing = TRUE)[1:500], 
     type="l",lty=2,lwd=2, ylab="-log10(p-value)",ylim = c(0,10),
     xlab="Proteins ranked by p-values",
     col="purple")

lines(sort(-log10(ttest.res$adj.P.Val),decreasing = TRUE)[1:500], 
      lty=2,lwd=2,col="orange")

lines(sort(-log10(wtest.res$adj.P.Val),decreasing = TRUE)[1:500], 
      lty=2,lwd=2,col="green")

lines(sort(-log10(DEqMS.res$sca.adj.pval),decreasing = TRUE)[1:500], 
      lty=1,lwd=2,col="red")

lines(sort(-log10(msstats.res$adj.pvalue),decreasing = TRUE)[1:500], 
      lty=2,lwd=2,col="steelblue")





# Combine all model results into a single table --------------------------------------
# 1. Get ttest results
ttest_tbl <- ttest.res %>% 
  tibble::rownames_to_column("Protein.ID") %>% 
  select(Protein.ID, logFC_ttest = logFC, adj.P.ttest = adj.P.Val)

# 2. Get wtest results 
wtest_tbl <- wtest.res %>% 
  tibble::rownames_to_column("Protein.ID") %>% 
  select(Protein.ID, logFC_wtest = logFC, adj.P.wtest = adj.P.Val)


# 3. Get limma topTable
limma_tbl <- limma.res %>% 
  tibble::rownames_to_column("Protein.ID") %>%
  select(Protein.ID, logFC_limma = logFC, adj.P.limma = adj.P.Val)

# 4. Get DEqMS resuls
deqms_tbl <- DEqMS.res %>% 
  tibble::rownames_to_column("Protein.ID") %>% 
  select(Protein.ID, logFC_deqms = logFC, sca.adj.pval) 


# 5. Get MSstats 
msstas_tbl <- msstats.res %>% 
  select(Protein, logFC_msstats = log2FC, adj.P.ms = adj.pvalue) %>% 
  dplyr::rename(Protein.ID = Protein)


# 6. Combine tables 
pval_cutoff <- 0.05

models_merged <- full_join(ttest_tbl, wtest_tbl, by = "Protein.ID") %>% 
  full_join(deqms_tbl, by = "Protein.ID") %>% 
  full_join(limma_tbl, by = "Protein.ID") %>% 
  full_join(msstas_tbl, by = "Protein.ID") %>% 
  mutate(sig_ttest = adj.P.ttest < pval_cutoff, 
         sig_wtest = adj.P.wtest < pval_cutoff, 
         sig_limma = adj.P.limma < pval_cutoff, 
         sig_deqms = sca.adj.pval < pval_cutoff, 
         sigmsstats = adj.P.ms < pval_cutoff)

models_merged %>% head()


models_merged %>% 
  summarize(ttest = sum(sig_ttest, na.rm = T) ,
            wtest = sum(sig_wtest, na.rm = T), 
            limma = sum(sig_limma, na.rm = T), 
            deqms = sum(sig_deqms, na.rm = T), 
            msstats = sum(sigmsstats, na.rm = T)) %>% 
  pivot_longer(cols = 1:5, names_to = "test", values_to = "n.significant") %>% 
  ggplot(aes(x = test, y = n.significant)) + 
  geom_col() + 
  labs(title = "Number of Significant Proteins (adj.p < 0.05)", x = "Test", y = "Number of Proteins") +
  theme_minimal()

Comparaciones visuales de ADOLFO : 
  # --- 2) Diagrama de Venn con VennDiagram (Intento con ajustes de tamaño) ---
  ###### Distribución p values (todos los ajustados) (Adaptado para archivos de 2_Htesting) ########
# --- 4) Gráficos de líneas de -log10(p.value ajustado) rankeado ---




install.packages("VennDiagram")
library(VennDiagram)

models_merged %>% head()
sig_ttest <- models_merged %>% filter(sig_ttest) %>% pull(Protein.ID)
sig_limma <- models_merged %>% filter(sig_limma) %>% pull(Protein.ID)
sig_deqms <- models_merged %>% filter(sig_deqms) %>% pull(Protein.ID)

venn.plot <- venn.diagram(
  x = list(
    limma = sig_limma ,
    DEqMS = sig_deqms,
    ttestest = sig_ttest
  ),
  filename = NULL,
  col = "black",
  fill = c("lightblue", "salmon", "palegreen"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = c(-20, 20, 0),
  cat.dist = 0.05,
  main = "Venn Diagram of Significant Proteins by Method"
)

grid::grid.newpage()
grid::grid.draw(venn.plot)




# Comparative limma vs DEqMS
gg <- ggplot(models_merged, aes(x = -log10(adj.P.limma), y = -log10(sca.adj.pval))) +
  geom_point(aes(color = sig_limma | sig_deqms), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "DEqMS vs limma",
    x = "-log10(adj.P.Val) (limma)",
    y = "-log10(sca.adj.pval) (DEqMS)",
    color = paste0("Significant (p <", pval_cutoff, ")")
  ) +
  scale_color_manual(values = c("grey70", "red")) +
  theme_minimal()

print(gg)




