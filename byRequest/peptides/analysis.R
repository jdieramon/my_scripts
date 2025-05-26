# Dependencies -----------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(MSnSet.utils)


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


# Make an index for rows with >1 Protein.IDs
ids <- prots %>% pull(Protein.IDs )
idx <- str_which(string = ids, pattern = ";")

# Subset the dataframe that requires tidying
df <- prots %>% slice(idx)

# Make a tidy dataframe 
df = lapply(seq(nrow(df)), function(i) make_longer(df[i,]))
df <- bind_rows(df)

# Combine with rows with = 1Protein.IDs
df <- df %>% 
  bind_rows(prots %>% slice(-idx))


# Remove contaminats CON_ & 'señuelo' REV_
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

df %>% 
  filter(Protein.IDs %in% ids_exc.treat) %>% View()


# sanity check 
df %>% 
  filter(Protein.IDs %in% ids_exc.treat) %>% 
  select(3:14) %>% 
  transmute(nsamples_cold = rowSums(across(1:6, ~ as.integer(. > 0))), 
            nsamples_treat = rowSums(across(7:12, ~ as.integer(. > 0)))) %>% 
  View()


## Exclusive 'control' condition
ids_exc.control <- exclusive_control(df, treat_col1, treat_col2, control_col1, control_col2)

df %>% 
  filter(Protein.IDs %in% ids_exc.control) %>% View()

# sanity check 
df %>% 
  filter(Protein.IDs %in% ids_exc.control) %>% 
  select(3:14) %>% 
  transmute(nsamples_cold = rowSums(across(1:6, ~ as.integer(. > 0))), 
            nsamples_treat = rowSums(across(7:12, ~ as.integer(. > 0)))) %>% 
  View()



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

# Make some tidy 
#rm(prots, ids, idx, names_tidy, target, df,edata,fdata,pdata)


load(file = "data/msnset.RData")

dim(exprs(m))
dim(fData(m))
dim(pData(m))
# 
edata = exprs(m)
edata[1:10, 1:6]
# 
pdata <- pData(m)
fdata <- fData(m)




# # EDA --------------------------------------------------------------------------
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











# # Pre-processing --------------------------------------------------------------
# # Zeros, crean problemas : imputo por el valor medio del grupo 
# 
# # remove protein id. 
# tmp <- proteins[, -1]
# 
# # impute 
# proteins <- impute_mean(tmp[,1:6]) %>% 
#   bind_cols(impute_mean(tmp[,7:12])) %>% 
#   mutate(Protein.IDs = proteins$Protein.IDs) %>% 
#   select(Protein.IDs, everything()) %>% 
#   tibble()


# # Re-Build MSnSet object :
# 
# # Re-Build expression dataset for DA
# edata <- proteins %>% 
#   select(-Protein.IDs) %>% 
#   as.matrix()
# 
# rownames(edata) <- proteins$Protein.IDs
# 
# 
# # Re-Build feature dataset for DA
# fdata <- data.frame(GenBank = proteins$Protein.IDs)
# rownames(fdata) <- proteins$Protein.IDs
# 
# # Check the factorial :1 x 2
# table(pdata[, c("treat")])
# 
# # Build MSnSet object 
# m <- MSnSet(edata, fdata, pdata)
# sampleNames(m)
# 
# # Save MSnSet (and 'proteins' dataset)
# save(m, proteins, file = "data/msnset.RData", compress = TRUE)
# load(file = "data/msnset.RData")
# 
# 
# # Make some tidy
# rm(df_pr, tmp)






# Differential Analysis ---------------------------------------------
load(file = "data/msnset.RData")

# dim(exprs(m))
# dim(fData(m))
# dim(pData(m))
# 
# edata = exprs(m)
# edata[1:10, 1:6]
edata <- log(edata+1)
edata[1:10, 1:6]

# 
# pdata <- pData(m)



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


#BiocManager::install("genefilter")
library(genefilter)


g <- rep(c(1,0), each = 6) # define columns per 'group'
ttest_res <- rowttests(edata,factor(g))

#Benjamini & Hochberg correction
fdr <- p.adjust(ttest_res$p.value, method="fdr") 

ttest_res$adj.P.Val <- fdr
head(ttest_res)

# estimate logFC : RT / Cold !!
logFC <- proteins %>% 
  transmute(mean_RT = rowMeans(across(3:8), na.rm = TRUE), 
            mean_Cold = rowMeans(across(9:14), na.rm = TRUE),
            logFC = mean_RT / mean_Cold) %>% 
  pull(logFC)


ttest_res$logFC <- logFC
head(ttest_res)


## Model summary logFC > |1|
#model_summary(lm_res, m, col.str = "treat")
model_summary(ttest_res, m, col.str = "treat")

# Features having a significant linear relationship with 'treatment' : top10
ttest_res %>% 
  filter(adj.P.Val < .05) %>% 
  arrange(adj.P.Val) %>% 
  head(10)

# top 3 features sorted by adjusted p-value 
#plot3top(proteins, c("P12787", "O08997", "P20108"))
plot3top(proteins, c("Q9CQ60", "D3Z4X1", "Q8CBG6"))

proteins %>% filter(Protein.IDs == 'Q9CQ60') %>% 
  transmute(mean_RT = rowMeans(across(3:8), na.rm = TRUE), 
            mean_Cold = rowMeans(across(9:14), na.rm = TRUE),
            logFC = mean_RT / mean_Cold) 



### COLD deberia ser mas pequeno. 
## como estoy haciendo el contrastre 
### aclarar este punto !!!!! 













# Features having a significant linear relationship with 'treatment'sorted by logFC
lm_res %>% 
  filter(adj.P.Val < .05, abs(logFC) > 1 ) %>% 
  arrange(desc(logFC)) %>% head()

plot3top(proteins, c("P40936", "Q60866", "A2AUR3"))




## t-test : one-comparison ----------------------------------------------------
# differences between the “cold” and “RT” groups
t_res1 <- limma_a_b(eset = m, model.str = "~ treat", 
                    coef.str = "treat")

table(t_res1$adj.P.Val < 0.05) 

t_res1  %>% arrange(adj.P.Val) %>% head() #top 6 rows sorted by adjusted P-val
#logFC, difference in means between  "cold" (reference) and "RT" groups
#AveExpr, overall mean

levels(as.factor(m$treat))
#to change the baseline (reference) level and use "RT" as first level (reference): 
# m$treat <- relevel(as.factor(m$treat), ref = "RT")

## Model summary logFC > |1|
model_summary(t_res1, m, col.str = "treat")


# Features having a significant difference between the "cold" and "RT" groups: top15
t_res1 %>% 
  filter(adj.P.Val < .05, abs(logFC) > 1 ) %>% 
  arrange(adj.P.Val) %>% 
  head(15)

# top 3 features sorted by adjusted p-value
plot3top(proteins, c("P12787", "O08997", "P20108"))


# Features having a significant difference between the "cold" and "RT" groups sorted by logFC
t_res1 %>% 
  filter(adj.P.Val < .05, abs(logFC) > 1 ) %>% 
  arrange(desc(logFC)) %>% head()

plot3top(proteins, c("P40936", "Q60866", "A2AUR3"))







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




## ANOVA : significant difference between the means of three or more groups. 


## Differential expression analysis with limma --------------------------------


### Pre-procesing for limma -------------------------------------------------------------

# Distribution of peptides abundance for all of your samples 
fdata <- fData(m)
edata <- exprs(m)
pdata <- pData(m)
limma::plotDensities(edata, legend = F)
limma::plotDensities(edata, group = pdata[, "treat"], legend = "topright")

# almost all densities lie ~ 0
# only a subset of the detected peptides are relevant to the case study 

# because most of the data lies near zero with a very long right tail, the data set 
# likely contains measurements for many peptides that are not relevant for the study 
# and should be removed. 


# step1 : log transform #######OJO! ya estaban log-transformados
# to view the entire distribution
edata <- log(edata+1)
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
limma::plotMDS(edata, labels = pData(m)[,"treat"], gene.selection = "common")
limma::plotMDS(edata, labels = pData(m)[,"rep"], gene.selection = "common")
limma::plotMDS(edata, labels = pData(m)[,"sample"], gene.selection = "common")

#The major sources of variation correlate with the variables of interest, so no need to correct for technical batch effects. 




## now perform the differential expression test  !

# Build MSnSet object 
m <- MSnSet(edata, fdata, pdata)
sampleNames(m)

# Save MSnSet (and 'proteins' dataset)
#save(m, df, proteins, file = "data/msnset_processed.RData", compress = TRUE)

# Remove all objects except for functions 
#rm(list = setdiff(ls(), lsf.str()))









#step1 : build the design matrix
group <- as.factor(pdata$treat)
#to change the baseline (reference) level and use "RT" as first level (reference):  
#group <- relevel(group, ref = "RT")

design = model.matrix(~0 + group) # 0 means no intercept for the linear model
colnames(design) = gsub("group","", colnames(design))
# count the number of samples modeled by each coefficient
colSums(design)


#step2 : construct the contrast matrix with `makeContrast`
x <- c("RT-cold")
contrast <- limma::makeContrasts(contrasts = x, levels = design)
contrast  # view the contrast matrix

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



## DEqMS analysis --------------------------------------------------------------

# Este método modela la varianza en función del número de péptidos identificados 
# por proteína, y para eso necesita trabajar con intensidades lineales, 
# no log-transformadas.


# read the input protein table 
proteins[1:3,] 
summary(proteins$MS.MS.count) ##el minimo PSM algunas es 0 !! ----SEGUIR @@@@@@@@@@@@@@
proteins[which(proteins$MS.MS.count == 0),]


psm.count.table = data.frame(count = proteins$MS.MS.count, 
                             row.names =  proteins$Protein.IDs)

# assign a extra variable `count` to `limma_res` object, telling how many PSMs are 
# quantifed for each protein
limma_res$count = psm.count.table[rownames(limma_res$coefficients),"count"]
summary(limma_res$count)

deqMS_res = DEqMS::spectraCounteBayes(limma_res)

# alfonso notas : 
# topTable(deqMS_res, adujunst.nethod = "fdr", n = Inf) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# # para calcular el minimum number of PSMs used for quantification no puede haber 0s
# # Zeros, crean problemas : imputo por el valor medio del grupo 
# 
# # remove protein id. 
# tmp <- proteins[, -c(1:14)]
# 
# # impute 
# proteins_pr <- impute_mean(tmp[,1:6]) %>% 
#   bind_cols(impute_mean(tmp[,7:12])) %>% 
#   mutate(Protein.IDs = proteins$Protein.IDs) %>% 
#   select(Protein.IDs, everything()) %>% 
#   tibble()
# 
# proteins_pr[1:3, ]
# 
# count_columns = 2:13
# psm.count.table = data.frame(count = matrixStats::rowMins(
#   as.matrix(proteins_pr[,count_columns])), row.names =  proteins$Protein.IDs)




# visualize the fit curve : variance dependence on quantified PSM
# show only proteins quantified by <= 40 PSMs
DEqMS::VarianceBoxplot(deqMS_res, n=40, main="Alfonso dataset",xlab="PSM count")

DEqMS::VarianceScatterplot(deqMS_res, main = "Alfonso dataset")



# exctract the results as a data frame and save it 
head(deqMS_res$coefficients) # each col is a specific contrast
DEqMS.res = DEqMS::outputResult(deqMS_res,coef_col = 1)
head(DEqMS.res)






# Results visualzation --------------------------------------------------------
# Distribution of P-values from the collection of hypothesis tests : linear regression
hist(lm_res$P.Value, breaks = seq(0,1,0.05),
     main = "Histogram of P-values from t-test Results", 
     xlab = "P-value")

# Distribution of P-values from the collection of hypothesis tests : t test 
hist(t_res1$P.Value, breaks = seq(0,1,0.05),
     main = "Histogram of P-values from t-test Results", 
     xlab = "P-value")
#There is a peak around 0 that indicates the null hypothesis is false for some of the tests. 


# Normally, the adjusted p-values would be used, thought if they are high in the results, the unadjusted p-values are commonly used.
plot_volcano(df = lm_res, logFC = "logFC", 
             pvals = "adj.P.Val", sig_threshold = 0.05)

# plot_volcano(df = t_res1, logFC = "logFC", 
#              pvals = "P.Value", sig_threshold = 0.05)


# label top features
t_res1$GenBank <- rownames(t_res1)
plot_volcano(df = t_res1, logFC = "logFC", 
             pvals = "adj.P.Val", sig_threshold = 0.05, 
             label = "GenBank", 
             num_features = 5)



# Compare p-values from DEqMS to other tests ----------------------------------
head(t_res1)
dim(t_res1)
table(t_res1$adj.P.Val < 0.05)

# extract limma results using `topTable` function
# coef = 1 allows you to extract the specific contrast
# n=Inf output all rows 

limma_res$contrasts
limma.res <- limma::topTable(limma_res, coef = 1, n= Inf)
head(limma.res)
dim(limma.res)

table(limma.res$ adj.P.Val  < 0.05)


head(DEqMS.res)
table(DEqMS.res$adj.P.Val < 0.05)


# Save results (models) ---------------------------
#save(t_res1, limma.res, DEqMS.res, file = "res/models.RData")
load("res/models.RData")


# Number of significant proteins 
table(t_res1$adj.P.Val < 0.05)
table(limma.res$ adj.P.Val  < 0.05)
table(DEqMS.res$adj.P.Val < 0.05)



# visualize the distribution of p-values by different analysis
plot(sort(-log10(limma.res$P.Value),decreasing = TRUE), 
     type="l",lty=2,lwd=2, ylab="-log10(p-value)",ylim = c(0,15),
     xlab="Proteins ranked by p-values",
     col="purple")

lines(sort(-log10(DEqMS.res$sca.P.Value),decreasing = TRUE), 
      lty=1,lwd=2,col="red")

lines(sort(-log10(t_res1$P.Value),decreasing = TRUE), 
      lty=2,lwd=2,col="orange")


legend("topright",legend = c("Limma","DEqMS","t.test"),
       col = c("purple","red","orange"),lty=c(2,1,2),lwd=2)


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
     type="l",lty=2,lwd=2, ylab="-log10(p-value)",ylim = c(0,12),
     xlab="Peptides/proteins ranked by p-values",
     col="purple")

lines(sort(-log10(DEqMS.res$sca.adj.pval),decreasing = TRUE)[1:500], 
      lty=1,lwd=2,col="red")

lines(sort(-log10(t_res1$adj.P.Val),decreasing = TRUE)[1:500], 
      lty=2,lwd=2,col="orange")

legend("topright",legend = c("Limma","DEqMS","t.test"),
       col = c("purple","red","orange"),lty=c(2,1,2),lwd=2)

abline(h = - log10(0.05), lty = 2)