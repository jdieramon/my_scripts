# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub", version = "3.8")


# Load library
library(AnnotationHub)

help("AnnotationHub")


## create an AnnotationHub object
ah = AnnotationHub()

## Summary of available records
ah

## Detail for a single record
ah[1]

## and what is the date we are using?
snapshotDate(ah)

## how many resources?
length(ah)

## from which resources, is data available?
head(sort(table(ah$dataprovider), decreasing=TRUE))

## from which species, is data available ? 
head(sort(table(ah$species),decreasing=TRUE)) 

## what web service and local cache does this AnnotationHub point to?
hubUrl(ah)
hubCache(ah)

## search the hub 
cicer <- query(ah, "Cicer arietinum")

# extract information
cicer$dataprovider
cicer$genome
cicer$maintainer

# "AnnotationHub" object
class(cicer)

# retrieve (=download) the element
cicer <- cicer[[1]]

#OrgDb object
class(cicer)
cicer

# accessing data :
## method columns
columns(cicer)

#### method keytypes
keytypes(cicer)

## method keys
head(keys(cicer, keytype = "GENENAME"))
head(keys(cicer, keytype = "ENTREZID"))
head(keys(cicer, keytype = "SYMBOL"))
head(keys(cicer, keytype = "SYMBOL", pattern = "ARF"))


# method select
select(cicer, keys = "101509359", keytype = "ENTREZID", 
       columns = c("GO", "ONTOLOGY", "SYMBOL"))

select(cicer, keys = "LOC101509359", keytype = "ALIAS", 
       columns = c("GO", "CHR", "GENENAME"))


library(dplyr)
select(cicer, keys = "101509359", keytype = "ENTREZID", 
       columns = columns(cicer)) %>% 
  display() 


