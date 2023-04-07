---
title: "GFF"
author: "Jose V. Die"
date: "4/7/2023"
output: 
  html_document:
    keep_md: true
---



## GFF

The General Feature Format ([GFF](https://www.ensembl.org/info/website/upload/gff.html)) file has nine columns:  
  
  
| colname | description |
|---|---|
| seqname | The name of the sequence; must be a chromosome or scaffold. |
| source | The program that generated this feature. |
| feature | The name of this type of feature, e.g. "CDS", "start_codon", "stop_codon", and "exon" |
| start | The starting position of the feature in the sequence; the first base is numbered 1. |
| end | The ending position of the feature (inclusive). |
| score | A score between 0 and 1000. |
| strand | Valid entries include "+", "-", or ".". |
| frame | If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be ".". |
  
Dependencies

```r
library(rtracklayer)
```

### Load GFF file 

```r
gff_file <- "../dat/Saccharomyces_cerevisiae_genome.gff3"
```
Available molecular features in the GFF file 


```r
gff_all = import.gff(con = gff_file ,
                      colnames = c("type", "mycolums"))
table(gff_all$type)
```

```
## 
##                chromosome                      gene                      mRNA 
##                        17                      6600                      6600 
##                      exon                       CDS                ncRNA_gene 
##                      7507                      6913                       424 
##                     ncRNA                      tRNA                    snoRNA 
##                        18                       299                        77 
## transposable_element_gene      transposable_element                pseudogene 
##                        91                        91                        12 
##    pseudogenic_transcript                     snRNA            five_prime_UTR 
##                        12                         6                         4 
##                      rRNA 
##                        24
```

Load gene models from the GFF file. 

```r
gff_chromosome = import.gff(con = gff_file ,
                     feature.type = "chromosome",
                     colnames = c("type", "mycolums"))

gff_gene = import.gff(con = gff_file ,
                     feature.type = "gene",
                     colnames = c("type", "mycolums"))

gff_exon = import.gff(con = gff_file ,
                     feature.type = "exon",
                     colnames = c("type", "mycolums"))
```

Look the data

```r
gff_chromosome
```

```
## GRanges object with 17 ranges and 2 metadata columns:
##        seqnames    ranges strand |       type    mycolums
##           <Rle> <IRanges>  <Rle> |   <factor> <character>
##    [1]        I  1-230218      * | chromosome        <NA>
##    [2]       II  1-813184      * | chromosome        <NA>
##    [3]      III  1-316620      * | chromosome        <NA>
##    [4]       IV 1-1531933      * | chromosome        <NA>
##    [5]       IX  1-439888      * | chromosome        <NA>
##    ...      ...       ...    ... .        ...         ...
##   [13]      XII 1-1078177      * | chromosome        <NA>
##   [14]     XIII  1-924431      * | chromosome        <NA>
##   [15]      XIV  1-784333      * | chromosome        <NA>
##   [16]       XV 1-1091291      * | chromosome        <NA>
##   [17]      XVI  1-948066      * | chromosome        <NA>
##   -------
##   seqinfo: 17 sequences from an unspecified genome; no seqlengths
```

```r
gff_gene
```

```
## GRanges object with 6600 ranges and 2 metadata columns:
##          seqnames        ranges strand |     type    mycolums
##             <Rle>     <IRanges>  <Rle> | <factor> <character>
##      [1]        I       335-649      + |     gene        <NA>
##      [2]        I       538-792      + |     gene        <NA>
##      [3]        I     1807-2169      - |     gene        <NA>
##      [4]        I     2480-2707      + |     gene        <NA>
##      [5]        I     7235-9016      - |     gene        <NA>
##      ...      ...           ...    ... .      ...         ...
##   [6596]      XVI 939922-941136      + |     gene        <NA>
##   [6597]      XVI 943032-943896      + |     gene        <NA>
##   [6598]      XVI 943880-944188      + |     gene        <NA>
##   [6599]      XVI 944603-947701      + |     gene        <NA>
##   [6600]      XVI 946856-947338      - |     gene        <NA>
##   -------
##   seqinfo: 17 sequences from an unspecified genome; no seqlengths
```

```r
gff_exon
```

```
## GRanges object with 7507 ranges and 2 metadata columns:
##          seqnames        ranges strand |     type    mycolums
##             <Rle>     <IRanges>  <Rle> | <factor> <character>
##      [1]        I       335-649      + |     exon        <NA>
##      [2]        I       538-792      + |     exon        <NA>
##      [3]        I     1807-2169      - |     exon        <NA>
##      [4]        I     2480-2707      + |     exon        <NA>
##      [5]        I     7235-9016      - |     exon        <NA>
##      ...      ...           ...    ... .      ...         ...
##   [7503]      XVI 943032-943050      + |     exon        <NA>
##   [7504]      XVI 943199-943896      + |     exon        <NA>
##   [7505]      XVI 943880-944188      + |     exon        <NA>
##   [7506]      XVI 944603-947701      + |     exon        <NA>
##   [7507]      XVI 946856-947338      - |     exon        <NA>
##   -------
##   seqinfo: 17 sequences from an unspecified genome; no seqlengths
```

Genes on chr 'I' and strand 'forward'

```r
gff_gene[seqnames(gff_gene) == "I" & strand(gff_gene) == "+"]
```

```
## GRanges object with 60 ranges and 2 metadata columns:
##        seqnames        ranges strand |     type    mycolums
##           <Rle>     <IRanges>  <Rle> | <factor> <character>
##    [1]        I       335-649      + |     gene        <NA>
##    [2]        I       538-792      + |     gene        <NA>
##    [3]        I     2480-2707      + |     gene        <NA>
##    [4]        I   10091-10399      + |     gene        <NA>
##    [5]        I   12046-12426      + |     gene        <NA>
##    ...      ...           ...    ... .      ...         ...
##   [56]        I 221049-221660      + |     gene        <NA>
##   [57]        I 222406-222891      + |     gene        <NA>
##   [58]        I 225460-226863      + |     gene        <NA>
##   [59]        I 227742-228953      + |     gene        <NA>
##   [60]        I 228844-229317      + |     gene        <NA>
##   -------
##   seqinfo: 17 sequences from an unspecified genome; no seqlengths
```

