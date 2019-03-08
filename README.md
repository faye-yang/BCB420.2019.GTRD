# `BCB420.2019.GTRD`

#### (GTRD data annotatation of human genes)

&nbsp;

###### [Yufei Yang]

----

# 1 About this package:


This package describes the workflow to download binding sites of transcription factor that is  identified from ChIP-seq experiments from the GTRD database, how to clean up the dataset, how to annotate the example gene set, and provides examples of computing database statistics. Ultimately, given a gene symbol, the package can compute other symbols are regulated by the same transcription factors as the given gene.



```text
 --BCB420.2019.GTRD/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.GTRD.Rproj
   |__DESCRIPTION
   |__dev/
      |__rptTwee.R
      |__toBrowser.R               # display .md files in your browser
   |__inst/
      |__extdata/
         |__
      |__img/
         |__workflow.png                  # image sources for .md document
      |__scripts/
         |__GTRDwork.R           # get TFs bind to the given gene's promotor region 
   |__LICENSE
   |__NAMESPACE
   |__R/
      |__zzz.R
   |__README.md                    # this file

```


----

# 2 GTRD Data semantics
Reference from http://wiki.biouml.org/index.php/GTRD
Raw ChIP-seq data and experiment information were collected from:literature,GEO,SRA, ENCODE.
Initial ChIP-seq and DNase-seq raw data were uniformly processed using specially developed workflow (pipeline) for the BioUML platform. ChIP-seq processing pipeline included the following steps:
1. sequenced reads were aligned to the corresponding reference genome using Bowtie2;
peaks were identified using MACS, SISSR, GEM and PICS peak callers;
2. peaks computed for the same TF and peak calling method, but different experiment conditions (e.g., cell line, treatment, etc.) were joined into clusters;
3. clusters for the same TF revealed by different peak calling methods were joined into metaclusters.
DNase-seq processing pipeline included the following steps:
1.sequenced reads were aligned to the corresponding reference genome using Bowtie2;
2.regions of open chromatin were identified using MACS2 and Hotspot2;
3. de novo putative protein-DNA interactions were revealed using a digital genomic footprinting tool Wellington.

GTRD is a database of f transcription factor (TF) binding sites identified from ChIP-seq experiments that were systematically collected and uniformly processed using a special workflow (pipeline) for a BioUML platform . 
![GTRD Workflow](/inst/img/workflow.png)
(http://www.biouml.org)
Chromosomal coordinates reference: Under Data processing workflow: alignment of reads—we used Bowtie2 (version 2.2.3) (12) to align ChIP-seq reads to the reference human (GRCh38) and mouse (GRCm38) genomes. (Reference:(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5210645/))



&nbsp;

# 3 Data download and cleanup

1. Navigate to the [**GTRD** database](http://gtrd.biouml.org/) and follow the link to the [download section](http://gtrd.biouml.org/downloads/18.06/).
2. Choose "Homo sapiens_meta_clusters_interval". (only download one data file since the database is too large)
3. Download data files. Warning: large.

* `Homo sapiens_meta_clusters_interval.gz` (725M)	gene binding site data (postions where TF bind);

4. Uncompress the files and place them into a sister directory of your working directory which is called `data`. (It should be reachable with `file.path("..", "data")`). **Warning:**  `../data/Homo sapiens_meta_clusters_interval` is 7.68 GB;
5. Keep only important information: chromosome number, position start, position end, TF title(gene symbol).
6. Only keep one chromosome to analyze since the data set is too lage. Using for loop, all chromosomes can be analyzed.
7. Remove duplicated row and update outdated gene symbol.
&nbsp;
```R

#check correction of the peak
indentifier<-c()
for(i in 1:nrow(chr_data)){
  x<-c(as.integer(chr_data[i,"start"]),as.integer(chr_data[i,"end"]))
  v<-as.integer(chr_data[i,"start"])+as.integer(chr_data[i,"summit"])
  indentifier<-c(indentifier,sum(findInterval(x, v)) ==1)
  
}

# all peak summit within the interval
all(indentifier) #all peaks fine are regulatory region
#remove peak that
chr_data<-chr_data[indentifier,]
```

Find outdated gene symbol
```R
sel <- ( ! (chr_data$sym %in% HGNC$sym)) & ( ! (is.na(chr_data$sym)))
sum(sel) 
uSym<-unique(chr_data$sym) 

unkSym <- data.frame(unk = chr_data$sym[sel],
                     new = NA,
                     stringsAsFactors = FALSE)
# grep() for the presence of the symbols in either HGNC$prev or
# HGNC$synonym. If either is found, that symbol replaces NA in unkSym$new
for (i in seq_len(nrow(unkSym))) {
  iPrev <- grep(unkSym$unk[i], HGNC$prev)[1] # take No. 1 if there are several
  if (length(iPrev) == 1) {
    unkSym$new[i] <- HGNC$sym[iPrev]
  } else {
    iSynonym <- which(grep(unkSym$unk[i], HGNC$synonym))[1]
    if (length(iSynonym) == 1) {
      unkSym$new[i] <- HGNC$sym[iSynonym]
    }
  }
}
sum( unique(is.na(unkSym$new)))#1  
sum( !is.na(unkSym$new))#51486
indexNA<-which(!is.na(unkSym$new))
#update the outdated gene symbol 
chr_data$sym[sel] <-unkSym$new

sum(is.na(chr_data$sym) ) #0
#do not need to remove rows that have NA
```

Retrieve Transcript information
&nbsp;
```R 
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
#get all transcripts of gene sybol 
transcript_data <- biomaRt::getBM(filters = "hgnc_symbol",
                                  attributes = c("hgnc_symbol",
                                                 "transcript_appris",
                                                 "chromosome_name",
                                                 "transcript_start",
                                                 "transcript_end"),
                                  values = HGNC$sym,
                                  mart = myMart)



```


# 4 A script for annotating gene sets

 ```R
#annotation
gene_list<-data.frame(gene=c(),TF=c(),stringsAsFactors = FALSE)  #sym -> TF
for(chr in chr_list){
  transcrip_subset_index<-which(transcript_data$chromosome_name==chr)
  chr_name<-paste ("chr",chr, sep = "", collapse = NULL)
  #for 1 chromosome TTS
  TFs<-unique(chr_data$sym[chr_data$chr==chr_name])
  
  for(tf in TFs){
    index<-which(chr_data$sym==tf&chr_data$chr==chr_name)
    site_tf<-chr_data$start[index]+chr_data$summit[index]
    for(site in site_tf){
      index1=which((transcript_data$transcript_start-1001)<=site)
      index2=which(site<=transcript_data$transcript_start)
      inter=intersect(index2,index1)
      if(length(inter)!=0){
        pro_gene<-transcript_data$hgnc_symbol[inter]
        temp=data.frame(pro_gene,tf)
        gene_list=rbind(gene_list, temp)
        print(gene_list)
        
      }
    }
  }
}
gene_list=unique(gene_list)
unique_tf=unique(gene_list$tf)
unique_gene=unique(gene_list$pro_gene)


#basic statistics
#number of TFs that binds to some gene's promoter region
length(gene_list)
#number of gene that each tf bind to their promotor region
length(tf_list)
 ```


# 5 Annotating the example set
```R
# The specification of the sample set is copy-paste from the 
# BCB420 resources project.

xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")


ex_set<-list()
#associated TF
for (i in xSet) {
  temp=gene_list$tf[which(gene_list$pro_gene==i)]
  ex_set[[i]]$tf=factor(temp)
}
ex_set
#$TPPP
#$TPPP$TFs
# [1] "MYC"   "ESR1"  "PGR"   "PPARG" "CEBPD" "MAX"   "FLI1"  "SP2"   "RELA"  "SMC3"
#ommiting other results

ex_genelist<-list()
#associated gene for a given gene
for(gene in ex_set){
  for(tf in ex_set[[gene]]$TFs){
    ex_genelist[[gene]]<-c(ex_genelist[[gene]],gene_list[[tf]]$genes)
  }
}

$DCTN1
$DCTN1$tf
 [1] MEIS2   ZNF263  CLOCK   NCOR1   MYC     ESR1    MYCN    TP53    JUN    
[10] PGR     SP1     ZBTB48  AR      EGR2    ERG     ETS1    ATF2    CREB1  
[19] JUNB    JUND    TAL1    ATF7    SPI1    ATF1    ATF3    WT1     RXRA   
[28] TBP     TAF1    USF1    RFX1    NFYA    GATA2   YY1     ARNT    OTX2   
[37] PPARG   STAT5A  CTCF    MXI1    STAT5B  ZNF143  FOXA1   MAZ     MAX    
[46] E2F1    RUNX1   FLI1    POU5F1  SPIB    SP2     SP4     CREM    KMT2A  
[55] RELA    EGR3    FOXM1   ZBTB17  SP140   NEUROD1 USF2    CRY1    NRF1   
[64] ZNF554  GLIS1   GRHL3   TFAP2C  TCF12   SCRT1   ZFHX2   NANOG   TCF7L1 
[73] PRDM9   ZNF629  ZNF639 

$MAP1LC3C
$MAP1LC3C$tf
 [1] FOS    NR3C1  JUN    ZBTB48 FOSL1  FOSL2  GATA1  CREB1  ZNF30  JUND   TAL1  
[12] CEBPB  SPI1   ATF3   RXRA   YY1    MZF1   STAT5A NKX2-1 CEBPA  STAT5B ZNF143
[23] FOSB   MAZ    MAX    RUNX1  FLI1   RELA   KLF9   NFE2   ZNF770 ZNF549 ZNF600
[34] ZNF554 GLIS1  RBAK  

$STX6
$STX6$tf
 [1] MEIS2   ZNF263  ZEB2    E2F6    ZBTB7A  FOS     MYC     ESR1    NR3C1  
[10] PGR     AR      NR2F6   NR2F1   ERG     VDR     ESRRA   IRF2    GATA1  
[19] JUNB    JUND    TAL1    CEBPB   SPI1    WT1     RXRA    BCL3    NR4A1  
[28] GATA2   NR2F2   SOX6    PPARG   PBX3    NR2C2   CEBPA   CEBPD   ASCL1  
[37] ZNF143  FOXA1   PKNOX1  MAX     FOXK2   RUNX1   ID3     CREM    KMT2A  
[46] RELA    MTA1    IKZF1   HES1    TEAD4   TWIST1  SMAD1   NFE2    ZNF766 
[55] ARID1B  MNT     NEUROG2 TCF7L2  KDM5B   SMC3    FOXA2 

```

&nbsp;

# 6 Further reading

* Yevshin I, Sharipov R, Valeev T, Kel A, Kolpakov F. GTRD: a database of transcription factor binding sites identified by ChIP-seq experiments. Nucleic Acids Res. 2016;45(D1):D61-D67.
* Ivan Yevshin, Ruslan Sharipov, Semyon Kolmykov, Yury Kondrakhin, Fedor Kolpakov; GTRD: a database on gene transcription regulation—2019 update, Nucleic Acids Research, Volume 47, Issue D1, 8 January 2019, Pages D100–D105

&nbsp;

# 7 Acknowledgements

Thanks to Simon Kågedal's [PubMed to APA reference tool](http://helgo.net/simon/pubmed/).



&nbsp;

<!-- END -->
