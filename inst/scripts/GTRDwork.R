
# fetch HGNC reference data
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame

chr_list<-c(1:22,"X","Y","MT")
file<-"data/Homo_sapiens_meta_clusters.interval"
file.rename(file,paste(file,".txt",sep=""))
chr_data <- readr::read_tsv("data/Homo_sapiens_meta_clusters.interval.txt",
                            col_names = c("chr","start","end","summit",
                                          "UniProt","sym"),skip = 1,
                            col_types = "ciiicc-------")


#check correction of the peak
indentifier<-c()
for(i in 1:nrow(chr_data)){
  x<-c(as.integer(chr_data[i,"start"]),as.integer(chr_data[i,"end"]))
  v<-as.integer(chr_data[i,"start"])+as.integer(chr_data[i,"summit"])
  indentifier<-c(indentifier,sum(findInterval(x, v)) ==1)
  
}

# all peak summit within the interval
all(indentifier) #all peaks fine. 


#find outdated gene symbol
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



#tf site for each TF
tfSites <- list()
# initialize all transcription factors

myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
#get all transcripts of gene sybol 
transcript_data2 <- biomaRt::getBM(filters = "hgnc_symbol",
                                   attributes = c("hgnc_symbol",
                                                  "transcript_appris",
                                                  "chromosome_name",
                                                  "transcript_start",
                                                  "transcript_end"),
                                   values = HGNC$sym,
                                   mart = myMart)
#remove duplicated row based on start position, chromosome, gene symbol 
transcript_data<-unique(transcript_data[,c(1,3,4)])
transcript_data$hgnc_symbol[1]
head(transcript_data$hgnc_symbol)


#annotation
transcript_data$transcript_start<-transcript_data$transcript_start-1001
head(transcript_data$transcript_start)
gene_list<-list()  #sym -> TF
for(chr in chr_list){
  transcrip_subset_index<-which(transcript_data$chromosome_name==chr)
  chr_name<-paste ("chr",chr, sep = "", collapse = NULL)
  #for 1 chromosome TTS
  nrow(transcript_data$transcript_start)
  promo_reg_start<-transcript_data$transcript_start[transcrip_subset_index]
  promo_reg_end<-promo_reg_start+1000
  TFs<-unique(chr_data$sym[chr_data$chr==chr_name])
  
  for(tf in TFs){
    index<-which(chr_data$sym==tf&chr_data$chr==chr_name)
    site_tf<-chr_data$start[index]+chr_data$summit[index]
    for(site in site_tf){
      check<-(promo_reg_start<=site & site<=promo_reg_end)
      #check binding gene promoter region
      if(any(check)){
        index_index<-which(check)
        index_bind<-transcrip_subset_index[index_index]
        pro_gene<-transcript_data$hgnc_symbol[index_bind]
        gene_list[[tf]]$genes<-unique(c(gene_list[[tf]]$genes,pro_gene))
      }
      
    }
    
  }
}

save(gene_list,file="data/gene_list.RData")

tf_list<-list()
for(gene_sym in HGNC$sym){
  for(tf_sym in names(gene_list)){
    if(gene_sym %in% gene_list[[tf_sym]]$genes){
      tf_list[[gene_sym]]$tfs<-c(tf_list[[gene_sym]]$tfs,tf_sym)
    }
  }
}

for (i in tf_list) {
  i$tfs<-unique(i$tfs)
  
}

save(tf_list,file="data/tf_List.RData")




#basic statistics
#number of TFs that binds to some gene's promoter region
length(gene_list)
#number of gene that each tf bind to their promotor region
length(tf_list)


#example annotation

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
  ex_set[[i]]$TFs<-tf_list[[i]]$tfs
}
ex_set

ex_genelist<-list()
#associated gene for a given gene
for(gene in ex_set){
  for(tf in ex_set[[gene]]$TFs){
    ex_genelist[[gene]]<-c(ex_genelist[[gene]],gene_list[[tf]]$genes)
  }
}




#[END]

