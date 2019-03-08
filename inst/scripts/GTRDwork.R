
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

save(gene_list,file="data/gene_df.RData")




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
  temp=gene_list$tf[which(gene_list$pro_gene==i)]
  ex_set[[i]]$tf=factor(temp)
}
ex_set

df=data.frame(gene=c(),num=c())
for(i in xSet){
  if(length(ex_set[[i]]$tf)>=1){
    temp=data.frame(i,num=length(ex_set[[i]]$tf))
    df=rbind(df,temp)
  }
  
}

library(ggplot2)
p<-ggplot(data=df, aes(x=i, y=num))+
  geom_bar(stat="identity")+
  theme(legend.direction = "vertical")+
  theme(axis.text.x = element_text(angle = -90))
+theme(legend.position = "bottom") 
p
barplot(df$num,
        ylim=c(0,100), breaks = 10,col = "#B22222",
        main = "Number of transcription factor of gene distribution",
        xlab = "gene name",
        ylab = "Counts",names.arg =df$i )



#[END]

