
library(GEOquery)
library(devtools)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5)
setwd("path/to//dataset/Richard_NatureGenetics_2015/")
options('download.file.method.GEOquery' = 'libcurl')# set this options when can't open the website

#gse1 <- GEOquery::getGEO("GSE71729",destdir=wd,GSEMatrix=TRUE)
#sub_gse <- gse1[[2]] # prostate tumor/normal exon sample in list2 annotated with GLPdata

sub_gse <- getGEO(filename="./file/GSE71729_series_matrix.txt.gz")

# convert epressionSet gse2 to probe expression matrix
gse_expr_string <- exprs(sub_gse)
gse_expr <- apply(gse_expr_string, 2 ,as.numeric)
rownames(gse_expr) <- rownames(gse_expr_string)
dim(gse_expr)
# subset sample into different group
sam.group.infor <- read.table("sample_group_infor.txt",header=F,sep="\t")
sam.group.infor <- sam.group.infor[,-2]
colnames(sam.group.infor) <- c("ID","group","tissue")

dat.met <- sam.group.infor[sam.group.infor$group=="Met",]
dim(dat.met)

dat.primary.PDAC <- sam.group.infor[sam.group.infor$group=="Primary",]
dim(dat.primary.PDAC)

dat.normal.PDAC <- sam.group.infor[sam.group.infor$group=="Normal" & sam.group.infor$tissue=="Pancreas",]
dim(dat.normal.PDAC)

dat.normal.other <- sam.group.infor[sam.group.infor$group=="Normal" & sam.group.infor$tissue!="Pancreas",]
dim(dat.normal.other)

library(ggplot2)
library(ggpubr) 
library(gridExtra)
all.sample <- colnames(gse_expr)
setwd("path/to//plot/Richard_NatureGenetics_2015")
#met, primaty, normal.PDAC, normal.other
for(i in 1:4){
  if(i==1){
    subsam <- dat.met[,1]
    group <- "metastatic.PDAC_"
  }else if(i==2){
    subsam <- dat.primary.PDAC[,1]
    group <- "primary.PDAC_"
  }else if(i==3){
    subsam <- dat.normal.PDAC[,1]
    group <- "pancreas.normal_"
  }else{
    subsam <- dat.normal.other[,1]
    group <- "other.normal_"
  }
  dat.expr <- gse_expr[,all.sample %in% subsam]
  expt.trans <- t(dat.expr)
  all.expr.gene <- colnames(expt.trans)
  dat.expr[,1] <- rownames(dat.expr)
  colnames(dat.expr)[1] <- "GeneName"
  #write.table(dat.expr,paste0(group,"expt.data.txt"),sep="\t",
  #            quote=F,row.names = F,col.names = T)
  
  
  NKCell <- read.table("path/to//geneset/NK_Cell_Signatures/NK_Cell_Signatures.txt",
                       sep="\t",header=T)
  TCell <- read.table("path/to//geneset/T_Cell_Signatures/T_Cell_Signatures.txt",
                      sep="\t",header=T)
  EZH2_genesets <- read.table("path/to//geneset/EZH2_Genesets/EZH2_Geneset.txt",
                              sep="\t",header=T)
  NF_KB <- read.table("path/to//geneset/NF-KB/NF-KB Signature.txt",
                      sep="\t",header=T)
  SASP <- read.table("path/to//geneset/Senescence/SASP.txt",
                     sep="\t",header=T)
  Interferon <- read.table("path/to//geneset/Interferon/Interferon_signature.txt",
                           sep="\t",header=T)
  IRF8 <- read.table("path/to//geneset/IRF8/IRF8_Signature.txt",
                     sep="\t",header=T)
  
  indi.gene.set12 <- c("IRF1", "IRF2", "IRF3", "IRF7", "IRF8","CD274","CCL2","CCL7", "CCL8", "IL15", "CXCL9", "CXCL10")
  # print gene not list in expression dataset
  if(any(!indi.gene.set12 %in% all.expr.gene)){
    warning(paste0(indi.gene.set12[!indi.gene.set12 %in% all.expr.gene],"is not in the expression dataset"))
  }
  indi.gene.set12 <- indi.gene.set12[indi.gene.set12 %in% all.expr.gene]
  indi.gene12 <- t(as.data.frame(indi.gene.set12))
  colnames(indi.gene12) <- indi.gene12
  
  indi.gene.set6 <- c("CCL2","CCL7", "CCL8", "IL15", "CXCL9", "CXCL10")
  # print gene not list in expression dataset
  if(any(!indi.gene.set6 %in% all.expr.gene)){
    warning(paste0(indi.gene.set6[!indi.gene.set6 %in% all.expr.gene],"is not in the expression dataset"))
  }
  indi.gene.set6 <- indi.gene.set6[indi.gene.set6 %in% all.expr.gene]
  indi.gene6 <- t(as.data.frame(indi.gene.set6))
  colnames(indi.gene6) <- indi.gene6
  
  
  indi.gene.set7 <- c("CD274","CCL2","CCL7","CCL8","IL15", "CXCL9", "CXCL10")
  # print gene not list in expression dataset
  if(any(!indi.gene.set7 %in% all.expr.gene)){
    warning(paste0(indi.gene.set7[!indi.gene.set7 %in% all.expr.gene],"is not in the expression dataset"))
  }
  indi.gene.set7 <- indi.gene.set7[indi.gene.set7 %in% all.expr.gene]
  indi.gene7 <- t(as.data.frame(indi.gene.set7))
  colnames(indi.gene7) <- indi.gene7
  
  #vs1.data.type="IRF8"
  for(vs1.data.type in c("NK Cell Signature","T Cell Signature","EZH2 Genesets",
                         "EZH2","SUZ12","IRF8 Signature","IRF8","NF-KB signatures")){
    
    vs2.set <- list()
    vs2.type <- list()
    
    ##::::::::::::::::::::::::::::::
    # nk sigmature
    ##::::::::::::::::::::::::::::::
    if(vs1.data.type == "NK Cell Signature"){
      vs1.data <- NKCell
      vs2.set[1] <- list(SASP)
      vs2.set[2] <- list(indi.gene6)
      
      vs2.type[1] <- list("SASP")
      vs2.type[2] <- list("indi.gene.set6")
    }
    
    ##::::::::::::::::::::::::::::::
    # T Cell sigmature
    ##::::::::::::::::::::::::::::::
    if(vs1.data.type=="T Cell Signature"){
      vs1.data <- TCell
      vs2.set[1] <- list(SASP)
      vs2.set[2] <- list(indi.gene6)
      
      vs2.type[1] <- list("SASP")
      vs2.type[2] <- list("indi.gene.set6")
    }
    
    ##::::::::::::::::::::::::::::::
    # EZH2 genesets
    ##::::::::::::::::::::::::::::::
    #NK cell signatures, b) T cell signatures, c) interferon signatures, d) IRF8 signatures, e) SASP and senescence signatures, f) NF-KB signatures, g) IRF1, IRF1, IRF2, IRF3, IRF7, IRF8, PDL-1, CCL2, CCL7, CCL8, IL-15, CXCL9, CXCL10 individual gene expression
    
    if(vs1.data.type=="EZH2 Genesets"){
      vs1.data <- EZH2_genesets
      vs2.set[1] <- list(NKCell)
      vs2.set[2] <- list(TCell)
      vs2.set[3] <- list(Interferon)
      vs2.set[4] <- list(IRF8)
      vs2.set[5] <- list(SASP)
      vs2.set[6] <- list(NF_KB)  
      vs2.set[7] <- list(indi.gene12) 
      
      vs2.type[1] <- list("NK Cell Signature")
      vs2.type[2] <- list("T Cell Signature")
      vs2.type[3] <- list("Interferon")
      vs2.type[4] <- list("IRF8 Signature")
      vs2.type[5] <- list("SASP")
      vs2.type[6] <- list("NF_KB")
      vs2.type[7] <- list("indi.gene.set12")
    }
    
    ##::::::::::::::::::::::::::::::
    # EZH2
    ##::::::::::::::::::::::::::::::
    if(vs1.data.type=="EZH2"){
      vs1.data <- as.data.frame("EZH2")
      vs2.set[1] <- list(NKCell)
      vs2.set[2] <- list(TCell)
      vs2.set[3] <- list(Interferon)
      vs2.set[4] <- list(IRF8)
      vs2.set[5] <- list(SASP)
      vs2.set[6] <- list(NF_KB)  
      vs2.set[7] <- list(indi.gene12)
      
      vs2.type[1] <- list("NK Cell Signature")
      vs2.type[2] <- list("T Cell Signature")
      vs2.type[3] <- list("Interferon")
      vs2.type[4] <- list("IRF8 Signature")
      vs2.type[5] <- list("SASP")
      vs2.type[6] <- list("NF_KB")
      vs2.type[7] <- list("indi.gene.set12")
    }
    
    ##::::::::::::::::::::::::::::::
    # SUZ12
    ##::::::::::::::::::::::::::::::
    if(vs1.data.type=="SUZ12"){
      vs1.data <- as.data.frame("SUZ12")
      vs2.set[1] <- list(NKCell)
      vs2.set[2] <- list(TCell)
      vs2.set[3] <- list(Interferon)
      vs2.set[4] <- list(IRF8)
      vs2.set[5] <- list(SASP)
      vs2.set[6] <- list(NF_KB)  
      vs2.set[7] <- list(indi.gene12)
      
      vs2.type[1] <- list("NK Cell Signature")
      vs2.type[2] <- list("T Cell Signature")
      vs2.type[3] <- list("Interferon")
      vs2.type[4] <- list("IRF8 Signature")
      vs2.type[5] <- list("SASP")
      vs2.type[6] <- list("NF_KB")
      vs2.type[7] <- list("indi.gene.set12")
    }
    
    ##::::::::::::::::::::::::::::::
    # IRF8
    ##::::::::::::::::::::::::::::::
    if(vs1.data.type=="IRF8"){
      vs1.data <- as.data.frame("IRF8")
      vs2.set[1] <- list(NKCell)
      vs2.set[2] <- list(TCell)
      vs2.set[3] <- list(SASP)
      vs2.set[4] <- list(indi.gene7)
      
      vs2.type[1] <- list("NK Cell Signature")
      vs2.type[2] <- list("T Cell Signature")
      vs2.type[3] <- list("SASP")
      vs2.type[4] <- list("indi.gene.set7")
    }
    
    ##::::::::::::::::::::::::::::::
    # IRF8 Signature
    ##::::::::::::::::::::::::::::::
    if(vs1.data.type=="IRF8 Signature"){
      vs1.data <- IRF8
      vs2.set[1] <- list(NKCell)
      vs2.set[2] <- list(TCell)
      vs2.set[3] <- list(SASP)
      vs2.set[4] <- list(indi.gene7) 
      
      vs2.type[1] <- list("NK Cell Signature")
      vs2.type[2] <- list("T Cell Signature")
      vs2.type[3] <- list("SASP")
      vs2.type[4] <- list("indi.gene.set7")
    }
    
    ##::::::::::::::::::::::::::::::
    # NF-KB signatures
    ##::::::::::::::::::::::::::::::
    if(vs1.data.type=="NF-KB signatures"){
      vs1.data <- NF_KB
      vs2.set[1] <- list(NKCell)
      vs2.set[2] <- list(TCell)
      vs2.set[3] <- list(SASP)
      vs2.set[4] <- list(indi.gene7) 
      
      vs2.type[1] <- list("NK Cell Signature")
      vs2.type[2] <- list("T Cell Signature")
      vs2.type[3] <- list("SASP")
      vs2.type[4] <- list("indi.gene.set7")
    }
    
    #vs2.index=1
    for(vs2.index in 1:length(vs2.set)){
      vs2 <- vs2.set[[vs2.index]]
      ## correlation plot
      p <- list()
      vs1.data.not.list.dataset <- list()
      #i=1
      for(i in 1:ncol(vs1.data)){
        vs1.data.sub <- vs1.data[,i]
        vs1.data.sub <- toupper(vs1.data.sub)
        vs1.data.sub.rmNA <- vs1.data.sub[!is.na(vs1.data.sub)]
        d1 <- data.frame(expt.trans[,colnames(expt.trans) %in% vs1.data.sub.rmNA])
        vs1.data.not.list.dataset[colnames(vs1.data)[i]] <- list(colnames(d1))
        d1.rowmean <- data.frame(rowMeans(d1))
        for(j in 1:ncol(vs2)){
          #j=1
          vs2.sub <- vs2[,j]
          vs2.sub <- toupper(vs2.sub)
          vs2.sub.rmNA <- vs2.sub[!is.na(vs2.sub)]
          d2 <- expt.trans[,colnames(expt.trans) %in% vs2.sub.rmNA]
          d2 <- as.data.frame(d2)
          d2.rowmean <- data.frame(rowMeans(d2))
          corrdat <- data.frame(d2.rowmean,d1.rowmean)
          if(ncol(vs1.data)==1){
            yname=vs1.data.type
          }else{
            yname=paste0(vs1.data.type,i)
          }
          colnames(corrdat) <- c(colnames(vs2)[j],yname)
          
          # Scatter plot with correlation coefficient
          
          if(i==ncol(vs1.data)){
            if(j == 1){
              sp <- ggscatter(corrdat, x = colnames(vs2)[j], y = yname,
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
              )
            }else{
              sp <- ggscatter(corrdat, x = colnames(vs2)[j], y = yname,
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE, # Add confidence interval
                              ylab = ""
              )
            }
          }else{
            if(j == 1){
              sp <- ggscatter(corrdat, x = colnames(vs2)[j], y = yname,
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE, # Add confidence interval
                              xlab=""
              )
            }else{
              sp <- ggscatter(corrdat, x = colnames(vs2)[j], y = yname,
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE, # Add confidence interval
                              ylab = "",
                              xlab=""
              )
            }
          }
          
          # Add correlation coefficient
          p[(i-1)*ncol(vs2)+j] <- list(sp + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),method = "pearson"))
        }
      }
      
      if(length(p)>0){
        nocol <- length(p)/ncol(vs1.data)
        pdf(file = paste0(getwd(),"/plot_R/",group,vs1.data.type,"_VS_",vs2.type[[vs2.index]],".pdf"),height = 5*ncol(vs1.data),width = 5*nocol)
        print(do.call('grid.arrange',c(p, nrow = ncol(vs1.data))))
        dev.off()
      }
    }
  }
}












