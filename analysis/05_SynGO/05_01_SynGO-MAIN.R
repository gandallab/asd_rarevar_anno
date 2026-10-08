# Collaborator script (Bert, CMU), rerun in place here. Only the file paths differ from
# the version delivered to us: inputs now come from this project's data/ and outputs are
# written to outputs/analysis/05_SynGO/ (tables) and outputs/figures/supplement/ (plots).
rm(list=ls()); gc()
options(stringsAsFactors = F)

project_dir <- here::here()
data_dir     <- file.path(project_dir, "data")
syngo_dir    <- file.path(data_dir, "SynGO")
out_dir      <- file.path(project_dir, "outputs/analysis/05_SynGO")
supp_fig_dir <- file.path(project_dir, "outputs/figures/supplement")

require(data.table)
require(readxl)
require(WriteXLS)

# syngo data
syngo.genes=as.data.frame(read_excel(file.path(syngo_dir,"syngo_genes.xlsx")))  # 1602
rownames(syngo.genes)=syngo.genes$ensembl_id

syngo.annot=as.data.frame(read_excel(file.path(syngo_dir,"syngo_annotations.xlsx")))  # 4218
syngo.annot$ensembl_id=syngo.genes$ensembl_id[match(syngo.annot$hgnc_id,syngo.genes$hgnc_id)]

syngo.ontol=as.data.frame(read_excel(file.path(syngo_dir,"syngo_ontologies.xlsx")))  # 302
rownames(syngo.ontol)=syngo.ontol$id

# process syngo data
# create a list of unique genes for each compartment
syngo.ontol[,c("all.genes","child.genes","unique.genes")]=NA
for(compartment in rownames(syngo.ontol)){
  all.genes=unique(unlist(strsplit(syngo.ontol[compartment,"ensembl_id"],", ",fixed=T)))
  child.genes=unique(unlist(strsplit(syngo.ontol[which(syngo.ontol[,"parent_id"] == compartment),"ensembl_id"],", ",fixed=T)))
  unique.genes=setdiff(all.genes,child.genes)
  
  syngo.ontol[compartment,"all.genes"]=paste(all.genes,collapse=";")
  syngo.ontol[compartment,"child.genes"]=paste(child.genes,collapse=";")
  syngo.ontol[compartment,"unique.genes"]=paste(unique.genes,collapse=";")
}
  


# asc data 
asc=fread(file.path(data_dir,"full_results_wcounts_2025-10-08.txt"),header=T,data.table=F)
asc=asc[,c("gene","gene_id","FDR","LOEUF","PTV_Proband","PTV_Sibling","Mis2_Proband","Mis2_Sibling","Mis1_Proband","Mis1_Sibling","Del_proband","Del_sibling","Dup_proband","Dup_sibling")]
rownames(asc)=asc$gene_id
# select the genes that are in the syngo files
asc=asc[rownames(asc) %in% syngo.genes$ensembl_id,] # 1533
# combined counts for both probands and siblings
asc$dnProband=rowSums(asc[,grep("proband",tolower(colnames(asc)))])
asc$dnSibling=rowSums(asc[,grep("sibling",tolower(colnames(asc)))])
colSums(asc[,-c(1:4)])

# PTV_Proband  PTV_Sibling Mis2_Proband Mis2_Sibling Mis1_Proband Mis1_Sibling  Del_proband  Del_sibling  Dup_proband  Dup_sibling    dnProband    dnSibling 
#        1031           94          683           42          760          119           78            4           30            6         2582          265 

# for each syngo compartment determine the count of denovo proband and sinling events
syngo.df=syngo.ontol[,c("id","domain","name","shortname","parent_id")]
syngo.df$is.parent=syngo.df$id %in% syngo.df$parent_id
syngo.df$orig.order=1:nrow(syngo.df)
rownames(syngo.df)=syngo.df$id

syngo.df[,c("n.all","dnPro.all","dnSib.all","n.chl","dnPro.chl","dnSib.chl","n.unq","dnPro.unq")]=NA
# add the counts
for(id in rownames(syngo.df)){
  
  # for all genes
  i=which(asc$gene_id %in% unique(unlist(strsplit(syngo.ontol[id,"all.genes"],";",fixed=T))))
  syngo.df[id,"n.all"]=length(i)
  syngo.df[id,"dnPro.all"]=sum(asc[i,"dnProband"])
  syngo.df[id,"dnSib.all"]=sum(asc[i,"dnSibling"])
  
  # for genes seen in children
  i=which(asc$gene_id %in% unique(unlist(strsplit(syngo.ontol[id,"child.genes"],";",fixed=T))))
  syngo.df[id,"n.chl"]=length(i)
  syngo.df[id,"dnPro.chl"]=sum(asc[i,"dnProband"])
  syngo.df[id,"dnSib.chl"]=sum(asc[i,"dnSibling"])
  
  # for genes in parents that are not seen in children
  i=which(asc$gene_id %in% unique(unlist(strsplit(syngo.ontol[id,"unique.genes"],";",fixed=T))))
  syngo.df[id,"n.unq"]=length(i)
  syngo.df[id,"dnPro.unq"]=sum(asc[i,"dnProband"])
  syngo.df[id,"dnSib.unq"]=sum(asc[i,"dnSibling"])
}

colnames(syngo.df)
nPro=38680 
nSib=9567
syngo.df$ratePro.all = syngo.df$dnPro.all/nPro
syngo.df$rateSib.all = (syngo.df$dnSib.all+0.01)/(nSib)
syngo.df$rateRatio.all = syngo.df$ratePro.all/syngo.df$rateSib.all
syngo.df$ratePro.chl = syngo.df$dnPro.chl/nPro
syngo.df$rateSib.chl = (syngo.df$dnSib.chl+0.01)/(nSib)
syngo.df$rateRatio.chl = syngo.df$ratePro.chl/syngo.df$rateSib.chl

propTestUsingZdist <- function(n1,n2,N1,N2){
  p1=n1/N1
  p2=n2/N2
  p=(n1+n2)/(N1+N2)
  Z=(p1-p2)/sqrt(p*(1-p)*(1/N1+1/N2))
  return(Z)
}
syngo.df[,c("Z.all","P.all")]=NA
for(i in 1:nrow(syngo.df)){
  syngo.df$Z.all[i]=propTestUsingZdist(n1=syngo.df$dnPro.all[i],n2=syngo.df$dnSib.all[i],N1=nPro,N2=nSib)
  syngo.df$P.all[i]=pnorm(syngo.df$Z.all[i],lower.tail=F)
}
syngo.df[,c("adjZ.all")]=NA
domain="CC"
summary(lm(Z.all~n.all,data=syngo.df[syngo.df$domain == domain,],na.action = "na.exclude"))$r.squared # 0.615
syngo.df$adjZ.all[syngo.df$domain == domain]=residuals(lm(Z.all~n.all,data=syngo.df[syngo.df$domain == domain,],na.action = "na.exclude"))

domain="BP"
summary(lm(Z.all~n.all,data=syngo.df[syngo.df$domain == domain,],na.action = "na.exclude"))$r.squared # 0.419
syngo.df$adjZ.all[syngo.df$domain == domain]=residuals(lm(Z.all~n.all,data=syngo.df[syngo.df$domain == domain,],na.action = "na.exclude"))

# add the average LOEUF
syngo.df$ave.LOEUF=NA
for(i in 1:nrow(syngo.df)){
  id=rownames(syngo.df)[i]
  ensg=unlist(strsplit(syngo.ontol[id,"all.genes"],";"))
  syngo.df$ave.LOEUF[i]=mean(asc[ensg,"LOEUF"],na.rm=T)
}

# rearrange the columns
syngo.df=syngo.df[,c("id","domain","name","shortname","parent_id","is.parent","orig.order","n.all","dnPro.all","dnSib.all","n.chl","dnPro.chl","dnSib.chl",
                    "n.unq","dnPro.unq","dnSib.unq","ratePro.all","rateSib.all","rateRatio.all","ratePro.chl","rateSib.chl","rateRatio.chl","Z.all","ave.LOEUF",
                    "P.all","adjZ.all","n.gen","SYNGO.tree")]

for(domain in c("BP","CC")){
    
  ave.loeuf=NULL
  for(i in 1:nrow(SYNGO.DF[[domain]])){
    id=rownames(SYNGO.DF[[domain]])[i]
    ensg=unlist(strsplit(syngo.ontol[id,"all.genes"],";"))
    ave.loeuf=c(ave.loeuf,mean(asc[ensg,"LOEUF"],na.rm=T))
  }
  plot(SYNGO.DF[[domain]]$Z.all,ave.loeuf)
  plot(SYNGO.DF[[domain]]$rateRatio.all,ave.loeuf,xlim=c(0,10))
  
  SYNGO.DF[[domain]]=cbind.data.frame(SYNGO.DF[[domain]][,1:31],data.frame(ave.LOEUF=ave.loeuf),SYNGO.DF[[domain]][,32:40])
}

# now collect the ancestry information
syngo.df[,c("n.gen","SYNGO.tree")]=NA
for(go.ids in rownames(syngo.df)){
  while(sum(is.na(go.ids)) == 0){
    go.ids=unique(c(go.ids,syngo.df[go.ids,"parent_id"]))
  }
  go.ids=go.ids[!is.na(go.ids)]
  syngo.df[go.ids[1],"n.gen"]=length(go.ids)
  syngo.df[go.ids[1],"SYNGO.tree"]=paste(go.ids,collapse="; ")
}

### split by CC and BP
SYNGO.DF=list()
SYNGO.DF[["CC"]]=syngo.df[syngo.df$domain == "CC",]
SYNGO.DF[["BP"]]=syngo.df[syngo.df$domain == "BP",]
SYNGO.DF[["CC"]]$orig.order=1:nrow(SYNGO.DF[["CC"]])
SYNGO.DF[["BP"]]$orig.order=1:nrow(SYNGO.DF[["BP"]])

head(SYNGO.DF$BP)


WriteXLS(SYNGO.DF,file.path(out_dir,"05_01_syngo-asc-summary.xlsx"))

for(domain in c("BP","CC")){
  
  ave.loeuf=NULL
  for(i in 1:nrow(SYNGO.DF[[domain]])){
    id=rownames(SYNGO.DF[[domain]])[i]
    ensg=unlist(strsplit(syngo.ontol[id,"all.genes"],";"))
    ave.loeuf=c(ave.loeuf,mean(asc[ensg,"LOEUF"],na.rm=T))
  }
  plot(SYNGO.DF[[domain]]$Z.all,ave.loeuf)
  plot(SYNGO.DF[[domain]]$rateRatio.all,ave.loeuf,xlim=c(0,10))
  
  SYNGO.DF[[domain]]=cbind.data.frame(SYNGO.DF[[domain]][,1:31],data.frame(ave.LOEUF=ave.loeuf),SYNGO.DF[[domain]][,32:40])
}

WriteXLS(SYNGO.DF,file.path(out_dir,paste0("05_01_syngo-asc-summary-incl-AVE-LOEUF-",Sys.Date(),".xlsx")))


pdf(file.path(supp_fig_dir,paste0("05_01_Z-all-vs-ave-loeuf-",Sys.Date(),".pdf")),height=8,width=8)
par(mfrow=c(2,1))
par(mar=c(5,5,1,1))
plot(SYNGO.DF[["BP"]]$Z.all,SYNGO.DF[["BP"]]$ave.LOEUF,xlab="Z.all",ylab="Average LOEUF",pch=20,col=viridis::viridis(9)[3],las=1,xlim=c(-5,15))
legend("topright",legend="BP",bty="n",cex=1.5)
plot(SYNGO.DF[["CC"]]$Z.all,SYNGO.DF[["CC"]]$ave.LOEUF,xlab="Z.all",ylab="Average LOEUF",pch=20,col=viridis::viridis(9)[6],las=1,xlim=c(-5,15))
legend("topright",legend="CC",bty="n",cex=1.5)
dev.off()

pdf(file.path(supp_fig_dir,"05_01_Z-score-vs-Entropy.pdf"),height=6,width=8)
par(mfrow=c(1,2))
plot(SYNGO.DF[["CC"]]$Z.all,SYNGO.DF[["CC"]]$HPro.all,pch=20,col="red",las=1,xlab="Z-score",ylab="Entropy",main="Cellular Components",ylim=c(0,6),xlim=c(0,14))
re=lm(SYNGO.DF[["CC"]]$HPro.all~SYNGO.DF[["CC"]]$Z.all)
abline(re$coef,lty=1,col="blue")
legend("bottomright",legend=c("a = 0.928","b = 0.386","R^2 = 0.647"))

plot(SYNGO.DF[["BP"]]$Z.all,SYNGO.DF[["BP"]]$HPro.all,pch=20,col="blue",las=1,xlab="Z-score",ylab="Entropy",main="Biological Processes",ylim=c(0,6),xlim=c(0,14))
re=lm(SYNGO.DF[["BP"]]$HPro.all~SYNGO.DF[["BP"]]$Z.all)
summary(re)
abline(re$coef,lty=1,col="red")
legend("bottomright",legend=c("a = 1.008","b = 0.324","R^2 = 0.348"))

dev.off()


SYNGO.DF[["CC"]][order(SYNGO.DF[["CC"]]$is.parent,SYNGO.DF[["CC"]]$Z.unq,decreasing=c(T,T))[1:10],]

