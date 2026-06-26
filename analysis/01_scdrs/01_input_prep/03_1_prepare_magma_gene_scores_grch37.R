library(data.table)
library(dplyr)
annot = fread(file = "/mnt/isilon/gandal_lab/liaoyd/tools/magma/NCBI37.3.gene.loc",header = F)
colnames(annot)  = c("entrezID", "chromosome", "start_position", "end_position", "strand", "gene")
out_dir <- Sys.getenv("out_dir")
trait   <- Sys.getenv("trait")

result = fread(paste0(out_dir,"/step2/",trait,".genes.out"))
result_rename = left_join(result %>% rename(entrezID=GENE),
                          annot)
fwrite(result_rename, file = paste0(out_dir,"/step2/",trait,".genes.out.withannot"),sep = "\t")

result_scDRS = result_rename %>% dplyr::select(gene,ZSTAT) %>% rename(GENE=gene,!!trait := ZSTAT)
fwrite(result_scDRS, file = paste0(out_dir,"/step2/",trait,".genes.tsv"),sep = "\t")

result_scDRS = result_rename %>% dplyr::select(gene,P) %>% rename(GENE=gene,!!trait := P)
fwrite(result_scDRS, file = paste0(out_dir,"/step2/",trait,".genes_pval.tsv"),sep = "\t")

