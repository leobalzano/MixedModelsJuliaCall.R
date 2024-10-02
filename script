################################################
########      MixedModelsJuliaCall.R     #######
################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)
# Created: June/13/2024
# Last update: June/13/2024

# This script is to show the procedure followed to generate the mixed-effects models published in the article called:
# "Unique Lymphocyte Transcriptomic Profiles in Septic Patients with Chronic Critical Illness".

################################################
# Libraries:
library(Seurat)
library(JuliaCall)
julia_setup(rebuild = TRUE)
find.package("JuliaCall")

# Install julia programming language
# https://julialang.org/downloads/ 
# In a terminal curl -fsSL https://install.julialang.org | sh


julia_setup(JULIA_HOME = "/Users/leobalzano/.juliaup/bin") # where it was installed for me locally. or on an HPC
julia_command('import Pkg; Pkg.add("Distributions")')
julia_command('import Pkg; Pkg.add("MixedModels")')
julia_command('import Pkg; Pkg.add("DataFrames")')
julia_command('using Pkg; Pkg.activate("~/Desktop/");') ## put this on hipergator folder..
julia_command("using Distributions, MixedModels, DataFrames;")

################################################
# Data:
setwd ("<path/to/your/desired/folder/>")
s <- readRDS ("<path/to/your/desired/RDS.file>") 

# Subset your data to the cell types you need to study
seu_mono_sub_comp1 <- subset(s, predicted.celltype.l1_Sample_Status %in% c("CD8 T_CCI","CD8 T_RAP"))

# Establish the reference and the comparison groups
reference_group = "CCI"
comparison_group = "RAP"


## gene filtering
allgenes <- rownames(seu_mono_sub_comp1)
allm <- apply(seu_mono_sub_comp1@assays$RNA@data, 1, mean)
allgenes <- names(which(allm > 0))
length(allgenes)

res_df <- list()
for(j in 1:length(allgenes)) {
  print(j)
  useG <- allgenes[j]
  gene_data <- data.frame(cell = colnames(seu_mono_sub_comp1), 
                          exp = unname(seu_mono_sub_comp1@assays$RNA@data[useG, ]), 
                          sample_label = as.character(seu_mono_sub_comp1$Sample_Status), 
                          sample = seu_mono_sub_comp1$sample, 
                          celltype = as.character(seu_mono_sub_comp1$predicted.celltype.l1_Sample_Status))
  
  julia_assign("gene_data", gene_data)
  julia_assign("ref_grp", reference_group)
  julia_command("model_fit = fit(MixedModel, @formula(exp ~ 1 + sample_label + (1 | sample)), 
                gene_data, contrasts= Dict(:sample_label => DummyCoding(base = ref_grp)));")
  model_coef <- julia_eval("DataFrame(coeftable(model_fit));")
  eff_df <- julia_eval("sum(leverage(model_fit));")
  ll <- julia_eval("loglikelihood(model_fit);")
  res_df[[j]] <- data.frame(
    gene=allgenes[j],
    beta0 = model_coef$Coef.[1], 
    beta1 = model_coef$Coef.[2], 
    se_beta0 = model_coef$`Std. Error`[1], 
    se_beta1 = model_coef$`Std. Error`[2], 
    zstat_beta0 = model_coef$z[1], 
    zstat_beta1 = model_coef$z[2], 
    pval_beta0 = model_coef$`Pr(>|z|)`[1],  
    pval_beta1 = model_coef$`Pr(>|z|)`[2], 
    effective_df = eff_df, 
    LL = ll, 
    note = NA_character_)
  
}
all_res <- do.call(rbind, res_df)


all_res$adj.pval <- p.adjust(all_res$pval_beta1, method = "fdr")

getwd()
#save(all_res, file = "all_res_CD8T_CCIvsCD8_RAP_JuliaMixedModel.RData")
# load (file = "all_res_CD8T_CCIvsCD8_RAP_JuliaMixedModel.RData") # adjpval <0
sum(all_res$adj.pval < 0.1)
checkSub <- subset(all_res, adj.pval < 0.1)
dim(checkSub)

head(checkSub)
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#DE
genes <- checkSub$gene
#df<-df[,-4]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
head(G_list)
DE<- merge (checkSub,G_list, by.x="gene", by.y="ensembl_gene_id")
DE<-DE[order(DE$adj.pval),]
DE
head(DE)

# write.table(DE, file = "DE_0p1_CD8T_CCIvsCD8_RAP_JuliaMixedModel.txt")


# END
