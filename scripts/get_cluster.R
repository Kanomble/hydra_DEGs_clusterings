library(readxl)
options(stringsAsFactors=FALSE)

# INPUT
cls.file <- file.path("./results/processedData/hvaep_cell_type_to_gene_cluster_table.tsf")

## OUTPUT
out.path <- file.path("./results/processedData/")
fig.path <- file.path("./results/figures/")

### PARAMETERS
## add super classes
add.super <- TRUE # FALSE #  

## read supplemental file
clsf <- read.csv(cls.file, sep="\t")
## investigate
table(clsf[,1])

clsf <- data.frame(clsf)

## get transcript IDs
tids <- sub("\\|.*","",clsf[,2])
sum(duplicated(tids))

## export
res <- cbind.data.frame(ID=tids, CLS=clsf[,1])
out.file <- file.path(out.path, "hvaep_clustering_normalized.tsv")
write.table(res, file=out.file, sep="\t",row.names=FALSE, quote=FALSE)


## analyze overlaps
tids <- unique(res$ID)
cids <- unique(res$CLS)

## generate T/F table
ctab <- matrix(FALSE, nrow=length(tids), ncol=length(cids))
rownames(ctab) <- tids
colnames(ctab) <- cids
for ( j in 1:ncol(ctab) ) 
  ctab[,j] <- tids%in%res$ID[res$CLS==cids[j]]

## add super-categories
if ( add.super ) {
  super <- list()
  
  super[["C_neuron"]] <- grep("N$",cids, value=TRUE) # 8 
  super[["C_germline"]] <- grep("GC",cids, value=TRUE) # 2
  super[["C_nematoblast"]] <- grep("NB$",cids, value=TRUE) # 3
  super[["C_nematocytes"]] <- c(grep("NC$", cids, value=TRUE), grep("I_EarlyNem",cids, value=TRUE)) # --> this is not in the superclass: I_EarlyNem --> it has to be included later
  super[["C_gland_cell"]] <- c(grep("GL$",cids, value=TRUE),grep("I_SpumMucGl",cids, value=TRUE),grep("I_GranGl",cids, value=TRUE),grep("I_ZymoGl",cids, value=TRUE)) # I_GranGl I_ZymoGl I_SpumMucGl | excluded for now: I_GlProgen 
  
  super[["C_progenitor"]] <- grep("I_GlProgen",cids, value=TRUE)
  super[["C_neuroprogenitor"]] <- grep("I_Neuro", cids, value=TRUE)
  super[["C_EC_foot"]] <- c(grep("Ec_BasalDisk", cids, value=TRUE), grep("Ec_Peduncle", cids, value=TRUE))
  super[["C_EN_foot"]] <- grep("En_Foot", cids, value=TRUE)
  super[["C_tentacle"]] <- grep("_Tentacle$", cids, value=TRUE)
  super[["L_foot"]] <- c(grep("Ec_BasalDisk", cids, value=TRUE), 
                         grep("Ec_Peduncle", cids, value=TRUE), grep("En_Foot", cids, value=TRUE))
  super[["C_head"]] <- grep("Head$", cids, value=TRUE)
  
  ## stem cell lineages # first level of hierachy
  super[["L_ENDODERM"]] <- grep("En_BodyCol/SC", cids, value=TRUE)
  super[["L_ECTODERM"]] <- grep("Ec_BodyCol/SC", cids, value=TRUE)
  super[["L_INTERSTITIAL"]] <- grep("I_ISC",cids, value=TRUE)
  
  # added drop = FALSE if there is just one column in x
  supm <- do.call(cbind,lapply(super, function(x) apply(ctab[, x, drop = FALSE], 1, any)))
  ctab <- cbind(ctab,supm)
}

out.file <- file.path(out.path,"hvaep_clustering_table.tsv")
clt <- apply(ctab,2,as.numeric)
clt[clt==0] <- NA
write.table(cbind.data.frame(ID=rownames(ctab),clt),
            file=out.file, sep="\t",row.names=FALSE, quote=FALSE, na="")
