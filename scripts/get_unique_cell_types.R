### ANALYZE DESEQ2 RESULTS ###

options(stringsAsFactors=FALSE)

### PATHS ###
## input paths
resp <- file.path("results/deseq2_rsem/")
datp <- file.path("results/processedData/") 

expf <- file.path("data/rsem_counts/sample_descriptions.csv")
rawf <- file.path("results/deseq2_rsem/hydra_all_counts_t_to_g.tsv")
datf <- file.path("results/deseq2_rsem/hydra_all_deseq2_lg2_t_to_g.tsv")
pvlf <- file.path("results/deseq2_rsem/hydra_all_deseq2_pvalues_t_to_g.tsv")

## annotation
cltf <- file.path("./results/processedData/hvaep_clustering_table.tsv")

### ANALYSIS PARAMETERS ###

## data set filter: remove some redundant comparison
## data order: 
selected.data <- c("Wild_vs_GF","Cvbct_vs_GF","Conventionalized_vs_GF",
                   "3OHC12_vs_control")# ,"Eco1KD_B8_vs_control_B8"
#
## SET THIS to FALSE if you use different combinations of controls!
## TODO: use this also for other plots!
short.class.names <- TRUE 

## data filter: cells with cell-type annotations, etc?
filter <- "gfzero" #celltype" "outlier" # 
filt.nap <- TRUE

## celltypes: only analyze unique gene-cell type combis
unique.celltype <- TRUE

## filter for up/down classes
lg2.cutoff <- 1
p.cutoff <- 0.05

## volcano plot cutoffs
maxx <- 4 # 8
maxy <- 10


## parse data
dat <- read.delim(datf)
pvl <- read.delim(pvlf)
clt <- read.delim(cltf) # sieber et al. cell types
got <- read.delim(gof)
#ann <- read.delim(anf)

# set all values lg2fold changes of dat < -15 and > 15 to NA
#identifier <- dat$ID
#dat[, -1][abs(dat[, -1]) >= 15] <- NA
#dat$ID <- identifier
#pvl[is.na(dat)] <- NA



## clean names: remove experiment names from column names
colnames(dat) <- sub(".*?_","",colnames(dat))
#colnames(pvl) <- sub(".*?_","",colnames(pvl))

## clean IDs in Siebert et al. supplement
#ann[,1] <- sub("\\|.*","", ann[,1])
#colnames(ann)[1] <- "ID"

## check identity
if ( any(dat$ID!=pvl$ID) ) stop("wrong IDs in p-value file")

## expand Annotation data to RNAseq data file
clte <- clt[match(dat$ID, clt$ID),]

## custom filter 
RM <- logical()
if ( filter=="ct_or_go" ) {
  ## rm: not present in any annotation file
  RM <- is.na(clte$ID) & is.na(gote$ID)
} else if ( filter=="celltype" ) {
  ## rm: not present in cell type annotation
  RM <- is.na(clte$ID)
} else if ( filter=="outlier" ) {
  ## rm: present in cell type annotation
  RM <- !is.na(clte$ID) 
  maxx <- 8
} else if ( filter=="gfzero" ) {
  
  ## get experiment sample IDs
  experiments <- read.csv(expf)
  gfid <- experiments$ID[experiments$treatment=="GF"]
  ## get raw data
  raw <- read.delim(rawf)
  raw <- raw[match(dat$ID, raw$gene_id),]
  ## transcript lengths
  taep <-raw[,1:3]
  ## germ free: GF experiments
  raw <- raw[,paste0(sub("-S1",".S1",gfid),"_expected_count")]
  RM <- apply(raw,1,sum)==0
}


## SELECTED DATA, incl. order 
datf <- dat[,selected.data]
colnames(pvl) <- colnames(dat)
pvlf <- pvl[,selected.data]
ids <- dat$ID

if ( any(RM) ) { #18638
  
  cat(paste("removing", sum(RM) , "features with filter:\t", filter,"\n",
            "KEEPING", sum(!RM), "features for cluster analysis\n\n"))
  
  ## filter all data
  ids <-   ids[!RM]
  datf <- datf[!RM,]
  pvlf <- pvlf[!RM,]
  clt <-  clte[!RM,]
}


if ( short.class.names ) {
  warning("using short class names; this can cause errors ",
          "when combining different controls")
  colnames(datf) <- sub("_vs_.*","", colnames(datf))
  colnames(pvlf) <- sub("_vs_.*","", colnames(pvlf))
}

## FURTHER FILTERS

## remove all with p.value
if ( filt.nap ) {
  
  rmp <- which(apply(is.na(pvlf),1,sum)>=1) #ncol(pvlf)-1 --> removing all columns with p_value = NA
  
  cat(paste("removing", length(rmp) , "features with NA p-values\n",
            "KEEPING", nrow(pvlf)-length(rmp),
            "features for cluster analysis\n\n\n\n"))
  
  ids <-   ids[-rmp]
  datf <- datf[-rmp,]
  pvlf <- pvlf[-rmp,]
  clt <-   clt[-rmp,]
}

### COMPARISON TO CELL TYPES ###

## FILTER CELL TYPE ANNOTATION:

## super categories present?
## TODO: do this smarter, e.g. store super clustering in separate table
## or use prefix super:?

## tag cells w/o annotation - TODO:after removing non-unique??
cna <- apply(clt[,2:ncol(clt)],1,function(x) !any(!is.na(x)))

## remove annotations that are NOT unique to a cell type, per hierachy level
cltu <- clt

## super clusters present?
scidx <- min(grep("^L_", colnames(clt)))
pridx <- min(grep("^C_", colnames(clt)))

## column indexes for separate categories
upos <- sort(unique(c(2,scidx,pridx,1+ncol(clt))))
urng <- list()
for ( j in seq_along(upos[-1]) ) 
  urng[[j]] <- upos[j]:(upos[j+1]-1)

# celltype "uniqueness" filter
for ( j in seq_along(urng) ) {
  ## NOT unique to cell type
  rmc <- apply(clt[,urng[[j]]],1,sum,na.rm=TRUE) > 2
  cltu[rmc,urng[[j]]] <- NA
}

if ( unique.celltype ) clt <- cltu

## test C_neuron cluster 
sum(cltu$C_neuron,na.rm=TRUE)


outfile_uc <- paste0(datp,"unique_celltypes.table")
write.table(x=cltu,file=outfile_uc,quote = FALSE, sep = ';')