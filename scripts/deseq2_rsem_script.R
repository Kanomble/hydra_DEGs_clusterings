library(tximport)
library(DESeq2)
library(DEGreport)
library(limma)
library(ggplot2)

# define input file paths
sample_file <- file.path("data/rsem_counts/sample_descriptions.csv") # "data/rsem_counts/sample_descriptions_without_temp_oxo.csv"
controls_file <- file.path("data/rsem_counts/sample_controls.csv") # "data/rsem_counts/sample_controls_without_temp_oxo.csv"

input_dir <- "data/rsem_counts/"
output_dir <- "results/deseq2_rsem/"

# settings
count_threshold <- 20

# read samples file
samples <- read.csv(sample_file)
controls <- read.csv(controls_file)

# extract unique experiments for looping
experiments <- unique(samples$EXP)
dds_objects <- list()

for(exp in experiments){
  # get experimental RSEM count file paths
  expSamples <- samples[samples$EXP == exp,]
  files <- file.path(input_dir, exp, paste0(expSamples$ID,".genes.results"))
  names(files) <- expSamples$ID
  # check if all files exist
  all(file.exists(files)) # files[!file.exists(files)]
  # recommended tximport vignette cmd
  txi.rsem <- tximport(files, type="rsem", txIn = FALSE, txOut = FALSE)
  # create sampleTable for deseq2
  sampleTable <- data.frame(condition = factor(expSamples$treatment), batch = factor(expSamples$batch))
  rownames(sampleTable) <- colnames(txi.rsem$counts)
  # deseq2 workflow
  dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
  
  if(exp=="HydraRecolonization"){
    cat("[*] Removing bad samples from HydraRecolonization data\n")
    bad_samples <- c("I27673-S1","I27674-S1","I27675-S1","I27678-S1","I27681-S1","I27672-S1")
    dds <- dds[, !colnames(dds) %in% bad_samples]
  } 
   
  keep <- rowSums(counts(dds)) >= count_threshold
  dds <- dds[keep,]
  
  if(exp == "HydraAHL"){
    cat("[*] Releveling To Real Reference\n")
    dds$condition <- relevel(dds$condition, ref = "control")
  } 
  
  if(exp == "EcoKD1"){
    cat("[*] Releveling To Real Reference\n")
    dds$condition <- relevel(dds$condition, ref = "control_B8")
    dds_corrected <- dds
    vsd <- vst(dds_corrected)
    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
    png(file=file.path(output_dir, "figures",paste0(exp,"_batch_corrected_pca.png")), width=800, height=550)
    pca <- plotPCA(vsd, intgroup="condition")
    print(pca)
    dev.off()
  }
  
  # design(dds) <- ~batch + condition
  dds <- DESeq(dds)
  
  dds_objects[[length(dds_objects) + 1]] <- dds
  
  # PCA
  rld <- rlog(dds, blind=TRUE)
  cat(paste("[*] PRODUCING PCA",file.path(output_dir,"figures",paste0(exp,"_pca.png")),"\n"))
  png(file=file.path(output_dir, "figures",paste0(exp,"_pca.png")), width=800, height=550)
  pca <- plotPCA(rld, intgroup="condition")
  print(pca)
  dev.off()
  
  if(exp == "HydraRecolonization"){
    dds_corrected <- dds
    vsd <- vst(dds_corrected)
    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
    png(file=file.path(output_dir, "figures",paste0(exp,"_batch_corrected_pca.png")), width=800, height=550)
    pca <- plotPCA(vsd, intgroup="condition")
    print(pca)
    dev.off()
  }
  
  cat("[*] DONE\n")
  
  cntrl <- controls[controls$EXP==exp,-1]
  ## obtaining results
  for ( i in 1:nrow(cntrl) ) {
    res <- results(dds, contrast=c("condition",
                                   cntrl$treatment.group[i],
                                   cntrl$control.group[i]))
    resOrdered <- res[order(res$pvalue),]
    
    ## writing output file
    write.csv(as.data.frame(resOrdered), 
              file=file.path(output_dir,"tables",paste0(exp,"_",
                                                        cntrl$treatment.group[i],"_vs_",
                                                        cntrl$control.group[i],
                                                        "_results.csv")))
    
    
    ## volcano plot
    
    cat(paste("[*] PRODUCING VOLCANO PLOT",file.path(output_dir,"figures",paste0(exp,"_",
                                                                                 cntrl$treatment.group[i],"_vs_",
                                                                                 cntrl$control.group[i],
                                                                                 "_volcano_plot.png")),"\n"))
    resOrdered[["id"]] <- row.names(resOrdered)
    show <- as.data.frame(resOrdered[1:10, c("log2FoldChange", "padj","id")])#"id"
    png(file=file.path(output_dir,"figures",paste0(exp,"_",
                                                   cntrl$treatment.group[i],"_vs_",
                                                   cntrl$control.group[i],
                                                   "_volcano_plot.png")),width=800, height=550)
    print(degVolcano(resOrdered[,c("log2FoldChange", "padj")],plot_text = show))
    cat("[*] DONE\n")
    dev.off()
  }
}