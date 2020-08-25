library(DESeq2)
library(biomaRt)
library(PharmacoGx)

mtx <- read.table("mRNA_counts.txt", header=T, row.names = "gene_symbols")
samples <- c("Clone1_N1", "Clone1_N2", "Clone1_N3",
             "Clone9_N1", "Clone9_N2", "Clone9_N3",
             "Control_N1", "Control_N2", "Control_N3")
colnames(mtx) <- samples
sample_vec <- c()
replicates_vec <- c()
sample_df <- as.data.frame(matrix(ncol=2, nrow=0))
for(i in samples){
  split <- strsplit(i, "\\_")
  sample <- split[[1]][1]
  rep <- split[[1]][2]
 sample_df[i,1] <- sample
  sample_df[i,2] <- rep
}
colnames(sample_df) <- c("Condition", "Replicate")
sample_df$Condition <- as.factor(sample_df$Condition)
sample_df$Replicate <- as.factor(sample_df$Replicate)
dds <- DESeqDataSetFromMatrix(mtx, colData = sample_df, design = ~ Replicate + Condition)
dds$Condition <- relevel(dds$Condition, ref = "Control")
  keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
dds <- DESeq(dds) 
cts <- counts(dds, normalized=T)

CMAP.sigs <- downloadPertSig("CMAP")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Function 1: Take Up-regulated | Down-reguated results table, filter by LFC +/- 2 and return ranked gene list. 
step_1 <- function(up, down){
  up_filt <- subset(up, Log2FC >= 2)
  down_filt <- subset(down, Log2FC <= -2)
  
  top <- as.data.frame(up_filt$Gene)
  top$direction <- up_filt$Log2FC
  top$padj <- up_filt$Adj.P.value
  colnames(top) <- c("feature", "direction", "padj")
  
  bottom <- as.data.frame(down_filt$Gene)
  bottom$direction <- down_filt$Log2FC
  bottom$padj <- down_filt$Adj.P.value
  colnames(bottom) <- c("feature", "direction", "padj")

  tmp <- rbind(top,bottom)
  tmp <- tmp[order(tmp$padj),]
  
  return(tmp)
}

# Function 2: Retrieve normalised expression values for Genes, calculate mean expression and filter out lower 50% quantile of mean expression. **PROVIDE SUBSETTED COUNTS MATRIX**. Retrieve Ensembl gene ID and append "_at" to string to generate final matrix for PharmacoGx.
step_2 <- function(sample_counts, de_genes){
    df <- as.data.frame(sample_counts[which(rownames(sample_counts) %in% de_genes$feature),])
    df$mean_exp <- rowMeans(df)
    
    quant <- quantile(df$mean_exp)
    quant_50 <- as.numeric(quant[3])
    
    df <- df[which(df$mean_exp > quant_50),]
    de_genes <- de_genes[which(de_genes$feature %in% rownames(df)),]
    rownames(de_genes) <- de_genes$feature
    
    ens <- getBM(attributes=c("ensembl_gene_id",
                             "external_gene_name"),
                        filters = c("external_gene_name"),
                        values = de_genes$feature,
                        mart = mart)
    # CMAP uses very strange probe names. ENSG000001_at
    ens$ensembl_gene_id <- lapply(ens$ensembl_gene_id, function(x) paste(x,"_at"))
    ens$ensembl_gene_id <- gsub(" ", "", ens$ensembl_gene_id)
    ens$ensembl_gene_id <- gsub(" ", "", ens$ensembl_gene_id)
    colnames(ens) <- c("ensembl", "feature")
    de_genes <- merge(de_genes, ens, by = "feature")
    rownames(de_genes) <- de_genes$ensembl
    de_genes <- de_genes[,c(4,2,3)]
    colnames(de_genes) <- c("feature", "direction", "padj")
    de_genes <- de_genes[order(de_genes$padj),]  
    return(de_genes)
}

# Fucntion 3: Extracting top | bottom hits for permutations:
step_3 <- function(df, x){
  tmp <- as.data.frame(head(df, n = x))
  return(tmp)
}


clone1_up <- read.table("Clone1_vs_Control_UP.txt", header =T, stringsAsFactors = F)

clone1_down <- read.table("Clone1_vs_Control_DOWN.txt", header =T, stringsAsFactors = F)
# run function 1:
clone1_v_control <- step_1(clone1_up, clone1_down)

# get samples from counts by indexing
key1 <- cts[,c(1,2,3,7,8,9)]

# run function 2:
clone1_v_control_PGX <- step_2(key1, clone1_v_control)

input <- step_3(clone1_v_control_PGX, 300)

  res <- apply(CMAP.sigs[,,c("tstat", "fdr")],
                2, function(x, input){
                            return(connectivityScore(x=x,
                                                    y=input[,2,drop=FALSE],
                                                    method="gsea", nperm=1000,
                                                    nthread = 8))
                                                    }, input=input)

  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  res <- as.data.frame(res[order(res[,1], decreasing=TRUE),])

write.table(res, "Clone1_vs_Control_n300.txt")




# Set up permuations

#grid <- seq(100,1250, by =50)
#rows <- length(grid)

#tracker <- matrix(ncol = 3, nrow = rows)
#colnames(tracker) <- c("iteration", "Connectivity_max", "n_sig")
#count = 1
#for(i in grid){
#  input <- step_3(clone1_v_control_PGX, i)
  
#  res <- apply(CMAP.sigs[,,c("tstat", "fdr")],
#                2, function(x, input){
#                            return(connectivityScore(x=x,
#                                                    y=input[,2,drop=FALSE],
#                                                    method="gsea", nperm=1000,
#                                                    nthread = 8))
#                                                    }, input=input)

#  rownames(res) <- c("Connectivity", "P Value")
#  res <- t(res)
#  res <- as.data.frame(res[order(res[,1], decreasing=TRUE),])

#  n_sig <- nrow(res[which(res$`P Value` < 0.05),])
#  x <- res[which(res$`P Value` < 0.05),]
#  conn_min <- min(x$Connectivity)
  
#  tracker[count,1] <- i
#  tracker[count,2] <- conn_min
#  tracker[count,3] <- n_sig

 # count <- count + 1
#}

#write.table(tracker, "Clone1_vs_Control_premutations_min.txt")
  
#clone1 vs clone 9 
clone19_up <- read.table("Clone1_vs_Clone9_UP.txt", header =T, stringsAsFactors = F)

clone19_down <- read.table("Clone1_vs_Clone9_DOWN.txt", header =T, stringsAsFactors = F)
# run function 1:
clone1_v_clone9 <- step_1(clone19_up, clone19_down)

# get samples from counts by indexing
key19 <- cts[,c(1,2,3,4,5,6)]

# run function 2:
clone1_v_clone9_PGX <- step_2(key19, clone1_v_clone9)

input <- step_3(clone1_v_clone9_PGX, 100)

  res <- apply(CMAP.sigs[,,c("tstat", "fdr")],
                2, function(x, input){
                            return(connectivityScore(x=x,
                                                    y=input[,2,drop=FALSE],
                                                    method="gsea", nperm=1000,
                                                    nthread = 8))
                                                    }, input=input)

  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  res <- as.data.frame(res[order(res[,1], decreasing=TRUE),])

write.table(res, "Clone1_vs_Clone9_n100.txt")




# Set up permuations

#grid <- seq(100,1150, by =50)
#rows <- length(grid)

#tracker <- matrix(ncol = 3, nrow = rows)
#colnames(tracker) <- c("iteration", "Connectivity_max", "n_sig")
#count = 1
#for(i in grid){
#  input <- step_3(clone1_v_clone9_PGX, i)
  
#  res <- apply(CMAP.sigs[,,c("tstat", "fdr")],
#                2, function(x, input){
#                            return(connectivityScore(x=x,
#                                                    y=input[,2,drop=FALSE],
#                                                    method="gsea", nperm=1000,
#                                                    nthread = 8))
#                                                    }, input=input)

#  rownames(res) <- c("Connectivity", "P Value")
#  res <- t(res)
#  res <- as.data.frame(res[order(res[,1], decreasing=TRUE),])
#
#  n_sig <- nrow(res[which(res$`P Value` < 0.05),])
#  x <- res[which(res$`P Value` < 0.05),]
#  conn_min <- min(x$Connectivity)
#  
#  tracker[count,1] <- i
#  tracker[count,2] <- conn_min
#  tracker[count,3] <- n_sig

#  count <- count + 1
#}

#write.table(tracker, "Clone1_vs_Clone9_permutations_min.txt")

clone9_up <- read.table("Clone9_vs_Control_UP.txt", header =T, stringsAsFactors = F)

clone9_down <- read.table("Clone9_vs_Control_DOWN.txt", header =T, stringsAsFactors = F)
# run function 1:
clone9_v_control <- step_1(clone9_up, clone9_down)

# get samples from counts by indexing
key9 <- cts[,c(4,5,6,7,8,9)]

# run function 2:
clone9_v_control_PGX <- step_2(key9, clone9_v_control)


input <- step_3(clone9_v_control_PGX, i)

  res <- apply(CMAP.sigs[,,c("tstat", "fdr")],
                2, function(x, input){
                            return(connectivityScore(x=x,
                                                    y=input[,2,drop=FALSE],
                                                    method="gsea", nperm=1000,
                                                    nthread = 8))
                                                    }, input=input)

  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  res <- as.data.frame(res[order(res[,1], decreasing=TRUE),])

# Set up permuations
write.table(res, "Clone9_vs_Control_n110")

#grid <- seq(100,144, by =10)
#rows <- length(grid)

#tracker <- matrix(ncol = 3, nrow = rows)
#colnames(tracker) <- c("iteration", "Connectivity_max", "n_sig")
#count = 1
#for(i in grid){
#  input <- step_3(clone9_v_control_PGX, i)
#
#  res <- apply(CMAP.sigs[,,c("tstat", "fdr")],
#                2, function(x, input){
#                            return(connectivityScore(x=x,
#                                                    y=input[,2,drop=FALSE],
#                                                    method="gsea", nperm=1000,
#                                                    nthread = 8))
#                                                    }, input=input)
#
#  rownames(res) <- c("Connectivity", "P Value")
#  res <- t(res)
#  res <- as.data.frame(res[order(res[,1], decreasing=TRUE),])
#
#  n_sig <- nrow(res[which(res$`P Value` < 0.05),])
#  x <- res[which(res$`P Value` < 0.05),]
#  conn_min <- min(x$Connectivity)
#  
#  tracker[count,1] <- i
#  tracker[count,2] <- conn_min
#  tracker[count,3] <- n_sig
#
#  count <- count + 1
#}
#
#write.table(tracker, "Clone9_vs_Control_permutation_min.txt")
