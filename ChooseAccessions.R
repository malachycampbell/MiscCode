#This function will try to find a set of accessions that maximizes variance for a set of loci.

#######################
##Preprocess SNP data##
#######################

#Acc.44k <- read.csv("~/Downloads/RiceDiversity.44K.germplasm.csv")
Acc700k <- read.delim("~/Downloads/HDRA-G6-4-RDP1-RDP2-NIAS/Misc/HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv", 
                      sep = "\t", header = F)
Acc700k <- Acc700k[grep("NSFTV", Acc700k$V3) ,]

FAM <- read.table("~/Downloads/HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.fam", sep = "", header = F)
FAM <- FAM[FAM$V1 %in% Acc700k$V1 ,]

write.table(FAM, "~/Downloads/Acc.txt", sep = "\t", quote = F, row.names = F)


##################################################
##PLINK command to extract phenotyped accessions##
##################################################
#plink --bfile HDRA-G6-4-RDP1-RDP2-NIAS --keep ~/Downloads/Acc.txt --maf 0.05 --make-bed --out RDP1



#######################
##PLINK to SNP matrix##
#######################
library(BGLR)
FAM <- read.table("~/Downloads/HDRA-G6-4-RDP1-RDP2-NIAS/RDP1.fam")[1:2]
Acc700k <- Acc700k[match(FAM$V1, Acc700k$V1) ,]
MAP <- read.table("~/Downloads/HDRA-G6-4-RDP1-RDP2-NIAS/RDP1.bim")
PED <- read_bed("/Users/malachycampbell/Downloads/HDRA-G6-4-RDP1-RDP2-NIAS/RDP1.bed", 
                "/Users/malachycampbell/Downloads/HDRA-G6-4-RDP1-RDP2-NIAS/RDP1.fam", 
                "/Users/malachycampbell/Downloads/HDRA-G6-4-RDP1-RDP2-NIAS/RDP1.bim")
m <- PED$p
n <- PED$n
PED <- PED$x
##SNPs in PED are coded as 0, 1, 2, 3. 2 is missing data. 1 are heterozygous, 0 and 3 are homozygous for 1/1 and 2/2 for major allele and minor allele respectively
PED[PED == 2] <- NA 
PED[PED == 0] <- 0
PED[PED == 1] <- 1
PED[PED == 3] <- 2

W <- t(matrix(PED, nrow=m, ncol=n, byrow = T))
rownames(W) <- MAP$V2
colnames(W) <- sub("NSFTV", "NSFTV_", Acc700k$V3)


################################################
#####Function to choose a set of accessions#####
##that maximizes the variance at a set of loci##
################################################


SNPofInt <- t(read.csv("~/Downloads/significant_snps.csv", row.names = 1))
BestSet <- function(Markers = NULL, Popsize = NULL, NoReps = NULL, SNPsOfInterest = NULL, SEED = NULL){
  #Places to store results
  AccList <- list()
  FreqList <- list()
  VarRes <- NULL
  
  #Subset markers for SNPs you want.
  W <- Markers[row.names(Markers) %in% SNPsOfInterest ,]
  for( i in 1:NoReps){
    #Sample for X number of Reps
    set.seed(SEED + i)
    Indx <- sample(1:ncol(W), size = Popsize)
    
    #Subset the marker matrix
    tmpW <- W[, Indx]
    
    #Get variance of allele counts
    tmp.var <- apply(tmpW, 1, var, na.rm = T)
    #If any loci have no variance then set the sum as 0. This ensures that the best set maximizes variance at all loci
    if(sum(tmp.var == 0) > 0){
      tmp.var = 0
    }else{
      tmp.var = sum(tmp.var)
    }
    #Compile a list of results
    tmpAccs <- colnames(W)[Indx]
    AccList[[i]] <- tmpAccs
    FreqList[[i]] <- rowMeans(tmpW)/2
    VarRes <- c(VarRes, tmp.var)
  }
  
  #Select the best
  BestSets <- which(VarRes == max(VarRes))
  print(BestSets)
  BestVar <- VarRes[BestSets]
  BestAcc <- AccList[[BestSets]]
  BestFreq <- FreqList[[BestSets]]
  
  return(list(Variance = BestVar,
              BestAccesssions = BestAcc,
              AlleleFreq = BestFreq))
}

BestSet(Markers = SNPofInt, Popsize = 15, NoReps = 100000, SNPsOfInterest = row.names(SNPofInt), SEED = 1)
