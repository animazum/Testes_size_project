
###
### PGLS de SSd con testis size 
### 31 july 2024
### We will going to use the data created in the SSD project already published and analyse it with Jeffs data using Testis size

### Prepariong for PGLS 

## install older version of dplyr
# devtools::install_url("http://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.4.tar.gz")

## for pre pgls analyses
library("logr")
library("dplyr")
library("ggpubr")

## for pgls
library("ape")
library("nlme")
library("parallel")

## for the plot 
library("reshape")
library("ggplot2")

# setwd("/Volumes/PortableSSD/bpm29/Orthologues_backup_130spp_mammalian/")
# setwd("/Users/mona/Dropbox/SSD/")
# setwd("/datoslab/BPM/SSD/")
# setwd("/Users/mona/Dropbox/TestisSize_SSD/Datos/")
# setwd("~/Dropbox/TestisSize_SSD/Datos/")

print("change always the ulimit of the computer $ulimit -s 21000 !!!!!!!!!!!!")
rm(list=ls())

path1<- "/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/Results/"
path1<- "/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/Results/"

SSDbias <- "MALE"
trait<- "longevity_brainmass"
Correcting<- "benjamini"

tmp <- file.path(getwd(), paste("GenefamilySize.",SSDbias,"longevity.",trait,".associated.genes.",Correcting,".correction.log", sep = ""))
lf <- log_open(tmp, traceback = F)

##### First we need to do the genome completion filtering. ####
## upload all the orthologs (gene families) per spp list. 
# Gene counts per orthogroup/Users/mona/Dropbox/SSD
# list_orthogroups_counts<- read.csv("~/Dropbox/kilili_shared_project/r", header = T, row.names = 1)
# list_orthogroups_counts<- read.table(paste(path1,"/Orthogroups.GeneCount.tsv", sep = ''), header = T, row.names = 1)
# log_print(paste("Species in the list:",dim(list_orthogroups_counts)[2]), hide_notes = T)
# colnames(list_orthogroups_counts)<- gsub(pattern = ".LS", replacement = "", x = colnames(list_orthogroups_counts))

# list_orthogroups_counts$Total<- NULL
### Reading the TESTIS data
# Nombre.file<-c(paste(path1,"/testes_mass_dataset.csv", sep = ""))
Nombre.file<-c(paste(path1,"/traits_albus_61spp.csv", sep = ""))
traits.filtered<-read.csv(Nombre.file, stringsAsFactors = F)


# tree <-read.tree("/Volumes/PortableSSD/Mammalian_orthologs_project/CDS_refseq_mammals_200816/Sequences_For_Orthofinder/OrthoFinder/Results_Apr25_3/Species_Tree/SpeciesTree_rooted_node_labels.txt")
tree <-read.tree(paste(path1,"/tree_albus_61spp.tree", sep = ""))
# tree$tip.label<-gsub(pattern = ".LS", replacement = "", x = tree$tip.label)
# 
# tree2<-drop.tip(tree,setdiff(tree$tip.label, traits$Spp.Name))

# list_orthogroups_counts<-list_orthogroups_counts[, colnames(list_orthogroups_counts) %in% tree2$tip.label]
# 
# traits<- traits[!is.na(traits$Combined.Testes.Mass..g.),]
# traits.rel.testes.size<-traits
# traits.rel.testes.size$Predicted.testes.size <- 0
# traits.rel.testes.size$Rel.testes.size <- 0


#*************   IT IS IMPORTANT TO HAVE THE SAME SPP ORDER BETWEEN TRAITS AND GENE NUMBERSS!!!!!!!   ********###
gene.numbers.filtered <- read.csv(paste(path1,"gene_numbers_albus_61spp.csv", sep = ""),header = T, row.names = 1)
dim(traits.filtered)
# total of gene counts per spp
sumgenes<-rowSums(t(gene.numbers.filtered))

# Note: from the 27801 gene families at the begining until this point, after the 20% of zeros, variance 
# equal to 0 and removing the orthogroups or gene families with only one gene among all the spp we ended up with
# 13753 orthogroups for 125 spp in total. 


PGLS.GF.size <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it, n_bootstrap = 100, percentTotalSample) {
  phenotypes <<- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  log_print("#########################", hide_notes = T)
  log_print("##### PGLS PGLS PGLS ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  
  # Initialize output storage for bootstrap results
  bootstrap_results <- vector("list", n_bootstrap)
  
  for (b in seq_len(n_bootstrap)) {
    log_print(paste("Bootstrap iteration:", b), hide_notes = T)
    
    # size of the data frame 
    Sizesample<-as.numeric(round(dim(traits.filtered)[1]*(percentTotalSample/100)))
    # Resample species with replacement
    resampled_species <- sample(traits[, spp.col.num], replace = F, size = Sizesample)
    # resampled_species <- sample(traits.filtered[, 1], replace = F, 60)
    
    # traits_resampled <- traits.filtered[traits.filtered[, 1] %in% resampled_species, , drop = FALSE]
    traits_resampled <- traits[traits[, spp.col.num] %in% resampled_species, , drop = FALSE]
    
    # Drop tips from the tree that are not in the resampled dataset
    tree_resampled <- drop.tip(tree, setdiff(tree$tip.label, traits_resampled[, 1]))
    
    # Drop species in gene numbers
    gene.numbers_resampled <- gene.numbers[,colnames(gene.numbers) %in% resampled_species , drop = FALSE]
    # gene.numbers_resampled <- gene.numbers.filtered[,colnames(gene.numbers.filtered) %in% resampled_species , drop = FALSE]
    
    
    if (No.variables == 1) {
      ###################################
      log_print("#### PGLS ONE Variable ######", hide_notes = T)
      ###################################
      out <- as.data.frame(t(apply(X = gene.numbers_resampled, 1, FUN = function(x) {
        traits2 <<- cbind(traits_resampled[, c(phenotypes)], as.numeric(x)) # add phenotypes as needed
        names(traits2)[dim(traits2)[2]] <<- "GFS"
        
        # Remove species with NA gene numbers
        tree_sub <- drop.tip(tree_resampled, names(x)[is.na(x)])
        
        # Create the formula for the model for 1 variable
        formilin <<- as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        formilin2 <- as.formula(paste("~", phenotypes[1], sep = ""))
        pglsModel <<- try(gls(model = formilin, correlation = corBrownian(phy = tree_sub, form = formilin2), method = "ML", data = traits2))
        
        if (inherits(pglsModel, "try-error"))
          return(rep(NA, 12)) # Adjusted to match output dimensions
        else
          return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2],
                   summary(pglsModel)$tTable[1, 1], summary(pglsModel)$tTable[2, 1], summary(pglsModel)$tTable[1, 2], summary(pglsModel)$tTable[2, 2],
                   summary(pglsModel)$tTable[1, 3], summary(pglsModel)$tTable[2, 3], summary(pglsModel)$tTable[1, 4], summary(pglsModel)$tTable[2, 4]))
      })))
      
      colnames(out) <<- c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2], sep = " "), paste("p_value", phenotypes[2], sep = " "),
                          "Coefficient Intercept", paste("Coefficient", phenotypes[2], sep = " "), "Coef SE Intercept", paste("Coef SE", phenotypes[2], sep = " "),
                          "Coef tval Intercept", paste("Coef tval", phenotypes[2], sep = " "), "Coef pval Intercept", paste("Coef pval", phenotypes[2], sep = " "))
    }
    
    # Save bootstrap results
    bootstrap_results[[b]] <- out
  }
  
  # Combine bootstrap results
  combined_results <- do.call(rbind, bootstrap_results)
  
  # Save combined results
  nameout <- paste("PGLS_Bootstrap", paste(phenotypes[2:length(phenotypes)], collapse = "_"), "Results.csv", sep = ".")
  log_print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  write.csv(combined_results, paste(where.Save.it, nameout, sep = "/"))
}

dirSave <- getwd()

PGLS.GF.size(traits = traits.filtered, pheno.col.nums = c(10), No.variables = 1, 
             tree = tree, spp.col.num = 1, gene.numbers = gene.numbers.filtered, 
             where.Save.it = dirSave, n_bootstrap = 10, percentTotalSample = 90)

