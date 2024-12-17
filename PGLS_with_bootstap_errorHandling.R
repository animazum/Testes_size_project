
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

print("change always the ulimit of the computer $ulimit -s 21000 !!!!!!!!!!!!")
rm(list=ls())

path1<- "/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/Results/"

SSDbias <- "MALE"
trait<- "longevity_brainmass"
Correcting<- "benjamini"

tmp <- file.path(getwd(), paste("GenefamilySize.",SSDbias,"longevity.",trait,".associated.genes.",Correcting,".correction.log", sep = ""))
lf <- log_open(tmp, traceback = F)

Nombre.file<-c(paste(path1,"/traits_albus_61spp.csv", sep = ""))
traits.filtered<-read.csv(Nombre.file, stringsAsFactors = F)


tree <-read.tree(paste(path1,"/tree_albus_61spp.tree", sep = ""))

#*************   IT IS IMPORTANT TO HAVE THE SAME SPP ORDER BETWEEN TRAITS AND GENE NUMBERSS!!!!!!!   ********###
gene.numbers.filtered <- read.csv(paste(path1,"gene_numbers_albus_61spp.csv", sep = ""),header = T, row.names = 1)
dim(traits.filtered)
# total of gene counts per spp
sumgenes<-rowSums(t(gene.numbers.filtered))



PGLS.GF.size <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it, 
                         bootstrap = FALSE, n_bootstrap = 100, Sizesample = NULL) {
  # Initialize logging
  log_file <- file.path(where.Save.it, "pgls_analysis.log")
  if (!file.exists(log_file)) {
    file.create(log_file)
  }
  
  # Custom logging function that doesn't rely on futile.logger
  log_message <- function(msg, hide_notes = TRUE) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message <- paste(timestamp, msg, sep = " - ")
    cat(message, "\n", file = log_file, append = TRUE)
    if (!hide_notes) cat(message, "\n")
  }
  
  # Input validation
  if (!is.data.frame(traits)) {
    stop("Error: 'traits' must be a data frame")
  }
  if (!inherits(tree, "phylo")) {
    stop("Error: 'tree' must be a phylogenetic tree object")
  }
  if (bootstrap && (is.null(Sizesample) || Sizesample > nrow(traits))) {
    stop("Error: Invalid sample size for bootstrap")
  }
  
  # Initialize list to store bootstrap results
  bootstrap_results <- list()
  
  # Function to run single PGLS analysis
  run_single_pgls <- function(traits_data, tree_data, gene_numbers_data) {
    tryCatch({
      phenotypes <- colnames(traits_data)[c(spp.col.num, pheno.col.nums)]
      
      if (No.variables == 1) {
        out <- as.data.frame(t(apply(X = gene_numbers_data, 1, FUN = function(x) {
          traits2 <- tryCatch({
            temp <- cbind(traits_data[, c(phenotypes)], as.numeric(x))
            names(temp)[dim(temp)[2]] <- "GFS"
            temp
          }, error = function(e) {
            log_message(paste("Error in creating traits2:", e$message))
            return(NULL)
          })
          
          if (is.null(traits2)) return(rep(NA, 12))
          
          tree_pruned <- tryCatch({
            drop.tip(tree_data, names(x)[is.na(x)])
          }, error = function(e) {
            log_message(paste("Error in pruning tree:", e$message))
            return(NULL)
          })
          
          if (is.null(tree_pruned)) return(rep(NA, 12))
          
          formilin <- as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], 
                                                    collapse = " + "), sep = " ~ "))
          formilin2 <- as.formula(paste("~", phenotypes[1], sep = ""))
          
          pglsModel <- tryCatch({
            gls(model = formilin, 
                correlation = corBrownian(phy = tree_pruned, form = formilin2),
                method = "ML",
                data = traits2)
          }, error = function(e) {
            log_message(paste("PGLS model fitting failed:", e$message))
            return(NULL)
          })
          
          if (is.null(pglsModel)) return(rep(NA, 12))
          
          # Extract statistics
          stats <- tryCatch({
            c(anova(pglsModel)$`F-value`[1], 
              anova(pglsModel)$`p-value`[1],
              anova(pglsModel)$`F-value`[2], 
              anova(pglsModel)$`p-value`[2],
              summary(pglsModel)$tTable[1,1], 
              summary(pglsModel)$tTable[2,1],
              summary(pglsModel)$tTable[1,2], 
              summary(pglsModel)$tTable[2,2],
              summary(pglsModel)$tTable[1,3], 
              summary(pglsModel)$tTable[2,3],
              summary(pglsModel)$tTable[1,4], 
              summary(pglsModel)$tTable[2,4])
          }, error = function(e) {
            log_message(paste("Error extracting model statistics:", e$message))
            return(rep(NA, 12))
          })
          
          return(stats)
        })))
        
        # Set column names
        colnames(out) <- c("F_value Intercept", "p_value Intercept",
                           paste("F_value", phenotypes[2], sep=" "),
                           paste("p_value", phenotypes[2], sep=" "),
                           "Coefficient Intercept",
                           paste("Coefficient", phenotypes[2], sep=" "),
                           "Coef SE Intercept",
                           paste("Coef SE", phenotypes[2], sep=" "),
                           "Coef tval Intercept",
                           paste("Coef tval", phenotypes[2], sep=" "),
                           "Coef pval Intercept",
                           paste("Coef pval", phenotypes[2], sep=" "))
        
        return(out)
      } else {
        stop("Currently only implemented for No.variables = 1")
      }
    }, error = function(e) {
      log_message(paste("Error in PGLS analysis:", e$message))
      return(NULL)
    })
  }
  
  # Main analysis
  if (bootstrap) {
    valid_results <- 0
    
    for (i in 1:n_bootstrap) {
      log_message(paste("Starting bootstrap iteration", i))
      
      tryCatch({
        # Create copies of input data
        traits.cl <- traits
        tree.cl <- tree
        gene.numbers.cl <- gene.numbers
        
        # Sample species without replacement
        resampled.species <- sample(traits.cl[, spp.col.num], 
                                    replace = FALSE, 
                                    size = Sizesample)
        
        # Save sampled species
        write.csv(resampled.species, 
                  file.path(where.Save.it, 
                            paste0("Species_PGLSbootstrapNumber_", i, ".csv")))
        
        # Subset data
        traits.cl <- traits.cl[traits.cl[, spp.col.num] %in% resampled.species, , 
                               drop = FALSE]
        tree.cl <- drop.tip(tree.cl, 
                            setdiff(tree.cl$tip.label, traits.cl[, 1]))
        gene.numbers.cl <- gene.numbers.cl[, colnames(gene.numbers.cl) %in% 
                                             resampled.species, drop = FALSE]
        
        # Run PGLS
        result <- run_single_pgls(traits.cl, tree.cl, gene.numbers.cl)
        
        # Validate result
        if (!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
          bootstrap_results[[i]] <- as.matrix(result)  # Convert to matrix
          valid_results <- valid_results + 1
          
          # Save intermediate results
          saveRDS(bootstrap_results[[i]], 
                  file.path(where.Save.it, 
                            paste0("PGLS_bootstrap_", i, ".rds")))
          
          log_message(paste("Completed bootstrap iteration", i))
        } else {
          log_message(paste("Bootstrap iteration", i, "produced invalid results"))
        }
        
      }, error = function(e) {
        log_message(paste("Error in bootstrap iteration", i, ":", e$message))
        bootstrap_results[[i]] <- NULL
      })
    }
    
    # Check if we have enough valid results
    if (valid_results < 2) {
      stop(paste("Insufficient valid bootstrap results:", valid_results, 
                 "out of", n_bootstrap, "iterations succeeded"))
    }
    
    # Combine and summarize bootstrap results
    tryCatch({
      # Remove NULL results
      valid_results <- bootstrap_results[!sapply(bootstrap_results, is.null)]
      
      # Convert all results to matrices and stack them
      result_array <- simplify2array(lapply(valid_results, as.matrix))
      
      # Calculate summary statistics
      bootstrap_summary <- list(
        n_valid = length(valid_results),
        n_total = n_bootstrap,
        mean = apply(result_array, 1:2, mean, na.rm = TRUE),
        sd = apply(result_array, 1:2, sd, na.rm = TRUE)
      )
      
      # Add quantiles if possible
      tryCatch({
        bootstrap_summary$quantiles <- apply(result_array, 1:2, 
                                             function(x) quantile(x, probs = c(0.025, 0.975), 
                                                                  na.rm = TRUE))
      }, error = function(e) {
        log_message(paste("Warning: Could not calculate quantiles:", e$message))
      })
      
      # Save results
      saveRDS(bootstrap_summary, 
              file.path(where.Save.it, "PGLS_bootstrap_summary.rds"))
      
      log_message(sprintf("Analysis complete. Valid results: %d/%d", 
                          bootstrap_summary$n_valid, n_bootstrap))
      
      return(bootstrap_summary)
      
    }, error = function(e) {
      # Save debug information
      saveRDS(bootstrap_results, 
              file.path(where.Save.it, "bootstrap_results_debug.rds"))
      stop(paste("Error in summarizing bootstrap results:", e$message))
    })
    
  } else {
    # Run single PGLS analysis without bootstrap
    return(run_single_pgls(traits, tree, gene.numbers))
  }
}


result <- tryCatch({
  PGLS.GF.size(
    traits = traits.filtered,
    pheno.col.nums = 10,
    tree = tree,
    spp.col.num = 1,
    gene.numbers = gene.numbers.filtered,
    No.variables = 1,
    where.Save.it = dirSave,
    bootstrap = TRUE,
    n_bootstrap = 10,
    Sizesample = 50
  )
}, error = function(e) {
  message("Error in analysis: ", e$message)
  # Check the log file for detailed error information
  if (file.exists(paste(dirSave,"/pgls_analysis.log", sep = ""))) {
    cat("\nLog file contents:\n")
    cat(readLines(paste(dirSave,"/pgls_analysis.log", sep = "")), sep="\n")
  }
  return(NULL)
})


my_data_bootstrap <- readRDS("/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/PGLS_bootstrap_1.rds")
