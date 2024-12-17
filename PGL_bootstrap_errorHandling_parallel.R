library("parallel")
library("logr")
library("dplyr")
library("ape")
library("nlme")

rm(list=ls())

PGLS.GF.size.parallel <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it,
                                  bootstrap = FALSE, n_bootstrap = 100, Sizesample = NULL, cores = 1) {
  # Initialize logging
  log_file <- file.path(where.Save.it, "pgls_analysis.log")
  if (!file.exists(log_file)) {
    file.create(log_file)
  }
  
  log_message <- function(msg, hide_notes = TRUE) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message <- paste(timestamp, msg, sep = " - ")
    cat(message, "\n", file = log_file, append = TRUE)
    if (!hide_notes) cat(message, "\n")
  }
  
  # Input validation
  if (!is.data.frame(traits)) stop("Error: 'traits' must be a data frame")
  if (!inherits(tree, "phylo")) stop("Error: 'tree' must be a phylogenetic tree object")
  if (bootstrap && (is.null(Sizesample) || Sizesample > nrow(traits))) {
    stop("Error: Invalid sample size for bootstrap")
  }
  
  # Function to run single PGLS analysis
  run_single_pgls <- function(traits_data, tree_data, gene_numbers_data, pheno_cols) {
    tryCatch({
      if (No.variables != 1) stop("Currently only implemented for No.variables = 1")
      
      phenotypes <- colnames(traits_data)[c(spp.col.num, pheno_cols)]
      
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
          c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1],
            anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2],
            summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1],
            summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2],
            summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3],
            summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4])
        }, error = function(e) {
          log_message(paste("Error extracting model statistics:", e$message))
          return(rep(NA, 12))
        })
        
        return(stats)
      })))
      
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
    }, error = function(e) {
      log_message(paste("Error in PGLS analysis:", e$message))
      return(NULL)
    })
  }
  
  # Bootstrap analysis with parallel processing
  if (bootstrap) {
    log_message(paste("Starting parallel bootstrap analysis with", cores, "cores"))
    
    # Set up parallel processing
    nucleos <- min(detectCores() + cores - detectCores(), parallel::detectCores())
    cl <- makeCluster(nucleos)
    on.exit(stopCluster(cl))
    
    # Export necessary packages to all cores
    clusterEvalQ(cl, {
      library(ape)
      library(nlme)
      library(dplyr)
    })
    
    # Export necessary functions and data to all cores
    clusterExport(cl = cl,
                  varlist = c("run_single_pgls", "log_message", "No.variables",
                              "spp.col.num", "where.Save.it", "traits", "tree", 
                              "gene.numbers", "pheno.col.nums", "Sizesample"),
                  envir = environment())
    
    # Function to run single bootstrap iteration
    run_bootstrap <- function(i) {
      tryCatch({
        # Sample species without replacement
        resampled.species <- sample(traits[, spp.col.num],
                                    replace = FALSE,
                                    size = Sizesample)
        
        # Save sampled species
        write.csv(resampled.species,
                  file.path(where.Save.it,
                            paste0("Species_PGLSbootstrapNumber_", i, ".csv")))
        
        # Subset data
        traits.cl <- traits[traits[, spp.col.num] %in% resampled.species, , drop = FALSE]
        tree.cl <- drop.tip(tree, setdiff(tree$tip.label, traits.cl[, 1]))
        gene.numbers.cl <- gene.numbers[, colnames(gene.numbers) %in% resampled.species, drop = FALSE]
        
        # Run PGLS
        result <- run_single_pgls(traits.cl, tree.cl, gene.numbers.cl, pheno.col.nums)
        
        if (!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
          saveRDS(as.matrix(result),
                  file.path(where.Save.it,
                            paste0("PGLS_bootstrap_", i, ".rds")))
          return(as.matrix(result))
        }
        return(NULL)
      }, error = function(e) {
        log_message(paste("Error in bootstrap iteration", i, ":", e$message))
        return(NULL)
      })
    }
    
    # Run bootstrap iterations in parallel
    bootstrap_results <- parLapply(cl, 1:n_bootstrap, run_bootstrap)
    
    # Process results
    valid_results <- bootstrap_results[!sapply(bootstrap_results, is.null)]
    
    if (length(valid_results) < 2) {
      stop(paste("Insufficient valid bootstrap results:", length(valid_results),
                 "out of", n_bootstrap, "iterations succeeded"))
    }
    
    # Calculate summary statistics
    result_array <- simplify2array(valid_results)
    bootstrap_summary <- list(
      n_valid = length(valid_results),
      n_total = n_bootstrap,
      mean = apply(result_array, 1:2, mean, na.rm = TRUE),
      sd = apply(result_array, 1:2, sd, na.rm = TRUE),
      quantiles = apply(result_array, 1:2,
                        function(x) quantile(x, probs = c(0.025, 0.975),
                                             na.rm = TRUE))
    )
    
    saveRDS(bootstrap_summary,
            file.path(where.Save.it, "PGLS_bootstrap_summary.rds"))
    
    return(bootstrap_summary)
  } else {
    return(run_single_pgls(traits, tree, gene.numbers, pheno.col.nums))
  }
}


print("change always the ulimit of the computer $ulimit -s 21000 !!!!!!!!!!!!")

path1<- getwd()

Nombre.file<-c(paste(path1,"/Datos/traits_albus_61spp.csv", sep = ""))
traits.filtered<-read.csv(Nombre.file, stringsAsFactors = F)
tree <-read.tree(paste(path1,"/Datos/tree_albus_61spp.tree", sep = ""))

#*************   IT IS IMPORTANT TO HAVE THE SAME SPP ORDER BETWEEN TRAITS AND GENE NUMBERSS!!!!!!!   ********###
gene.numbers.filtered <- read.csv(paste(path1,"/Datos/gene_numbers_albus_61spp.csv", sep = ""),header = T, row.names = 1)
dirSave <- getwd()

result <- PGLS.GF.size.parallel(
  traits = traits.filtered,
  pheno.col.nums = 10,
  tree = tree,
  spp.col.num = 1,
  gene.numbers = gene.numbers.filtered,
  No.variables = 1,
  where.Save.it = dirSave,
  bootstrap = TRUE,
  n_bootstrap = 10,
  Sizesample = 50,
  cores = 4  # Adjust based on your system
)

