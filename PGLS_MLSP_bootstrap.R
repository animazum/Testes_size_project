############################################
######## PGLS     PGLS      PGLS ###########
############################################

library(ape)
library(caper)
library(nlme)
library(parallel)

# testgf<-head(x = gene.numbers.filtered, n = dim(gene.numbers.filtered)[2])  ### test gene number file with 100 Orthogroups.

# write.tree(tree2, paste("/Users/mona/Dropbox/SSD/mammalian.phylo.",tree2$Nnode,".ssp.SSDchapter.tre",sep = ""))
################
##### PGLS #####
### function ###
################

### calculation average bodymass, the results are in absolut values.

PGLS.GF.size<- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it){
  phenotypes<<- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  log_print("#########################", hide_notes = T)
  log_print("##### PGLS PGLS PGLS ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  if (No.variables == 1){
    ###################################
    log_print("#### PGLS ONE Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,1,FUN = function(x){
      
      traits2<<-cbind(traits[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      tree<- drop.tip(tree,names(x)[is.na(x)]) 
      
      ## creating the formula for the model for 1 variable
      formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin2<-as.formula(paste("~",phenotypes[1], sep = ""))
      pglsModel<<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      
      if (inherits(pglsModel, "try-error")) 
        return(c(NA,NA))
      else
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2], 
                 summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2], 
                 summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4]))
    })))
    colnames(out)<<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2] , sep=" "), paste("p_value", phenotypes[2], sep=" "), 
                      "Coefficient Intercept", paste("Coefficient", phenotypes[2], sep=" "), "Coef SE Intercept", paste("Coef SE", phenotypes[2], sep=" "), 
                      "Coef tval Intercept", paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept", paste("Coef pval", phenotypes[2], sep=" "))
    
  } 
  else if (No.variables == 2) {
    ###################################
    log_print("#### PGLS TWO Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,MARGIN = 1,FUN = function(x){
      pheno1<<-phenotypes[3]
      pheno2<<-phenotypes[2]
      traits2<<-cbind(traits[,c("Spp.Name", pheno2, pheno1)],as.numeric(x)) #add phenotypes as needed
      # traits2<<-cbind(traits[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      tree<- drop.tip(tree,names(x)[is.na(x)])
      ## creating the formula for the model for 2 variables
      # formilin<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin<<-as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
      #log_print((paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ ")))
      # formilin2<-as.formula(paste("~",colnames(traits)[spp.col.num], sep = ""))
      formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
      #log_print((paste("~",colnames(traits)[spp.col.num], sep = "")))
      # pglsModel<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      pglsModel<<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      
      if (inherits(pglsModel, "try-error"))
        return(c(NA,NA))
      else
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2],
                 anova(pglsModel)$`F-value`[3], anova(pglsModel)$`p-value`[3], summary(pglsModel)$tTable[1,1],
                 summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[3,1], summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2],
                 summary(pglsModel)$tTable[3,2], summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[3,3],
                 summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4], summary(pglsModel)$tTable[3,4]))
    })))
    colnames(out)<<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[3] , sep=" "), paste("p_value", phenotypes[3], sep=" "),
                      paste("F_value", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "), "Coefficient Intercept", paste("Coefficient", phenotypes[3], sep=" "),
                      paste("Coefficient", phenotypes[2], sep=" "),  "Coef SE Intercept", paste("Coef SE", phenotypes[3], sep=" "), paste("Coef SE", phenotypes[2], sep=" "),
                      "Coef tval Intercept", paste("Coef tval", phenotypes[3], sep=" "), paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept",
                      paste("Coef pval", phenotypes[3], sep=" "), paste("Coef pval", phenotypes[2], sep=" "))
  }
  else if (No.variables > 2){
    ###################################
    ## PGLS more than two Variables ###
    ###################################
    log_print("At the moment the function only works with upto 2 variables, sorry :P", hide_notes = T)
  }
  for (i in 1:No.variables) {
    ###########################
    ### R^2 calculation ###
    ###########################
    log_print("Calc. Rs", hide_notes = T)
    sps<-nrow(traits2)
    DeFr <- sps - (1 + No.variables)
    Coef.tval1<-toString(paste("Coef tval", phenotypes[i+1], sep=" "))
    R1 <- toString(paste("R t value", phenotypes[i+1], sep=" "))
    R_t_value1 <- out[[Coef.tval1]] / (sqrt(((out[[Coef.tval1]])^2) + DeFr)) # tvalue / square root of(tvalue^2 + DF))
    out[[R1]] <<- R_t_value1
    ### benjamini correction
    log_print("Calc. benjamini correction", hide_notes = T)
    out[[paste("p.adjusted GFS vs.", phenotypes[i+1], sep = "")]] <<- p.adjust(out[,paste("p_value", phenotypes[i+1], sep=" ")], method = "fdr")
  }
  
  nameout<- paste("PGLS", paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits2)[1], "Spp", dim(out)[1], "GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  log_print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  log_print(paste("object created:", paste( "out", sep = "")), hide_notes = T)
  write.csv(out, paste(where.Save.it, nameout, sep = "/"))
}

### PGLS function
# dirSave<-"/Volumes/PortableSSD/Mammalian_orthologs_project/PLGS"
# dirSave<-"/Users/mona/Dropbox/SSD/PGLSgenomicafuncional/"
dirSave<- "~/Dropbox/TestesSize_SSD/Datos/"
dirSave<-path1

log_print(paste("You have ",dim(gene.numbers.filtered)[1], " gene families to input the PGLS", sep = ""), hide_notes = T)

write.csv(traits.filtered, "traits.filtered_SSDLog10testessize_paper.csv")

PGLS.GF.size(traits = traits.filtered, pheno.col.nums = c(10), No.variables = 1, 
             tree = tree, spp.col.num = 1, gene.numbers = gene.numbers.filtered, 
             where.Save.it = dirSave)

colnames(out)
rownames(out)


# Calculate counts and percentages
summary_data <- out %>%
  filter(`p_value SSD` < 0.05) %>% ### CHANGE THIS NAME TO THE TRAIT YOU ARE INTERESTED IN
  summarise(
    Positive = sum(`R t value SSD` > 0), ### ALSO CHANGE R VALUES
    Negative = sum(`R t value SSD` < 0)  ### ALSO CHANGE R VALUES
  ) %>%
  tidyr::pivot_longer(everything(), 
                      names_to = "Direction", 
                      values_to = "Count") %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Create horizontal 100% stacked barplot
Nombeplot<- paste("Gene Families Sig.Associated (p < 0.05)_", trait, "_PGLS.pdf", sep = "")
pdf(Nombeplot, width = 8, height = 6, useDingbats = FALSE)

ggplot(summary_data, aes(y = "Significant Genes", x = Percentage, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Negative" = "indianred", "Positive" = "steelblue")) +
  geom_text(aes(label = sprintf("%.1f%% \n(n=%d)", Percentage, Count)), 
            position = position_stack(vjust = 0.5),
            size = 4) +
  labs(title = paste("Proportion of Gene Families with p < 0.05", trait, "PGLS"),
       x = "Percentage",
       y = "") +
  theme_minimal() +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  theme(legend.position = "top",
        legend.title = element_blank())
dev.off()

### output object called 'out'
#traits.filtered$Spp.Name
#out<- read.csv("PGLS.Av.Body.mass.125.Spp.13753.GeneFams.Rs_benjamini.2022-08-03.csv")
#perms<- read.csv("PGLS.Av.Body.mass.125.Spp.13753.GeneFams.Rs_benjamini.2022-08-03.csv")
