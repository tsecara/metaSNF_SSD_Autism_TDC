####### Running combat on all data types before running the metaSNF analysis 

##Install packages 
install.packages("devtools")
library(devtools)
install_github("jfortin1/neuroCombatData")
install_github("jfortin1/neuroCombat_Rpackage")
library(neuroCombat)

##########################1###############################
# Loading in the different data types
CT <- read.csv(".../cort_thick_combinedFINAL.csv", header=TRUE)
SA <- read.csv(".../cort_SAavg_combined_regressedICV_FINAL.csv", header=TRUE)
VOL <- read.csv(".../cort_vol_RL_combined_regressedICV_FINAL.csv", header=TRUE)
FA <- read.csv(".../FA_val_combinedFINAL.csv", header=TRUE)
MD <- read.csv(".../MD_val_combinedFINAL.csv", header=TRUE)
clinical_scog <- read.csv(".../scog_metaSNF_FINAL.csv", header=TRUE)
clinical_nc <- read.csv(".../ncog_metaSNF_FINAL.csv", header=TRUE)
EA <- read.csv(".../EA_task_tmaps_combined_FULLsample.csv", header=TRUE)
RS <- read.csv(".../RS_data_combined_noGSR.csv", header=TRUE)
orignial_demo <- read.csv(".../combined_demo_ALL_IRI_additioanal_SPASDcase.csv", header = TRUE)
orignial_demo <- orignial_demo[orignial_demo$age <= 35, ]
validation_data <- read.csv(".../demo_combined_reduced_FINAL3.csv", header = TRUE)

###Need to do this step so that there are no participants older than 35
clinical_scog <- clinical_scog[(clinical_scog$participant_id %in% orignial_demo$participant_id), ]

###Reducing the participants based on what is included in the SNF analysis 
sublist <- read.csv(".../sublist/FINAL_sample.csv")

#Based on the sublist, remove participants in preparation for parameter iterations: 
CT <- CT[(CT$participant_id %in% sublist$participant_id), ]
SA <- SA[(SA$participant_id %in% sublist$participant_id), ]
VOL <- VOL[(VOL$participant_id %in% sublist$participant_id), ]
FA <- FA[(FA$participant_id %in% sublist$participant_id), ]
MD <- MD[(MD$participant_id %in% sublist$participant_id), ]
clinical_scog <- clinical_scog[(clinical_scog$participant_id %in% sublist$participant_id), ]
clinical_nc <- clinical_nc[(clinical_nc$participant_id %in% sublist$participant_id), ]
EA <- EA[(EA$participant_id %in% sublist$participant_id), ]
RS <- RS[(RS$participant_id %in% sublist$participant_id), ]...
validation_data <- validation_data[(validation_data$participant_id %in% sublist$participant_id),]
validation_data <- validation_data[(validation_data$participant_id %in% RS$participant_id),]

###### Removing participant id columns 
rownames(validation_data) <- validation_data$participant_id
validation_data <- validation_data[,-1]
rownames(CT) <- CT$participant_id
CT <- CT[,-1]
rownames(SA) <- SA$participant_id
SA <- SA[,-1]
rownames(VOL) <- VOL$participant_id
VOL <- VOL[,-1]
rownames(FA) <- FA$participant_id
FA <- FA[,-1]
rownames(MD) <- MD$participant_id
MD <- MD[,-1]
rownames(clinical_scog) <- clinical_scog$participant_id
clinical_scog <- clinical_scog[,-1]
rownames(clinical_nc) <- clinical_nc$participant_id
clinical_nc <- clinical_nc[,-1]
rownames(EA) <- EA$participant_id
EA <- EA[,-1]
rownames(RS) <- RS$participant_id
RS <- RS[,-1]
rownames(IMOBS) <- IMOBS$participant_id
IMOBS <- IMOBS[,-1]

######################################2#############################################
# RUNNING COMBAT ON DATA 

#Must transpose so that the participants are in the columns and the features are in the rows  
CT <- t(as.matrix(CT)) #first convert to a matrix, and then transpose 
SA <- t(as.matrix(SA))
VOL <- t(as.matrix(VOL))
FA <- t(as.matrix(FA))
MD <- t(as.matrix(MD))
clinical_scog <- t(as.matrix(clinical_scog))
clinical_nc <- t(as.matrix(clinical_nc))
EA <- t(as.matrix(EA))
RS <- t(as.matrix(RS))
IMOBS <- t(as.matrix(IMOBS))

###Ensure that the validation data is in the right form 
class(validation_data$group)
validation_data$group <- as.factor(validation_data$group)
class(validation_data$scanner)
validation_data$scanner <- as.factor(validation_data$scanner)

# mod is a design matrix specifying biological covariates that should be protected
modcombat <- model.matrix(~group +sex + age + 
                            rmet_total +	er40_total +
                            tasit1_total + tasit2_sin +	tasit2_ssar +	tasit2_psar +	tasit3_lies +	tasit3_sar +
                            mccb_pspeed +	mccb_atvig + mccb_wmem + mccb_verblearn	+ mccb_vislearn	+ mccb_reasonps	+ mccb_scog	+
                            bsfs_total, data=validation_data)

# R run ComBat
# batch is a vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for
CT_combat <- neuroCombat(dat=CT, batch=c(validation_data$scanner), mod=modcombat)
SA_combat <- neuroCombat(dat=SA, batch=c(validation_data$scanner), mod=modcombat)
VOL_combat <- neuroCombat(dat=VOL, batch=c(validation_data$scanner), mod=modcombat)
FA_combat <- neuroCombat(dat=FA, batch=c(validation_data$scanner), mod=modcombat)
MD_combat <- neuroCombat(dat=MD, batch=c(validation_data$scanner), mod=modcombat)
clinical_scog_combat <- neuroCombat(dat=clinical_scog, batch=c(validation_data$scanner), mod=modcombat)
clinical_nc_combat <- neuroCombat(dat=clinical_nc, batch=c(validation_data$scanner), mod=modcombat)
EA_combat <- neuroCombat(dat=EA, batch=c(validation_data$scanner), mod=modcombat)
RS_combat <- neuroCombat(dat=RS, batch=c(validation_data$scanner), mod=modcombat)

# transpose the harmonized data matrix
CT <- as.data.frame(t(CT_combat$dat.combat))
SA <- as.data.frame(t(SA_combat$dat.combat))
VOL <- as.data.frame(t(VOL_combat$dat.combat))
FA <- as.data.frame(t(FA_combat$dat.combat))
MD <- as.data.frame(t(MD_combat$dat.combat))
clinical_scog <- as.data.frame(t(clinical_scog_combat$dat.combat))
clinical_nc <- as.data.frame(t(clinical_nc_combat$dat.combat))
EA <- as.data.frame(t(EA_combat$dat.combat))
RS <- as.data.frame(t(RS_combat$dat.combat))

# Install and load the dplyr package if not already installed
install.packages("dplyr")
library(dplyr)
library(tibble)

# Add row names as a new column named "participant_id"
CT <- CT %>%
  rownames_to_column(var = "participant_id")
SA <- SA %>%
  rownames_to_column(var = "participant_id")
VOL <- VOL %>%
  rownames_to_column(var = "participant_id")
FA <- FA %>%
  rownames_to_column(var = "participant_id")
MD <- MD %>%
  rownames_to_column(var = "participant_id")
clinical_scog <- clinical_scog %>%
  rownames_to_column(var = "participant_id")
clinical_nc <- clinical_nc %>%
  rownames_to_column(var = "participant_id")
EA <- EA %>%
  rownames_to_column(var = "participant_id")
RS <- RS %>%
  rownames_to_column(var = "participant_id")

#Export COMBAT data for further analysis 
write.csv(CT, ".../combat_data3/CT_S200.csv")
write.csv(SA, ".../combat_data3/SA.csv")
write.csv(VOL, ".../combat_data3/VOL.csv")
write.csv(FA, ".../combat_data3/FA.csv")
write.csv(MD, ".../combat_data3/MD.csv")
write.csv(clinical_scog, ".../combat_data3/clinical_scog.csv")
write.csv(clinical_nc, ".../combat_data3/clinical_nc.csv")
write.csv(EA, ".../combat_data3/EA.csv")
write.csv(RS, ".../combat_data3/RS_noGSR.csv")
