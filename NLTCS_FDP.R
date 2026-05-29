################################################################################
#                              IMPORT LIBRARIES                                #
################################################################################

library(progress)
library(Matrix)
library(testit)
library(FlexRL)
library(rcompanion)
library(fastLink)
library(BRL) # remotes::install_github("robachowyk/BRLrobachowyk")
library(mice)
library(synthpop)
library(arf)
source("Synth.R")

################################################################################
#                LAUNCH THE PROCEDURE POTENTIALLY WITH BLOCKING                #
################################################################################

gc()

################################################################################
#                         IMPORT ENCODED DATA NLTCS                            #
################################################################################

DataSource = "NLTCS"
Method = "arf"
Blocking = "REG_CODE"   # "None", "REG_CODE"
BlockingValue = 29   # "NA", "REG_CODE number"

source("NLTCS_Block.R")

gc()