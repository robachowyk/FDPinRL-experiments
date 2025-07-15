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
#               LAUNCH THE PROCEDURE POTENTIALLY WITH BLOCKING                 #
################################################################################

gc()

################################################################################
#                         IMPORT ENCODED DATA SHIW                             #
################################################################################

DataSource = "SHIW"
Method = "synthpop"
Blocking = "AREA3"   # "None", "IREG", "AREA3"
BlockingValue = 2   # "NA", "IREG number", "AREA3 number"

source("SHIW_Block.R")

gc()