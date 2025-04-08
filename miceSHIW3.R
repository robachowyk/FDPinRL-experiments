################################################################################
#                RUN APPLICATIONS FOR THE RECORD LINKAGE TASK                  #
#          AND COMPARE FLEXRL WITH A NAIVE APPROACH ON THE FULL DATA           #
################################################################################

################################################################################
#                                 IMPORTS                                      #
################################################################################
# library(progress)
# library(Matrix)
# library(testit)
# library(FlexRL)
# library(rcompanion)
# library(fastLink)
# remotes::install_github("robachowyk/BRLrobachowyk")
# library(mice)
# library(synthpop)
# library(arf)
# source("UntitledSynth.R")

################################################################################
#                       CREATE A TABLE TO SAVE RESULTS                         #
################################################################################

setwd("draftwork/")

DataSource = "SHIW"
Area = "3"
Method = "mice"

newDirectory = sprintf("%s-Area: %s-%s---%s", DataSource, Area, Method, Sys.time())
dir.create(newDirectory)

thef1score = data.frame(matrix(0,20,4))
colnames(thef1score) = c("Naive","BRL","FastLink","FlexRL")

thesensitivity = data.frame(matrix(0,20,4))
colnames(thesensitivity) = c("Naive","BRL","FastLink","FlexRL")

thetp = data.frame(matrix(0,20,4))
colnames(thetp) = c("Naive","BRL","FastLink","FlexRL")

thefp = data.frame(matrix(0,20,4))
colnames(thefp) = c("Naive","BRL","FastLink","FlexRL")

thefn = data.frame(matrix(0,20,4))
colnames(thefn) = c("Naive","BRL","FastLink","FlexRL")

thesynth = data.frame(matrix(0,20,4))
colnames(thesynth) = c("Naive","BRL","FastLink","FlexRL")

true_fdr = data.frame(matrix(0,20,4))
colnames(true_fdr) = c("Naive","BRL","FastLink","FlexRL")

prob_fdr = data.frame(matrix(0,20,4))
colnames(true_fdr) = c("Naive","BRL","FastLink","FlexRL")

estimated_fdr_1 = data.frame(matrix(0,20,4))
colnames(estimated_fdr_1) = c("Naive","BRL","FastLink","FlexRL")

estimated_fdr_2 = data.frame(matrix(0,20,4))
colnames(estimated_fdr_2) = c("Naive","BRL","FastLink","FlexRL")

estimated_fdr_3 = data.frame(matrix(0,20,4))
colnames(estimated_fdr_3) = c("Naive","BRL","FastLink","FlexRL")

estimated_fdr_4 = data.frame(matrix(0,20,4))
colnames(estimated_fdr_4) = c("Naive","BRL","FastLink","FlexRL")

assumption1_filled = data.frame(matrix(0,20,4))
colnames(assumption1_filled) = c("Naive","BRL","FastLink","FlexRL")

assumption2_filled = data.frame(matrix(0,20,4))
colnames(assumption2_filled) = c("Naive","BRL","FastLink","FlexRL")

ratio_duplicata = data.frame(matrix(0,20,1))
colnames(ratio_duplicata) = c("Naive")

################################################################################
#                              IMPORT ENCODED DATA                             #
################################################################################

df1 = read.csv("SHIWData/df2016.csv", row.names = 1)
df2 = read.csv("SHIWData/df2020.csv", row.names = 1)

df1 = df1[df1$AREA3==as.integer(Area), ]
df2 = df2[df2$AREA3==as.integer(Area), ]

################################################################################
#                               PREPARE THE DATA                               #
################################################################################

PIVs_config = list( ANASCI     = list(stable = TRUE),
                    SESSO      = list(stable = TRUE),
                    STUDIO     = list(stable = TRUE),
                    STACIV     = list(stable = TRUE),
                    IREG       = list(stable = TRUE) )
PIVs = names(PIVs_config)
PIVs_stable = sapply(PIVs_config, function(x) x$stable)
controlOnMistakes = c(TRUE, TRUE, FALSE, FALSE, FALSE)

### FILTER THE DATA ON THE INTERSECTING SUPPORT OF THE PIVS

idColumn = "ID"

### REMOVE DUPLICATES
df1 = df1[ !df1[,idColumn] %in% df1[,idColumn][duplicated(df1[,idColumn])] , ]
df2 = df2[ !df2[,idColumn] %in% df2[,idColumn][duplicated(df2[,idColumn])] , ]

### FILTER COLUMNS
for(i in 1:length(PIVs)){
  intersect_support_piv = intersect( unique(df1[,PIVs[i]]), unique(df2[,PIVs[i]]) )
  df1 = df1[df1[,PIVs[i]] %in% c(NA,intersect_support_piv),]
  df2 = df2[df2[,PIVs[i]] %in% c(NA,intersect_support_piv),]
}

# for(i in 1:length(PIVs)){
#   for(j in 1:length(PIVs)){
#     if(i>j){
#       pivsi = PIVs[i]
#       pivsj = PIVs[j]
#       contingency_table_A = table(df1[,pivsi], df1[,pivsj])
#       cramers_v_A <- cramerV(contingency_table_A)
#       contingency_table_B = table(df2[,pivsi], df2[,pivsj])
#       cramers_v_B <- cramerV(contingency_table_B)
#       cat("PIVs:", pivsi, "&", pivsj, "\nCramerV in A:", cramers_v_A, "\nCramerV in B:", cramers_v_B, "\n")
#     }
#   }
# }


# ENCODING
if(nrow(df2)>nrow(df1)){
  encodedA = df1
  encodedB = df2
}else{
  encodedB = df1
  encodedA = df2
}
levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(encodedA[,x], encodedB[,x])))))
for(i in 1:length(PIVs))
{
  encodedA[,PIVs[i]] = as.factor(as.numeric(factor(as.character(encodedA[,PIVs[i]]), levels=levels_PIVs[[i]])))
  encodedB[,PIVs[i]] = as.factor(as.numeric(factor(as.character(encodedB[,PIVs[i]]), levels=levels_PIVs[[i]])))
}
nvalues = sapply(levels_PIVs, length)
nbrRealRecordsB = nrow(encodedB)
nbrRealRecordsA = nrow(encodedA)
encodedA$localIndex = 1:nbrRealRecordsA
encodedB$localIndex = 1:nbrRealRecordsB
encodedA$source = "A"
encodedB$source = "B"

# IDENTIFICATION
linkedID = intersect(encodedA[,idColumn], encodedB[,idColumn])
Nlinks = length( linkedID )
### TRUE LINKAGE STRUCTURE
Delta = data.frame( matrix(0, nrow=0, ncol=2) )
colnames(Delta) = c("idxA", "idxB")
### AGREEMENTS IN THE TRUE LINKS
recapAgreementsTrueLinks = data.frame(matrix(0, nrow=0, ncol=length(PIVs)))
colnames(recapAgreementsTrueLinks) = PIVs
for (i in 1:Nlinks)
{
  id = linkedID[i]
  idA = which(encodedA[,idColumn] == id)
  idB = which(encodedB[,idColumn] == id)
  Delta = rbind(Delta, cbind(idxA=encodedA[idA,"localIndex"],idxB=encodedB[idB,"localIndex"]))
  recapAgreementsTrueLinks = rbind(recapAgreementsTrueLinks, encodedA[idA,PIVs] == encodedB[idB,PIVs])
}
recapAgreementsTrueLinks = colSums( recapAgreementsTrueLinks, na.rm = TRUE ) / nrow( recapAgreementsTrueLinks )
true_pairs = do.call(paste, c(Delta[,c("idxA","idxB")], list(sep="_")))

### MISSING IN THE WHOLE DATASETS (ENCODED WITH NA)
NA_A = ( colSums(is.na(encodedA)) / nrow(encodedA) )[PIVs]
NA_B = ( colSums(is.na(encodedB)) / nrow(encodedB) )[PIVs]

### MISSING IN THE TRUELINKS (ENCODED WITH NA)
NA_A_true = ( colSums(is.na(encodedA[encodedA[,idColumn] %in% linkedID,])) / Nlinks )[PIVs]
NA_B_true = ( colSums(is.na(encodedB[encodedB[,idColumn] %in% linkedID,])) / Nlinks )[PIVs]

### UNIQUE VALUES IN THE WHOLE DATASETS (ENCODED WITH NA)
unique_values = sapply( rbind(encodedA[,PIVs],encodedB[,PIVs]), function(x) length(unique(x[!(is.na(x))])) )

### STORY TELLING
df = data.frame( rbind(NA_A, NA_B, NA_A_true, NA_B_true, recapAgreementsTrueLinks, unique_values, nbrRealRecordsA, nbrRealRecordsB, Nlinks) )
colnames(df) = PIVs
rownames(df) = c("NaN in A", "NaN in B", "NaN in A true", "NaN in B true", "agreements btw A and B true links", "unique values", "size A", "size B", "Nlinks")
write.csv(df, file.path(newDirectory, "dataset_recapstory.csv"))

print(nbrRealRecordsA)
print(nbrRealRecordsB)
print(Nlinks)

for(iterator in 1:10){
  
  ### SYNTHESISE
  
  subsample_size = as.integer(0.15*nbrRealRecordsB)
  
  Newdata = synthesise("mice", encodedA, encodedB, PIVs, as.integer(0.15*nbrRealRecordsB))
  NewA = Newdata$NewA
  NewB = Newdata$NewB
  
  # Newdata = synthesise("arf", encodedA, encodedB, PIVs, as.integer(0.20*nbrRealRecordsB))
  # NewA = Newdata$NewA
  # NewB = Newdata$NewB
  # 
  # Newdata = synthesise("synthpop", encodedA, encodedB, PIVs, as.integer(0.20*nbrRealRecordsB))
  # NewA = Newdata$NewA
  # NewB = Newdata$NewB
  
  ### BRL
  print("BRL starts")
  myCompData <- BRL::compareRecords(NewB, NewA, flds=PIVs, types=rep("bi",length(PIVs)))
  chain <- BRL::bipartiteGibbs(myCompData, nIter=1000, a=1, b=1, aBM=1, bBM=1)
  names(chain)
  write.csv(chain$Z, file.path(newDirectory, sprintf("chain Z BRL_%s_%s_%s_%s.csv", iterator, nbrRealRecordsA, nbrRealRecordsB, Sys.time())))
  Zhat <- BRL::linkRecords(chain$Z[,setdiff(1:1000,seq_len(100)),drop=FALSE], n1=nrow(NewB), lFNM=1, lFM1=1, lFM2=2, lR=Inf)
  # lower fdr: λFNM < λFM1, λFM2
  idxA = which( Zhat <= nrow(NewB) )
  idxB = Zhat[ Zhat <= nrow(NewB) ]
  DeltaResult = data.frame(cbind(idxA=NewA[idxA,"localIndex"],idxB=NewB[idxB,"localIndex"]))
  if(nrow(DeltaResult)>1){
    linkedpairsA = DeltaResult[,"idxA"]
    linkedpairsB = DeltaResult[,"idxB"]
    real = (DeltaResult$idxB<=nbrRealRecordsB)&(DeltaResult$idxA<=nbrRealRecordsA)
    linkedpairsAreal = linkedpairsA[real]
    linkedpairsBreal = linkedpairsB[real]
    linked_real_pairs = do.call(paste, c(DeltaResult[real,c("idxA","idxB")], list(sep="_")))
    synthfalsepositive = sum(!real)
    truepositive = length( intersect(linked_real_pairs, true_pairs) )
    falsepositive = length( setdiff(linked_real_pairs, true_pairs) )
    falsenegative = length( setdiff(true_pairs, linked_real_pairs) )
    precision = truepositive / (truepositive + falsepositive)
    sensitivity = truepositive / (truepositive + falsenegative)
    f1score = 2 * (precision * sensitivity) / (precision + sensitivity)
    synthfdrhat1 = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat2 = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    synthfdrhat3 = (synthfalsepositive * (nbrRealRecordsA/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat4 = (synthfalsepositive * (1 + (nbrRealRecordsA/subsample_size))) / nrow(DeltaResult)
    assumption1 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    assumption2 = synthfalsepositive/min(c(nbrRealRecordsA,subsample_size)) - falsepositive/min(c(nbrRealRecordsA,nbrRealRecordsB))
    # SAVE RESULTS
    thef1score[iterator,"BRL"] = f1score
    thesensitivity[iterator,"BRL"] = sensitivity
    thetp[iterator,"BRL"] = truepositive
    thefp[iterator,"BRL"] = falsepositive
    thefn[iterator,"BRL"] = falsenegative
    thesynth[iterator,"BRL"] = synthfalsepositive
    true_fdr[iterator,"BRL"] = 1-precision
    prob_fdr[iterator,"BRL"] = NA
    estimated_fdr_1[iterator,"BRL"] = synthfdrhat1
    estimated_fdr_2[iterator,"BRL"] = synthfdrhat2
    estimated_fdr_3[iterator,"BRL"] = synthfdrhat3
    estimated_fdr_4[iterator,"BRL"] = synthfdrhat4
    assumption1_filled[iterator,"BRL"] = abs(assumption1)
    assumption2_filled[iterator,"BRL"] = abs(assumption2)
    write.csv(DeltaResult, file.path(newDirectory, sprintf("DeltaResult_BRL_%s_%s_%s_%s.csv", iterator, nbrRealRecordsA, nbrRealRecordsB, Sys.time())))
  }else{
    true_fdr[iterator,"BRL"] = "nolink"
  }
  print("BRL done")
  
  ### FastLink
  print("FastLink starts")
  matches.out = fastLink(NewA[,PIVs], NewB[,PIVs], PIVs, n.cores=8, threshold.match=0.5, tol.em=1e-06, return.all=TRUE)
  idxA = matches.out$matches$inds.a
  idxB = matches.out$matches$inds.b
  DeltaResult = data.frame(cbind(idxA=NewA[idxA,"localIndex"],idxB=NewB[idxB,"localIndex"],probaLink=matches.out$posterior))
  if(nrow(DeltaResult)>1){
    DeltaProba = DeltaResult$probaLink
    fdrhat = 1 - sum(DeltaProba[DeltaProba>0.5])/sum(DeltaProba>0.5)
    
    DeltaResult = DeltaResult[DeltaResult$probaLink>0.5,]
    
    linkedpairsA = DeltaResult[,"idxA"]
    linkedpairsB = DeltaResult[,"idxB"]
    real = (DeltaResult$idxB<=nbrRealRecordsB)&(DeltaResult$idxA<=nbrRealRecordsA)
    linkedpairsAreal = linkedpairsA[real]
    linkedpairsBreal = linkedpairsB[real]
    linked_real_pairs = do.call(paste, c(DeltaResult[real,c("idxA","idxB")], list(sep="_")))
    synthfalsepositive = sum(!real)
    truepositive = length( intersect(linked_real_pairs, true_pairs) )
    falsepositive = length( setdiff(linked_real_pairs, true_pairs) )
    falsenegative = length( setdiff(true_pairs, linked_real_pairs) )
    precision = truepositive / (truepositive + falsepositive)
    sensitivity = truepositive / (truepositive + falsenegative)
    f1score = 2 * (precision * sensitivity) / (precision + sensitivity)
    synthfdrhat1 = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat2 = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    synthfdrhat3 = (synthfalsepositive * (nbrRealRecordsA/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat4 = (synthfalsepositive * (1 + (nbrRealRecordsA/subsample_size))) / nrow(DeltaResult)
    assumption1 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    assumption2 = synthfalsepositive/min(c(nbrRealRecordsA,subsample_size)) - falsepositive/min(c(nbrRealRecordsA,nbrRealRecordsB))
    # SAVE RESULTS
    thef1score[iterator,"FastLink"] = f1score
    thesensitivity[iterator,"FastLink"] = sensitivity
    thetp[iterator,"FastLink"] = truepositive
    thefp[iterator,"FastLink"] = falsepositive
    thefn[iterator,"FastLink"] = falsenegative
    thesynth[iterator,"FastLink"] = synthfalsepositive
    true_fdr[iterator,"FastLink"] = 1-precision
    prob_fdr[iterator,"FastLink"] = fdrhat
    estimated_fdr_1[iterator,"FastLink"] = synthfdrhat1
    estimated_fdr_2[iterator,"FastLink"] = synthfdrhat2
    estimated_fdr_3[iterator,"FastLink"] = synthfdrhat3
    estimated_fdr_4[iterator,"FastLink"] = synthfdrhat4
    assumption1_filled[iterator,"FastLink"] = abs(assumption1)
    assumption2_filled[iterator,"FastLink"] = abs(assumption2)
    write.csv(DeltaResult, file.path(newDirectory, sprintf("DeltaResult_FastLink_%s_%s_%s_%s.csv", iterator, nbrRealRecordsA, nbrRealRecordsB, Sys.time())))
  }else{
    true_fdr[iterator,"FastLink"] = "nolink"
  }
  print("FastLink done")
  
  ### ENCODE THE MISSING VALUES
  NewA[,PIVs][ is.na(NewA[,PIVs]) ] = 0
  NewB[,PIVs][ is.na(NewB[,PIVs]) ] = 0
  
  ### LAUNCH SIMPLISTIC APPROACH
  print("Naive starts")
  DeltaResult = launchNaive(PIVs, NewA[,PIVs], NewB[,PIVs])
  if(nrow(DeltaResult)>1){
    linkedpairsA = DeltaResult[,"idxA"]
    linkedpairsB = DeltaResult[,"idxB"]
    real = (DeltaResult$idxB<=nbrRealRecordsB)&(DeltaResult$idxA<=nbrRealRecordsA)
    linkedpairsAreal = linkedpairsA[real]
    linkedpairsBreal = linkedpairsB[real]
    linked_real_pairs = do.call(paste, c(DeltaResult[real,c("idxA","idxB")], list(sep="_")))
    synthfalsepositive = sum(!real)
    truepositive = length( intersect(linked_real_pairs, true_pairs) )
    falsepositive = length( setdiff(linked_real_pairs, true_pairs) )
    falsenegative = length( setdiff(true_pairs, linked_real_pairs) )
    precision = truepositive / (truepositive + falsepositive)
    sensitivity = truepositive / (truepositive + falsenegative)
    f1score = 2 * (precision * sensitivity) / (precision + sensitivity)
    synthfdrhat1 = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat2 = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    synthfdrhat3 = (synthfalsepositive * (nbrRealRecordsA/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat4 = (synthfalsepositive * (1 + (nbrRealRecordsA/subsample_size))) / nrow(DeltaResult)
    assumption1 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    assumption2 = synthfalsepositive/min(c(nbrRealRecordsA,subsample_size)) - falsepositive/min(c(nbrRealRecordsA,nbrRealRecordsB))
    # SAVE RESULTS
    thef1score[iterator,"Naive"] = f1score
    thesensitivity[iterator,"Naive"] = sensitivity
    thetp[iterator,"Naive"] = truepositive
    thefp[iterator,"Naive"] = falsepositive
    thefn[iterator,"Naive"] = falsenegative
    thesynth[iterator,"Naive"] = synthfalsepositive
    true_fdr[iterator,"Naive"] = 1-precision
    prob_fdr[iterator,"Naive"] = NA
    estimated_fdr_1[iterator,"Naive"] = synthfdrhat1
    estimated_fdr_2[iterator,"Naive"] = synthfdrhat2
    estimated_fdr_3[iterator,"Naive"] = synthfdrhat3
    estimated_fdr_4[iterator,"Naive"] = synthfdrhat4
    assumption1_filled[iterator,"Naive"] = abs(assumption1)
    assumption2_filled[iterator,"Naive"] = abs(assumption2)
    ratio_duplicata[iterator,"Naive"] = length(linked_real_pairs) / (nbrRealRecordsA*nbrRealRecordsB)
    write.csv(DeltaResult, file.path(newDirectory, sprintf("DeltaResult_Naive_%s_%s_%s_%s.csv", iterator, nbrRealRecordsA, nbrRealRecordsB, Sys.time())))
  }else{
    true_fdr[iterator,"Naive"] = "nolink"
  }
  print("Naive done")
  
  ### LAUNCH FLEXRL
  print("FlexRL starts")
  data = list(A = NewA,
              B = NewB,
              Nvalues = nvalues,
              PIVs_config = PIVs_config,
              controlOnMistakes = controlOnMistakes,
              sameMistakes = TRUE,
              phiMistakesAFixed = FALSE,
              phiMistakesBFixed = FALSE,
              phiForMistakesA = c(NA,NA,NA,NA,NA),
              phiForMistakesB = c(NA,NA,NA,NA,NA))
  fit = stEM( data = data,
              StEMIter = 50,
              StEMBurnin = 25,
              GibbsIter = 40,
              GibbsBurnin = 20,
              musicOn = FALSE,
              newDirectory = newDirectory,
              saveInfoIter = FALSE )
  DeltaResult = fit$Delta
  colnames(DeltaResult) = c("idxA","idxB","probaLink")
  if(nrow(DeltaResult)>1){
    DeltaProba = DeltaResult$probaLink
    fdrhat = 1 - sum(DeltaProba[DeltaProba>0.5])/sum(DeltaProba>0.5)
    
    DeltaResult = DeltaResult[DeltaResult$probaLink>0.5,]
    
    linkedpairsA = DeltaResult[,"idxA"]
    linkedpairsB = DeltaResult[,"idxB"]
    real = (DeltaResult$idxB<=nbrRealRecordsB)&(DeltaResult$idxA<=nbrRealRecordsA)
    linkedpairsAreal = linkedpairsA[real]
    linkedpairsBreal = linkedpairsB[real]
    linked_real_pairs = do.call(paste, c(DeltaResult[real,c("idxA","idxB")], list(sep="_")))
    synthfalsepositive = sum(!real)
    truepositive = length( intersect(linked_real_pairs, true_pairs) )
    falsepositive = length( setdiff(linked_real_pairs, true_pairs) )
    falsenegative = length( setdiff(true_pairs, linked_real_pairs) )
    precision = truepositive / (truepositive + falsepositive)
    sensitivity = truepositive / (truepositive + falsenegative)
    f1score = 2 * (precision * sensitivity) / (precision + sensitivity)
    synthfdrhat1 = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat2 = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    synthfdrhat3 = (synthfalsepositive * (nbrRealRecordsA/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat4 = (synthfalsepositive * (1 + (nbrRealRecordsA/subsample_size))) / nrow(DeltaResult)
    assumption1 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    assumption2 = synthfalsepositive/min(c(nbrRealRecordsA,subsample_size)) - falsepositive/min(c(nbrRealRecordsA,nbrRealRecordsB))
    # SAVE RESULTS
    thef1score[iterator,"FlexRL"] = f1score
    thesensitivity[iterator,"FlexRL"] = sensitivity
    thetp[iterator,"FlexRL"] = truepositive
    thefp[iterator,"FlexRL"] = falsepositive
    thefn[iterator,"FlexRL"] = falsenegative
    thesynth[iterator,"FlexRL"] = synthfalsepositive
    true_fdr[iterator,"FlexRL"] = 1-precision
    prob_fdr[iterator,"FlexRL"] = fdrhat
    estimated_fdr_1[iterator,"FlexRL"] = synthfdrhat1
    estimated_fdr_2[iterator,"FlexRL"] = synthfdrhat2
    estimated_fdr_3[iterator,"FlexRL"] = synthfdrhat3
    estimated_fdr_4[iterator,"FlexRL"] = synthfdrhat4
    assumption1_filled[iterator,"FlexRL"] = abs(assumption1)
    assumption2_filled[iterator,"FlexRL"] = abs(assumption2)
    write.csv(DeltaResult, file.path(newDirectory, sprintf("DeltaResult_FlexRL_%s_%s_%s_%s.csv", iterator, nbrRealRecordsA, nbrRealRecordsB, Sys.time())))
  }else{
    true_fdr[iterator,"FlexRL"] = "nolink"
  }
  print("FlexRL done")
  
}

write.csv(thef1score, file.path(newDirectory, "thef1score.csv"))
write.csv(thesensitivity, file.path(newDirectory, "thesensitivity.csv"))
write.csv(thetp, file.path(newDirectory, "thetp.csv"))
write.csv(thefp, file.path(newDirectory, "thefp.csv"))
write.csv(thefn, file.path(newDirectory, "thefn.csv"))
write.csv(thesynth, file.path(newDirectory, "thesynth.csv"))
write.csv(true_fdr, file.path(newDirectory, "true_fdr.csv"))
write.csv(prob_fdr, file.path(newDirectory, "prob_fdr.csv"))
write.csv(estimated_fdr_1, file.path(newDirectory, "estimated_fdr_1.csv"))
write.csv(estimated_fdr_2, file.path(newDirectory, "estimated_fdr_2.csv"))
write.csv(estimated_fdr_3, file.path(newDirectory, "estimated_fdr_3.csv"))
write.csv(estimated_fdr_4, file.path(newDirectory, "estimated_fdr_4.csv"))
write.csv(assumption1_filled, file.path(newDirectory, "assumption1_filled.csv"))
write.csv(assumption2_filled, file.path(newDirectory, "assumption2_filled.csv"))
write.csv(ratio_duplicata, file.path(newDirectory, "ratio_duplicata.csv"))