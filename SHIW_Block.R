################################################################################
#                RUN APPLICATIONS FOR THE RECORD LINKAGE TASK                  #
#          AND COMPARE FLEXRL WITH A NAIVE APPROACH ON THE FULL DATA           #
################################################################################

################################################################################
#                       CREATE A TABLE TO SAVE RESULTS                         #
################################################################################

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

estimated_fdr_realdata = data.frame(matrix(0,20,4))
colnames(estimated_fdr_realdata) = c("Naive","BRL","FastLink","FlexRL")

estimated_fdr_alldata = data.frame(matrix(0,20,4))
colnames(estimated_fdr_alldata) = c("Naive","BRL","FastLink","FlexRL")

assumption_condition_4 = data.frame(matrix(0,20,4))
colnames(assumption_condition_4) = c("Naive","BRL","FastLink","FlexRL")

ratio_duplicata = data.frame(matrix(0,20,1))
colnames(ratio_duplicata) = c("Naive")

################################################################################
#                               PREPARE THE DATA                               #
################################################################################

df1 = read.csv("SHIWData/df2016.csv", row.names = 1)
df2 = read.csv("SHIWData/df2020.csv", row.names = 1)

# Blocking = "None", "IREG", "AREA3"
if(Blocking == "None"){
  print("No blocking.")
  print("BlockingValue is set to NA")
  BlockingValue = "NA"
  
  PIVs_config = list( ANASCI     = list(stable = TRUE),
                      SESSO      = list(stable = TRUE),
                      STUDIO     = list(stable = TRUE),
                      STACIV     = list(stable = TRUE),
                      IREG       = list(stable = TRUE) )
}else if(Blocking == "IREG"){
  print("Blocking IREG.")
  print("IREG values:")
  print(unique(c(df1$IREG,df2$IREG)))
  print("Blocking on:")
  print(BlockingValue)
  
  df1 = df1[df1$IREG==as.integer(BlockingValue), ]
  df2 = df2[df2$IREG==as.integer(BlockingValue), ]
  
  PIVs_config = list( ANASCI     = list(stable = TRUE),
                      SESSO      = list(stable = TRUE),
                      STUDIO     = list(stable = TRUE),
                      STACIV     = list(stable = TRUE) )
}else if(Blocking == "AREA3"){
  print("Blocking AREA3.")
  print("AREA3 values:")
  print(unique(c(df1$AREA3,df2$AREA3)))
  print("Blocking on:")
  print(BlockingValue)
  
  df1 = df1[df1$AREA3==as.integer(BlockingValue), ]
  df2 = df2[df2$AREA3==as.integer(BlockingValue), ]
  
  PIVs_config = list( ANASCI     = list(stable = TRUE),
                      SESSO      = list(stable = TRUE),
                      STUDIO     = list(stable = TRUE),
                      STACIV     = list(stable = TRUE),
                      IREG       = list(stable = TRUE) )
}else{
  print("Blocking should take value in None, IREG, AREA3")
}

newDirectory = sprintf("Data:%s-Blocking:%s-Value:%s-Method:%s-%s", DataSource, Blocking, BlockingValue, Method, Sys.time())
dir.create(newDirectory)

PIVs = names(PIVs_config)
PIVs_stable = sapply(PIVs_config, function(x) x$stable)
controlOnMistakes = c(TRUE, TRUE, TRUE, TRUE, TRUE)

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

### CHECK DEPENDENCIES AMONG PIVs
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

### ENCODING
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

write.csv(encodedA, file.path(newDirectory, "dataset_A_used.csv"))
write.csv(encodedB, file.path(newDirectory, "dataset_B_used.csv"))

for(iterator in 1:10){
  
  ### SYNTHESISE
  
  subsample_size = as.integer(0.10*nbrRealRecordsB)
  
  Newdata = synthesise(Method, encodedA, encodedB, PIVs, subsample_size)
  NewA = Newdata$NewA
  NewB = Newdata$NewB
  
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
    synthfdrhat_onrealdata = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat_onalldata = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    assumption_condition4 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    # SAVE RESULTS
    thef1score[iterator,"BRL"] = f1score
    thesensitivity[iterator,"BRL"] = sensitivity
    thetp[iterator,"BRL"] = truepositive
    thefp[iterator,"BRL"] = falsepositive
    thefn[iterator,"BRL"] = falsenegative
    thesynth[iterator,"BRL"] = synthfalsepositive
    true_fdr[iterator,"BRL"] = 1-precision
    prob_fdr[iterator,"BRL"] = NA
    estimated_fdr_realdata[iterator,"BRL"] = synthfdrhat_onrealdata
    estimated_fdr_alldata[iterator,"BRL"] = synthfdrhat_onalldata
    assumption_condition_4[iterator,"BRL"] = abs(assumption_condition4)
    # SAVE THE INDICES FROM FILE A AND B OF THE RECORDS WHICH GET LINKED
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
    synthfdrhat_onrealdata = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat_onalldata = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    assumption_condition4 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    # SAVE RESULTS
    thef1score[iterator,"FastLink"] = f1score
    thesensitivity[iterator,"FastLink"] = sensitivity
    thetp[iterator,"FastLink"] = truepositive
    thefp[iterator,"FastLink"] = falsepositive
    thefn[iterator,"FastLink"] = falsenegative
    thesynth[iterator,"FastLink"] = synthfalsepositive
    true_fdr[iterator,"FastLink"] = 1-precision
    prob_fdr[iterator,"FastLink"] = fdrhat
    estimated_fdr_realdata[iterator,"FastLink"] = synthfdrhat_onrealdata
    estimated_fdr_alldata[iterator,"FastLink"] = synthfdrhat_onalldata
    assumption_condition_4[iterator,"FastLink"] = abs(assumption_condition4)
    # SAVE THE INDICES FROM FILE A AND B OF THE RECORDS WHICH GET LINKED
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
    synthfdrhat_onrealdata = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat_onalldata = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    assumption_condition4 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    # SAVE RESULTS
    thef1score[iterator,"Naive"] = f1score
    thesensitivity[iterator,"Naive"] = sensitivity
    thetp[iterator,"Naive"] = truepositive
    thefp[iterator,"Naive"] = falsepositive
    thefn[iterator,"Naive"] = falsenegative
    thesynth[iterator,"Naive"] = synthfalsepositive
    true_fdr[iterator,"Naive"] = 1-precision
    prob_fdr[iterator,"Naive"] = NA
    estimated_fdr_realdata[iterator,"Naive"] = synthfdrhat_onrealdata
    estimated_fdr_alldata[iterator,"Naive"] = synthfdrhat_onalldata
    assumption_condition_4[iterator,"Naive"] = abs(assumption_condition4)
    ratio_duplicata[iterator,"Naive"] = length(linked_real_pairs) / (nbrRealRecordsA*nbrRealRecordsB)
    # SAVE THE INDICES FROM FILE A AND B OF THE RECORDS WHICH GET LINKED
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
              StEMIter = 5,
              StEMBurnin = 2,
              GibbsIter = 5,
              GibbsBurnin = 2,
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
    synthfdrhat_onrealdata = (synthfalsepositive * (nbrRealRecordsB/subsample_size)) / (nrow(DeltaResult) - synthfalsepositive)
    synthfdrhat_onalldata = (synthfalsepositive * (1 + (nbrRealRecordsB/subsample_size))) / nrow(DeltaResult)
    assumption_condition4 = synthfalsepositive/(nbrRealRecordsA*subsample_size) - falsepositive/(nbrRealRecordsA*nbrRealRecordsB)
    # SAVE RESULTS
    thef1score[iterator,"FlexRL"] = f1score
    thesensitivity[iterator,"FlexRL"] = sensitivity
    thetp[iterator,"FlexRL"] = truepositive
    thefp[iterator,"FlexRL"] = falsepositive
    thefn[iterator,"FlexRL"] = falsenegative
    thesynth[iterator,"FlexRL"] = synthfalsepositive
    true_fdr[iterator,"FlexRL"] = 1-precision
    prob_fdr[iterator,"FlexRL"] = fdrhat
    estimated_fdr_realdata[iterator,"FlexRL"] = synthfdrhat_onrealdata
    estimated_fdr_alldata[iterator,"FlexRL"] = synthfdrhat_onalldata
    assumption_condition_4[iterator,"FlexRL"] = abs(assumption_condition4)
    # SAVE THE INDICES FROM FILE A AND B OF THE RECORDS WHICH GET LINKED
    write.csv(DeltaResult, file.path(newDirectory, sprintf("DeltaResult_FlexRL_%s_%s_%s_%s.csv", iterator, nbrRealRecordsA, nbrRealRecordsB, Sys.time())))
  }else{
    true_fdr[iterator,"FlexRL"] = "nolink"
  }
  print("FlexRL done")
  
  write.csv(thef1score, file.path(newDirectory, "thef1score.csv"))
  write.csv(thesensitivity, file.path(newDirectory, "thesensitivity.csv"))
  write.csv(thetp, file.path(newDirectory, "thetp.csv"))
  write.csv(thefp, file.path(newDirectory, "thefp.csv"))
  write.csv(thefn, file.path(newDirectory, "thefn.csv"))
  write.csv(thesynth, file.path(newDirectory, "thesynth.csv"))
  write.csv(true_fdr, file.path(newDirectory, "true_fdr.csv"))
  write.csv(prob_fdr, file.path(newDirectory, "prob_fdr.csv"))
  write.csv(estimated_fdr_realdata, file.path(newDirectory, "estimated_fdr_realdata.csv"))
  write.csv(estimated_fdr_alldata, file.path(newDirectory, "estimated_fdr_alldata.csv"))
  write.csv(assumption_condition_4, file.path(newDirectory, "assumption_condition_4.csv"))
  write.csv(ratio_duplicata, file.path(newDirectory, "ratio_duplicata.csv"))
  
}