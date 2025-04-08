synthesise = function(method, encodedA, encodedB, PIVs, syntheticSample)
{
  
  if(method == "arf"){
    
    # missing values as NA not yet up-to-date
    # encodedBPIVsImputed = arf::impute(encodedB[,PIVs], m = 1)
    
    # missing values as value 0
    encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
    
    # synthesis
    arf = adversarial_rf(encodedB[,PIVs])
    psi = forde(arf, encodedB[,PIVs])
    syntheticNewB = forge(psi, syntheticSample)
    rownames(syntheticNewB) = 1:nrow(syntheticNewB)
    
    # encode missing as NA
    syntheticNewB[,PIVs][ syntheticNewB[,PIVs]==0 ] = NA

  }else if(method == "synthpop"){
    
    # missing values as NA, generates missing values
    
    for(p in 1:length(PIVs)){
      piv = PIVs[p]
      if(length(unique(encodedB[,piv]))>=60){
        encodedB[,piv] = as.numeric(encodedB[,piv])
      }
      if(length(unique(encodedA[,piv]))>=60){
        encodedA[,piv] = as.numeric(encodedA[,piv])
      }
    }
    
    syntheticNewB = syn(encodedB[,PIVs], k=syntheticSample)
    syntheticNewB = syntheticNewB$syn
    rownames(syntheticNewB) = 1:nrow(syntheticNewB)
    
  }else if(method == "mice"){
    
    # missing values as NA, does not generate missing values
    
    for(p in 1:length(PIVs)){
      piv = PIVs[p]
      if(length(unique(encodedB[,piv]))>=60){
        encodedB[,piv] = as.numeric(encodedB[,piv])
      }
      if(length(unique(encodedA[,piv]))>=60){
        encodedA[,piv] = as.numeric(encodedA[,piv])
      }
    }
    
    empty = matrix(NA, nrow=syntheticSample, ncol=length(PIVs)) 
    colnames(empty) = PIVs
    to_fill_synthetic = rbind(encodedB[,PIVs], empty)
    mice_synthesis = mice(to_fill_synthetic, 1)
    syntheticNewB = complete(mice_synthesis)
    rownames(syntheticNewB) = 1:nrow(syntheticNewB)
    syntheticNewB = syntheticNewB[(nrow(encodedB)+1):nrow(syntheticNewB),]
    
  }else{
    message("Only the following methods are supported: arf, synthpop and mice.")
  }
  
  columnsToAdd = colnames(encodedB)[(!colnames(encodedB) %in% PIVs) & (!colnames(encodedB) %in% c("localIndex","source"))]
  for(l in 1:length(columnsToAdd)){
    syntheticNewB[,columnsToAdd[l]] = rep(NA, nrow(syntheticNewB))
  }
  
  syntheticNewB$localIndex = 1:nrow(syntheticNewB)
  syntheticNewB$localIndex = syntheticNewB$localIndex + nrow(encodedB)
  syntheticNewB$source = "synthetic"
  
  encodedNewB = rbind(encodedB, syntheticNewB)
  
  for(p in 1:length(PIVs)){
    piv = PIVs[p]
    encodedNewB[,piv] = as.integer(encodedNewB[,piv])
    encodedA[,piv] = as.integer(encodedA[,piv])
  }
  
  list(NewA = encodedA, NewB = encodedNewB)
}