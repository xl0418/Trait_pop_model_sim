
process.data <- function(emdata.file, save.dir){
  source('phylo2L.R')
  source('pruneL.R')
  emdatadir = emdata.file  # NEXUS file
  dir = save.dir  # treedata directory
  
  emdata = read.nexus(emdatadir)
  
  
  # Pruned tree plot
  dropextinct = T
  species = phylo2L(emdata,error = 1e-5)
  L_ext = species$L
  extantspecieslabel = species$ESL

  L = L_ext
  
  prune = 1
  if (prune==1){
    L=pruneL(L)
  }
  
  # correct the crown age for baleen whales
  L[1,] = L[2,]
  L[1,2] = 0
  L[1,3] = -1
  L[2,3] = 2
  
  positive.clade = c(-2)
  do = TRUE
  while(do){
    negative.row = which(match(L[,2],positive.clade)>0)
    if(length(negative.row) == 0){
      break
    }
    L[negative.row,2] = - L[negative.row,2]
    positive.clade = L[negative.row,3]
    L[negative.row,3] = - L[negative.row,3]
    
  }
  
  time.list = c(sort(c(L[,1],L[which(L[,4]!= -1),4]),decreasing = TRUE),0)
  #total number of species
  num.species = nrow(L)
  trait.table = matrix(0,nrow = length(time.list)-1,ncol = nrow(L)+1)
  time.branching = match(L[,1],time.list)
  time.end = match(L[,4],time.list)
  time.end[is.na(time.end)] = length(time.list)
  timelist = as.data.frame(time.list)
  timebranch = as.data.frame(time.branching)
  timeend = as.data.frame(time.end)
  
  
  for(i in 1:num.species){
    
    trait.table[,i+1][time.branching[i]:(time.end[i]-1) ] = 1
  }
  trait.table = rbind(trait.table,trait.table[dim(trait.table)[1],])
  trait.table[,1] = time.list
  existing_species_table = trait.table[-1,-1]
  
  
  write.csv(timelist, file = paste0(dir,"/timelist.csv"))
  write.csv(timebranch, file = paste0(dir,"/timebranch.csv"))
  write.csv(timeend, file = paste0(dir,"/timeend.csv"))
  write.csv(existing_species_table, file = paste0(dir,"/traittable.csv"))
  write.csv(L, file = paste0(dir,"/Ltable.csv"))
  write.csv(extantspecieslabel, file = paste0(dir,"/extantspecieslabels.csv"))
}
