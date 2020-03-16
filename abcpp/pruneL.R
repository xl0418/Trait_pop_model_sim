pruneL=function (L,dropextinct=T)
{
  brts = NULL
  L = L[order(abs(L[, 3])), 1:4]
  age = L[1, 1]
  L[, 1] = age - L[, 1]
  L[1, 1] = -1
  notmin1 = which(L[, 4] != -1)
  L[notmin1, 4] = age - L[notmin1, 4]

  sall = which(L[, 4] == -1)
  tend = age

  L = L[, -4]
  linlist = cbind(data.frame(L[sall, ]))
  linlist_prun=linlist
  done = 0
  while (done == 0) {
    j = which.max(linlist[, 1])
    daughter = linlist[j, 3]
    parent = linlist[j, 2]
    parentj = which(parent == linlist[, 3])
    parentinlist = length(parentj)
    if (parentinlist == 1) {
      brts = c(brts, linlist[j, 1])
      linlist = linlist[-j, ]
    }
    else {
      linlist[j, 1:3] = L[which(L[, 3] == parent), 1:3]
      linlist_prun[which(linlist_prun[,3]==daughter), 1:3] = L[which(L[, 3] == parent), 1:3]
      linlist_prun[which(linlist_prun[,2]==daughter),2]=parent
    }
    if (nrow(linlist) == 1) {
      done = 1
    }
  }
  brts = rev(sort(age - brts))
  brts=c(age,brts)
  L_p=cbind(linlist_prun[order(linlist_prun[,1]),1:3],-1)
  L_p[,1]=brts
  L_temp=L_p
  species_index=c(1:nrow(L_p))
  L_p[,3]=species_index*(-1)^(L_p[,3]<0)
  for (i in 2:nrow(L_p)){
    pos=which(L_temp[,3] %in% L_temp[i,2])
    L_p[i,2]=L_p[pos,3]
  }
  L_p = as.matrix(L_p)
  dimnames(L_p) = NULL
  return(L_p)
}

