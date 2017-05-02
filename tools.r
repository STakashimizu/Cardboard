parse_crv <- function(filename) {
  d_tmp = read.csv(filename, skip = 1) # ignore the 1st row
  rownames(d_tmp) = NULL
  colnames(d_tmp) = NULL
  
  # Change data structure for convenience
  num_jikken = (length(d_tmp)-1)/3
  d0 = matrix(NA, nrow = 0, ncol = 4)
  ds = list()
  for(jikken_no in 1:num_jikken) {
    col = (jikken_no - 1)*3 + 2
    tmp = d_tmp[,col:(col+2)]
    d0 = rbind(d0, cbind(jikken_no, tmp))
    ds[[jikken_no]] = tmp
  }
  d0 = d0[!is.na(d0[,2]),]
  colnames(d0) = c("No", "Hennikei", "Kajuu", "Idouryou") # where "No" is index of a experiment
  d = d0
  return(d)
}

part = function(result, maxNo = Inf) {
  maxNo = if(maxNo == Inf) max(result$No) else maxNo
  return(result[1:max(which(result$No == maxNo)),])
}

crop = function(result, by = "Idouryou", lower = 0, upper = Inf) {
  return(result[lower < result[by] & result[by] < upper,])
}

# From: http://stackoverflow.com/questions/29214932/split-a-file-path-into-folder-names-vector
split_path <- function(path, mustWork = FALSE, rev = FALSE) {
  output <- c(strsplit(dirname(normalizePath(path,mustWork = mustWork)),
                       "/|\\\\")[[1]], basename(path))
  ifelse(rev, return(rev(output)), return(output))
}