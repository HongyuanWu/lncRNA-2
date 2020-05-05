
GetGeneList <- function(directory){
  FileList <- list.files(directory, pattern = ".txt")
  ResList <- data.frame()
  for (file in FileList){
    path <- paste0(directory,file)
    Filedata <- as.list(read.delim(path, sep = "\n"))
    ResList <- append(ResList, Filedata)
    

  }
  
  return(ResList)
}




