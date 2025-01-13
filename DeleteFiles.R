DeleteFiles<-function(dir, patterns) {
  files_to_delete = c()
  for (p in patterns) {
    files_to_delete<-c(files_to_delete, list.files(dir, pattern=p, full.names=TRUE))
  }
  for (f in files_to_delete) {
    file.remove(f)
  }
}