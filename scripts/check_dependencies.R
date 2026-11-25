f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'));release_require_dependencies()
cat(R.version.string,'\nCore dependencies are installed. No installation was attempted.\n')
for(p in c('survival','matrixStats','mvtnorm','statmod','MASS','jsonlite','digest'))cat(p,as.character(packageVersion(p)),'\n')
