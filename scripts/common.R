release_script_root <- function() {
  f<-grep('^--file=',commandArgs(FALSE),value=TRUE)
  if(length(f)!=1L)stop('Run this entry with Rscript --vanilla scripts/<name>.R')
  normalizePath(file.path(dirname(gsub('~+~',' ',sub('^--file=','',f),fixed=TRUE)),'..'),mustWork=TRUE)
}
release_args <- function(args=commandArgs(TRUE)) {
  z<-list();i<-1L
  flags<-c('smoke','dry-run','generate-only','list')
  while(i<=length(args)) {
    key<-sub('^--','',args[i]);if(key==args[i])stop('Expected --option: ',args[i])
    if(key %in% flags){z[[key]]<-TRUE;i<-i+1L} else {
      if(i==length(args))stop('Missing value for ',args[i]);z[[key]]<-args[i+1L];i<-i+2L
    }
  };z
}
release_positive_integer <- function(x,name) {
  v<-suppressWarnings(as.numeric(x))
  if(length(v)!=1L||!is.finite(v)||v<1||v>2147483646||v!=floor(v))stop('Supply a positive integer --',name)
  as.integer(v)
}
release_require_dependencies <- function() {
  for(p in c('survival','matrixStats','mvtnorm','statmod','MASS','jsonlite','digest'))
    if(!requireNamespace(p,quietly=TRUE))stop('Missing installed dependency: ',p,'; see README.md')
}
