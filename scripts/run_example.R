f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'));root<-release_script_root()
a<-release_args();out<-if(is.null(a$out))'results/example' else a$out
cmd<-c('--vanilla',file.path(root,'scripts/run_simulation.R'),'--scenario','Same_N200_r1_D0',
  '--reps','1','--seed','902104','--workers','1','--smoke','--out',out)
quit(status=system2(file.path(R.home('bin'),'Rscript'),shQuote(cmd)))
