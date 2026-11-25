f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'));root<-release_script_root();a<-release_args()
if(is.null(a$analysis)||is.null(a$kind)||!a$kind %in% c('weights','cv5','calibration'))
  stop('Supply --analysis new_trial/analysis.rds --kind weights|cv5|calibration --seed N --out new_directory')
seed<-release_positive_integer(a$seed,'seed')
source(file.path(root,'R/load.R'));e<-release_load(root)
if(isTRUE(a[['dry-run']])) {cat('Entry parsed:',a$kind,'; no saved input read, RNG draw or fit.\n');quit(status=0L)}
x<-readRDS(a$analysis);e$release_focus_scope(x)
if(is.null(a$out)||file.exists(a$out))stop('Supply a new --out directory')
dir.create(a$out,recursive=TRUE);RNGkind("L'Ecuyer-CMRG",'Inversion','Rejection');set.seed(seed)
tryCatch(withCallingHandlers({
  ans<-switch(a$kind,weights=setNames(lapply(c('original','equal','normalized'),function(r)e$release_reweight(x,r)),c('original','equal','normalized')),
    cv5=e$release_cv5(x),calibration=e$release_calibration_sensitivity(x))
  saveRDS(ans,file.path(a$out,'sensitivity.rds'))
},warning=function(w)stop('Warning: ',conditionMessage(w))),error=function(er){
  saveRDS(list(error=conditionMessage(er),retry=FALSE),file.path(a$out,'failure.rds'));stop(er)
})
