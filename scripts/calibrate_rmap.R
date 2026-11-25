f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'));root<-release_script_root();a<-release_args()
source(file.path(root,'R/load.R'));e<-release_load(root)
if(isTRUE(a[['dry-run']])){cat('C-RMAP entry: exact 1000-score null bank, rank 950, inclusive >=; no scores read.\n');quit(status=0L)}
if(is.null(a$null)||is.null(a$evaluation)||is.null(a[['null-bank']])||is.null(a[['evaluation-bank']])||is.null(a$out))
  stop('Supply --null scores.csv --null-bank ID --evaluation scores.csv --evaluation-bank ID --out new_directory')
n<-read.csv(a$null,stringsAsFactors=FALSE);v<-read.csv(a$evaluation,stringsAsFactors=FALSE)
n<-n[n$bank==a[['null-bank']],];v<-v[v$bank==a[['evaluation-bank']],]
stopifnot(!anyDuplicated(n$replication_id),!anyDuplicated(v$replication_id),nrow(n)==1000L,
  all(n$status=='COMPLETED'),all(v$status=='COMPLETED'),nrow(v)>0L,!dir.exists(a$out))
ans<-e$release_crmap(n$score,v$score);dir.create(a$out,recursive=TRUE)
saveRDS(ans,file.path(a$out,'calibration.rds'))
v$cutoff<-ans$cutoff;v$reject_CRMAP<-ans$evaluation_rejections
write.csv(v,file.path(a$out,'decisions.csv'),row.names=FALSE)
cat('Cutoff:',sprintf('%.17g',ans$cutoff),'Evaluation rejections:',sum(ans$evaluation_rejections),'/',nrow(v),'\n')
