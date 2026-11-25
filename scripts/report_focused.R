f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'));root<-release_script_root();a<-release_args()
if(is.null(a$kind)||!a$kind %in% c('full-vs-split','coefficients'))stop('Supply --kind full-vs-split|coefficients --input RUN_DIRECTORY --out NEW_DIRECTORY')
if(isTRUE(a[['dry-run']])){cat('Saved-output focused reporting entry parsed:',a$kind,'; no records read or calculations run.\n');quit(status=0L)}
source(file.path(root,'R/load.R'));e<-release_load(root)
sys.source(file.path(root,'R/reporting/focused.R'),e)
e$release_report_focused(a$input,a$out,a$kind)
