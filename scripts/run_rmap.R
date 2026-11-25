f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'));root<-release_script_root();a<-release_args()
source(file.path(root,'R/load.R'));e<-release_load(root)
reps<-release_positive_integer(a$reps,'reps');seed<-release_positive_integer(a$seed,'seed')
workers<-release_positive_integer(a$workers,'workers');stopifnot(reps<=1000L)
roster<-e$release_roster()
if(!is.null(a$scenario))selected<-list(e$release_scenario(a$scenario)) else {
  if(is.null(a$family)||!a$family %in% c('main','H1-H5','A-D','all'))stop('Supply --scenario or --family')
  selected<-Filter(function(s)a$family=='all'||s$family==a$family,roster)
}
key<-function(s)if(s$family=='A-D'&&startsWith(s$id,'C'))paste0('C_common_',if(s$delta==0)'D0' else 'Dm05')else s$canonical_id
keys<-sort(unique(vapply(selected,key,'')));universe<-sort(unique(vapply(roster,key,'')))
if(isTRUE(a[['dry-run']])) {cat(length(keys),'canonical R-MAP banks;',length(keys)*reps,'fits; no execution.\n');quit(status=0L)}
if(is.null(a[['model-file']]))stop('Supply authorized --model-file; the study JAGS model is not redistributed')
model<-normalizePath(a[['model-file']],mustWork=TRUE)
required_model<-jsonlite::read_json(file.path(root,'external/RMAP_MODEL_REQUIRED.json'))
stopifnot(digest::digest(file=model,algo='sha256')==required_model$sha256)
for(pkg in c('rjags','coda'))if(!requireNamespace(pkg,quietly=TRUE))stop('Missing optional dependency: ',pkg)
input_run<-if(is.null(a[['input-run']]))NULL else normalizePath(a[['input-run']],mustWork=TRUE)
if(!is.null(input_run)) {
  primary_manifest<-readRDS(file.path(input_run,'run_manifest.rds'))
  stopifnot(!isTRUE(primary_manifest$smoke_only))
}
if(is.null(a$out)||file.exists(a$out))stop('Supply new --out directory')
dir.create(a$out,recursive=TRUE);out<-normalizePath(a$out)
RNGkind("L'Ecuyer-CMRG",'Inversion','Rejection');set.seed(seed)
roots<-setNames(vector('list',length(universe)),universe);state<-.Random.seed
for(id in universe){roots[[id]]<-state;state<-parallel::nextRNGStream(state)}
jobs<-list()
for(k in keys) {
  s<-selected[[which(vapply(selected,key,'')==k)[1]]]
  if(startsWith(k,'C_common'))s<-e$release_scenario(if(s$delta==0)'C1_D0' else 'C1_Dm05')
  if(s$family=='H1-H5'&&s$canonical_id!=s$id)s<-e$release_scenario(s$canonical_id)
  state<-roots[[k]]
  for(i in seq_len(reps)) {
    jobs[[length(jobs)+1L]]<-list(scenario=s,bank=k,replication_id=i,state=state,
      jags_seed=as.integer((as.double(seed)+match(k,universe)*100000+i)%%2147483646+1))
    state<-parallel::nextRNGSubStream(state)
  }
}
saveRDS(list(seed=seed,jobs=jobs,memberships=selected,input_run=input_run),file.path(out,'rmap_manifest.rds'))
run<-function(job,root,out,model,input_run=NULL) {
  source(file.path(root,'R/load.R'));e<-release_load(root)
  leaf<-file.path(out,job$bank,sprintf('rep_%04d',job$replication_id));dir.create(leaf,recursive=TRUE)
  assign('.Random.seed',job$state,.GlobalEnv)
  ans<-tryCatch({
    if(is.null(input_run)) {ctx<-e$release_context(job$scenario);input<-e$release_generate(job$scenario,ctx)} else {
      source_file<-file.path(input_run,job$scenario$canonical_id,sprintf('rep_%04d',job$replication_id),'inputs.rds')
      input<-readRDS(source_file)
      saveRDS(list(path=source_file,sha256=digest::digest(file=source_file,algo='sha256')),file.path(leaf,'shared_input_identity.rds'))
    }
    saveRDS(input,file.path(leaf,'inputs.rds'));v<-e$release_rmap(input,model,job$jags_seed)
    saveRDS(v,file.path(leaf,'fit.rds'));data.frame(bank=job$bank,replication_id=job$replication_id,
      score=v$score,estimate=v$estimate,reject_RMAP=v$reject_RMAP,status='COMPLETED')
  },error=function(er){saveRDS(list(error=conditionMessage(er),retry=FALSE),file.path(leaf,'failure.rds'))
    data.frame(bank=job$bank,replication_id=job$replication_id,score=NA_real_,estimate=NA_real_,reject_RMAP=NA,status='FAILED')})
  ans
}
if(workers==1L)result<-lapply(jobs,run,root=root,out=out,model=model,input_run=input_run)else {
  cl<-parallel::makePSOCKcluster(workers)
  result<-tryCatch(parallel::parLapply(cl,jobs,run,root=root,out=out,model=model,input_run=input_run),finally=parallel::stopCluster(cl))
}
tab<-do.call(rbind,result);write.csv(tab,file.path(out,'scores.csv'),row.names=FALSE,na='NA')
if(any(tab$status=='FAILED'))quit(status=1L)
