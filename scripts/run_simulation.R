f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'))
root<-release_script_root();a<-release_args();release_require_dependencies()
source(file.path(root,'R/load.R'));e<-release_load(root);roster<-e$release_roster()
if(isTRUE(a$list)) {
  print(do.call(rbind,lapply(roster,function(s)data.frame(scenario=s$id,family=s$family,canonical=s$canonical_id))),row.names=FALSE)
  quit(status=0L)
}
if(!is.null(a$scenario) && !is.null(a$family))stop('Choose --scenario or --family')
if(!is.null(a$scenario))selected<-list(e$release_scenario(a$scenario)) else {
  if(is.null(a$family)||!a$family %in% c('main','H1-H5','A-D','all'))stop('Supply --scenario ID or --family main|H1-H5|A-D|all; --list lists IDs')
  selected<-Filter(function(s)a$family=='all'||s$family==a$family,roster)
}
reps<-release_positive_integer(a$reps,'reps');seed<-release_positive_integer(a$seed,'seed')
workers<-release_positive_integer(a$workers,'workers')
if(is.null(a$out)||!nzchar(a$out))stop('Supply --out directory')
if(isTRUE(a$smoke)&&length(selected)!=1L)stop('Smoke is one scenario only')
canonical<-sort(unique(vapply(selected,`[[`,'','canonical_id')))
universe<-sort(unique(vapply(roster,`[[`,'','canonical_id')))
plan<-data.frame(canonical_id=rep(canonical,each=reps),replication_id=rep(seq_len(reps),length(canonical)))
print(list(scenario_memberships=length(selected),distinct_scenarios=length(canonical),
  planned_replications=nrow(plan),workers=workers,smoke_only=isTRUE(a$smoke),generation_only=isTRUE(a[['generate-only']])))
if(isTRUE(a[['dry-run']]))quit(status=0L)
out<-path.expand(a$out)
if(file.exists(out))stop('Output exists; choose a new directory: ',out)
dir.create(out,recursive=TRUE);out<-normalizePath(out)
RNGkind("L'Ecuyer-CMRG",'Inversion','Rejection');set.seed(seed)
roots<-setNames(vector('list',length(universe)),universe);state<-.Random.seed
for(id in universe){roots[[id]]<-state;state<-parallel::nextRNGStream(state)}
jobs<-lapply(seq_len(nrow(plan)),function(i) {
  r<-plan[i,];state<-roots[[r$canonical_id]]
  if(r$replication_id>1L)for(j in seq_len(r$replication_id-1L))state<-parallel::nextRNGSubStream(state)
  list(scenario=e$release_scenario(r$canonical_id),state=state,replication_id=r$replication_id,
    out=file.path(out,r$canonical_id,sprintf('rep_%04d',r$replication_id)),
    smoke=isTRUE(a$smoke),generation_only=isTRUE(a[['generate-only']]))
})
saveRDS(list(seed=seed,RNGkind=RNGkind(),plan=plan,states=lapply(jobs,`[[`,'state'),
  memberships=selected,smoke_only=isTRUE(a$smoke)),file.path(out,'run_manifest.rds'))
write.csv(plan,file.path(out,'planned_trials.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'SESSION.txt'))
runjob<-function(job,root) {
  source(file.path(root,'R/load.R'));e<-release_load(root)
  do.call(e$release_run_one,job)
}
if(workers==1L)records<-lapply(jobs,runjob,root=root) else {
  cl<-parallel::makePSOCKcluster(workers)
  records<-tryCatch(parallel::parLapply(cl,jobs,runjob,root=root),finally=parallel::stopCluster(cl))
}
saveRDS(records,file.path(out,'execution_records.rds'))
bad<-vapply(records,function(r)r$status=='FAILED',logical(1))
cat('Saved',length(records),'identities;',sum(bad),'failed without retry. Output:',out,'\n')
if(any(bad))quit(status=1L)
