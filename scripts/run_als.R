f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
source(file.path(dirname(f),'common.R'));root<-release_script_root();a<-release_args()
if(is.null(a$profile)||!a$profile %in% c('modern','historical')||is.null(a$sampling)||!a$sampling %in% c('random','biased')||
   is.null(a$scenario)||!a$scenario %in% c('NULL','ALTERNATIVE','DIAGNOSTIC'))
  stop('Supply --profile modern|historical --sampling random|biased --scenario NULL|ALTERNATIVE|DIAGNOSTIC --input supplied.csv --reps N --seed N --out new_directory')
reps<-release_positive_integer(a$reps,'reps');seed<-release_positive_integer(a$seed,'seed')
if(isTRUE(a[['dry-run']])) {cat('ALS entry parsed:',a$profile,a$sampling,a$scenario,reps,'identities; input excluded and not read.\n');quit(status=0L)}
release_require_dependencies()
if(is.null(a$input)||!file.exists(a$input))stop('Separately supplied author-approved synthetic CSV required; none is bundled')
expected<-'9b6e65cb614142d21509ffd5f9567931651284307ef9401769254f0dc1ce5563'
if(digest::digest(file=a$input,algo='sha256')!=expected)stop('Input differs from author-confirmed synthetic input; do not substitute another DGM')
if(is.null(a$out)||file.exists(a$out))stop('Supply new --out directory')
dir.create(a$out,recursive=TRUE);out<-normalizePath(a$out)
e<-new.env(parent=globalenv());e$Surv<-survival::Surv;e$coxph<-survival::coxph;e$basehaz<-survival::basehaz
modern<-file.path(root,'modules/ALS/modern')
for(p in c('config.R','R/als_core.R','R/reporting.R','smd_config.R','R/smd_reporting.R'))sys.source(file.path(modern,p),e)
sys.source(file.path(root,'modules/ALS/reporting_rules.R'),e)
df<-read.csv(a$input,stringsAsFactors=FALSE,check.names=FALSE)
RNGkind("L'Ecuyer-CMRG",'Inversion','Rejection');set.seed(seed)
context<-e$als_synthetic_context('separately_supplied_author_confirmed_synthetic_input')
if(a$profile=='modern') {
  prepared<-e$als_prepare_fixture(df,paste0('synthetic_row_',seq_len(nrow(df))),context)
  gen<-e$als_fit_generator_once(e$als_new_generator_state(),prepared,context)
  engine<-e$als_inherited_engine(modern,context)
} else {
  for(p in c('foreach','matrixStats'))if(!requireNamespace(p,quietly=TRUE))stop('Missing optional dependency: ',p)
  foreach::registerDoSEQ();e$foreach<-foreach::foreach;e$`%dopar%`<-foreach::`%dopar%`
  sys.source(file.path(modern,'R/inherited_engine_pinned.R'),e)
  sys.source(file.path(root,'modules/ALS/historical',paste0(a$sampling,'_profile.R')),e)
  covs<-c('SEX_NUM','AGE','BMI','FVC',if(a$sampling=='random')'HEIGHT','DISDUR','TOTAL','SLOPE')
  cache<-e$build_real_survival_cache(df,covs)
}
score_rows<-list()
for(i in seq_len(reps)) {
  leaf<-file.path(out,sprintf('rep_%04d',i));dir.create(leaf)
  tryCatch(withCallingHandlers({
    if(a$profile=='modern') {
      allocation<-e$als_sample_assignments(prepared,a$sampling,context)
      all<-e$als_paired_scenarios(prepared,allocation,gen,context)
      scenario<-all[[a$scenario]];X<-prepared$X[allocation$trial,,drop=FALSE]
      sources<-lapply(seq_along(allocation$external),function(j)e$als_fit_external(scenario$external[[j]],
        prepared$X[allocation$external[[j]],,drop=FALSE],context))
      ans<-e$als_recast_analysis(scenario$trial,X,lapply(sources,coef),engine$fit_target,engine$predict_draws,context)
      ans$reporting_rules<-if(a$scenario=='DIAGNOSTIC')e$als_release_rules(ans$final$oe$p_value,ans$split_p_values) else
        list(original_lower_tail=list(p=ans$final$oe$p_value,cutoff=ans$threshold,reject=ans$reject,
          uncalibrated_reject=ans$final$oe$p_value<.05))
      ans$cox<-lapply(c(FALSE,TRUE),function(ad)e$als_release_two_sided_cox(e$als_cox_fit(scenario$trial,X,ad,context)))
      ans$balance<-e$als_baseline_summary_v002(prepared,allocation,i,context)
      inputs<-list(trial=list(data=scenario$trial),external=lapply(scenario$external,function(z)list(data=z)))
    } else {
      extra<-if(a$sampling=='biased')list(trial_sampling='younger_shorter_disdur',trial_age_bias=.75,trial_disdur_bias=.75,trial_weight_cap=20)else list()
      tr<-do.call(e$simulate_survival_data,c(list(n=300,source='trial',r=2,cache=cache,spike_treatment=FALSE,exact_randomization=TRUE,total_time=60),extra))
      used<-tr$used_ids;external<-list()
      for(j in 1:3) {external[[j]]<-e$simulate_survival_data(c(500,1000,990)[j],source='external',cache=cache,used_ids=used,total_time=60);used<-external[[j]]$used_ids}
      if(a$scenario!='NULL')tr<-do.call(e$simulate_survival_data,c(list(n=300,source='trial',r=2,cache=cache,
        sampled_ids=tr$sampled_ids,fixed_trt=tr$trt,spike_treatment=TRUE,delta=if(a$scenario=='ALTERNATIVE')-.5 else 0,total_time=60),extra))
      sources<-e$fit_external_models(external)
      ans<-e$one_simulation(tr,sources,1:60,2)
      ans$profile_qualification<-'Historical-profile candidate; exact historical Table 3 producer is unauthenticated'
      inputs<-list(trial=tr,external=external)
    }
    if(!is.null(a[['model-file']])) {
      source(file.path(root,'R/load.R'));re<-release_load(root)
      ans$rmap<-re$release_rmap(inputs,a[['model-file']],as.integer((as.double(seed)+i)%%2147483646+1))
      score_rows[[i]]<-data.frame(bank=paste(a$profile,a$scenario,a$sampling,sep='_'),replication_id=i,
        score=ans$rmap$score,estimate=ans$rmap$estimate,reject_RMAP=ans$rmap$reject_RMAP,status='COMPLETED')
    }
    saveRDS(list(profile=a$profile,sampling=a$sampling,scenario=a$scenario,replication_id=i,result=ans),file.path(leaf,'analysis.rds'))
  },warning=function(w)stop('Scientific warning; no retry: ',conditionMessage(w))),error=function(er){
    saveRDS(list(error=conditionMessage(er),retry=FALSE),file.path(leaf,'failure.rds'));stop(er)
  })
}
if(length(score_rows))write.csv(do.call(rbind,score_rows),file.path(out,'scores.csv'),row.names=FALSE)
