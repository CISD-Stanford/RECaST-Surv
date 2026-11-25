f<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),fixed=TRUE)
root<-normalizePath(file.path(dirname(f),'..'));source(file.path(root,'R/load.R'));e<-release_load(root)
for(p in list.files(root,pattern='\\.R$',recursive=TRUE,full.names=TRUE))parse(p)
had<-exists('.Random.seed',.GlobalEnv,inherits=FALSE);prior<-if(had).Random.seed else NULL
r<-e$release_roster();stopifnot(length(r)==108L,!anyDuplicated(vapply(r,`[[`,'','id')))
stopifnot(identical(as.integer(table(vapply(r,`[[`,'','family'))),c(58L,14L,36L)))
main<-Filter(function(s)s$family=='main',r)
keys<-vapply(main,function(s)paste(s$spec$family,s$spec$trial_N,s$spec$r,s$spec$delta),'')
expected<-expand.grid(family=c('Same','Diff'),N=c(200,300),r=1:3,delta=c(0,-.25,-.5),stringsAsFactors=FALSE)
expected_keys<-vapply(seq_len(nrow(expected)),function(i)paste(expected$family[i],expected$N[i],expected$r[i],expected$delta[i]),'')
stopifnot(setequal(keys,expected_keys))
for(s in r) {
  ctx<-e$release_context(s);c<-ctx$config
  stopifnot(c$n_calib==100L,c$calib_steps==1000L,c$calib_burnin==200L,c$calib_pred_draws==50L,
    c$final_steps==5000L,c$final_burnin==1000L,c$n_thin_draws==1000L,identical(c$times_grid,1:60),
    all(is.finite(c$external_beta)),is.numeric(c$external_beta),is.function(ctx$engine$simulate_survival_data))
  if(s$family=='A-D'&&startsWith(s$id,'C'))stopifnot(identical(c$init_par,c(0,-10,0,-2)))
}
stopifnot(identical(e$release_scenario('H3_null')$canonical_id,'Same_N300_r2_D0'))
# Fixed-input score-set boundary cases from the approved A05 specification.
stopifnot(e$invert_cutoffs(0,5,0,2)$status=='empty')
c<-e$invert_cutoffs(0,5,-2,1);stopifnot(c$status=='nonempty',c$lower==0,c$upper==.8,c$zero_boundary)
stopifnot(e$construct_ci(1,2,rep(1,99),rep(2,99))$status=='unavailable')
stopifnot(is.na(e$release_reference(e$release_scenario('H3_delayed'),'Cox-A')$value))
stopifnot(abs(e$release_reference(e$release_scenario('H3_delayed'),'median')$value-.78987812521002365)<1e-15)
# Inclusive learned-cutoff ties and the strict original R-MAP boundary.
z<-e$release_crmap(c(rep(.5,949),rep(1,51)),c(.9,1))
stopifnot(z$cutoff==1,identical(z$evaluation_rejections,c(FALSE,TRUE)))
stopifnot(identical(had,exists('.Random.seed',.GlobalEnv,inherits=FALSE)),if(had)identical(prior,.Random.seed)else TRUE)
cat('PASS: 108 memberships, 105 distinct scenarios, production defaults, source configuration types, interval boundaries, H3 references, C-RMAP ties; no RNG or fit.\n')
