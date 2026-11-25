release_roster <- function() {
  m <- jsonlite::read_json(file.path(release_root,'config/main/scenarios.json'),simplifyVector=FALSE)
  main <- lapply(m,function(s)list(id=s$scenario_id,family='main',canonical_id=s$scenario_id,spec=s))
  h <- read.csv(file.path(release_root,'config/H1-H5/roster.csv'),stringsAsFactors=FALSE)
  hh <- lapply(seq_len(nrow(h)),function(i) {
    r <- h[i,]; canonical <- r$cell_id
    if(r$cell_id %in% c('H1_m1_D0','H3_null'))canonical <- 'Same_N300_r2_D0'
    if(r$cell_id=='H1_m1_Dm05')canonical <- 'Same_N300_r2_Dm05'
    list(id=r$cell_id,family='H1-H5',canonical_id=canonical,spec=as.list(r))
  })
  a <- read.csv(file.path(release_root,'config/A-D/cases.csv'),stringsAsFactors=FALSE)
  aa <- unlist(lapply(seq_len(nrow(a)),function(i)lapply(c(0,-.5),function(d) {
    id <- paste0(a$case[i],'_',if(d==0)'D0' else 'Dm05')
    list(id=id,family='A-D',canonical_id=id,spec=as.list(a[i,]),delta=d)
  })),recursive=FALSE)
  c(main,hh,aa)
}

release_scenario <- function(id) {
  all <- release_roster(); matches <- which(vapply(all,function(x)x$id==id,logical(1)))
  if(length(matches)!=1L)stop('Unknown or duplicate scenario: ',id)
  all[[matches]]
}

release_context <- function(scenario,smoke=FALSE) {
  cfg <- A03a_config(); s <- scenario$spec
  engine <- new.env(parent=environment(release_context))
  for(n in c('simulate_survival_data','fit_external_models','compare_predicted_to_observed',
             'mcmc_survival_ph','posterior_predict_srv_prob_ph_optimized'))engine[[n]] <- get(n)
  f <- environment(A03a_stage)
  if(scenario$family=='main') {
    for(n in c('trial_N','r','delta','lambda','nu','total_time','dropout','family'))cfg[[n]] <- s[[n]]
    cfg$trial_beta <- as.numeric(unlist(s$trial_beta))
    cfg$external_sizes <- as.integer(unlist(s$external_sizes))
    cfg$external_beta <- do.call(rbind,lapply(s$external_beta,function(z)as.numeric(unlist(z))))
    cfg$omitted_external_beta <- as.numeric(unlist(s$omitted_external_beta))
  } else if(scenario$family=='H1-H5') {
    configs <- dget(file.path(release_root,'config/H1-H5',paste0(s$scenario,'.R')))
    h <- configs[[match(s$variant,vapply(configs,function(z)z$variant,character(1)))]]
    H_validate_config(h);cfg$h_dgm <- h;cfg$family <- 'Same';cfg$delta <- h$effect$Delta
    cfg$trial_beta <- h$populations$target$beta
    cfg$external_beta <- do.call(rbind,lapply(h$populations[c('S1','S2','S3')],function(z)z$beta))
  } else {
    cfg$family <- 'Same';cfg$delta <- scenario$delta
    cfg$trial_N <- as.integer(s$trial_N);cfg$r <- 2
    cfg$trial_beta <- jsonlite::fromJSON(s$trial_beta)
    cfg$external_sizes <- as.integer(jsonlite::fromJSON(s$external_N))
    cfg$external_beta <- jsonlite::fromJSON(s$external_beta)
    cfg$AD <- list(case=s$case,
      rho=as.double(unlist(s[c('trial_rho_X1_X2',paste0('external_',1:3,'_rho_X1_X2'))])),
      support=as.character(unlist(s[c('trial_X1_support',paste0('external_',1:3,'_X1_support'))])),
      corruption=jsonlite::fromJSON(s$historical_corruption_map,simplifyVector=FALSE))
    if(startsWith(s$case,'C'))cfg$init_par <- c(0,-10,0,-2)
    gen <- new.env(parent=environment(release_context));gen$mvrnorm <- MASS::mvrnorm
    sys.source(file.path(release_root,'R/simulation',if(startsWith(s$case,'D'))'ad_restricted.R' else 'ad_general.R'),gen)
    engine$simulate_survival_data <- gen[[if(startsWith(s$case,'D'))'simulate_survival_data_biased' else 'simulate_survival_data']]
    f <- new.env(parent=environment(A03a_stage))
    sys.source(file.path(release_root,'R/analysis/median_prediction.R'),f)
  }
  cfg$smoke_only <- smoke
  if(smoke) {
    z <- jsonlite::read_json(file.path(release_root,'examples/smoke.json'),simplifyVector=TRUE)
    for(n in names(z$computation))cfg[[n]] <- as.integer(z$computation[[n]])
  }
  list(config=cfg,engine=engine,functions=f)
}

# Derived from the accepted H runtime helper: only its consumed HPC capability
# gate is omitted. Draw order, populations, treatment and survival arithmetic agree.
release_h_population <- function(n,source,hcfg,population) {
  p<-H_population(hcfg,population);stopifnot(n==p$n)
  X1<-rnorm(n);X2<-rnorm(n);X3<-rnorm(n);X<-cbind(X1,X2,X3)
  A<-if(source=='trial')rbinom(n,1,hcfg$treatment_probability)else rep(0,n)
  U<-runif(n);T<-H_event_time_from_uniform(hcfg,population,U,X,A)
  D<-rexp(n,rate=H_dropout_rate(hcfg,population));obs<-H_observation(hcfg,population,T,D)
  df<-cbind(data.frame(id=seq_len(n),source=source,trt=A,time=obs$time,status=obs$status),as.data.frame(X))
  list(data=df,X_mat=X,trt=A,time=obs$time,status=obs$status,t_time=T)
}

release_generate <- function(scenario,ctx) {
  cfg<-ctx$config;engine<-ctx$engine
  if(scenario$family=='A-D')return(AD_generate(engine,cfg))
  if(scenario$family=='H1-H5')return(list(
    trial=release_h_population(cfg$trial_N,'trial',cfg$h_dgm,'target'),
    external=lapply(1:3,function(j)release_h_population(cfg$external_sizes[j],'external',cfg$h_dgm,paste0('S',j)))))
  rate <- -log(1-cfg$dropout)/cfg$total_time
  trial <- engine$simulate_survival_data(cfg$trial_N,source='trial',betas=cfg$trial_beta,
    lambda=cfg$lambda,nu=cfg$nu,r=cfg$r,delta=cfg$delta,total_time=cfg$total_time,random_cens_rate=rate)
  external <- lapply(1:3,function(j)engine$simulate_survival_data(cfg$external_sizes[j],
    source='external',ext_model_type=if(cfg$family=='Same')'same' else 'different',
    betas=if(cfg$family=='Same')cfg$external_beta[j,] else c(cfg$external_beta[j,],cfg$omitted_external_beta[j]),
    lambda=cfg$lambda,nu=cfg$nu,total_time=cfg$total_time,random_cens_rate=rate))
  list(trial=trial,external=external)
}
