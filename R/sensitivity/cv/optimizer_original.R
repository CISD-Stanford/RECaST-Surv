# Selectively extracted definitions; no historical runner is sourced.
# See docs/PROVENANCE.csv for source identity and extraction scope.
lse <- function(x) {a<-max(x);if(a==-Inf)return(-Inf);a+log(sum(exp(x-a)))}

stack_optimize <- function(logcontrib) {
 stopifnot(is.matrix(logcontrib),ncol(logcontrib)==3,all(is.finite(logcontrib)))
 shift<-apply(logcontrib,1,max); A<-exp(logcontrib-shift); n<-nrow(A)
 if(any(A==0))stop('Optimizer row-scaled source contrast underflow; no silent zero, floor or fallback')
 obj<-function(w) {v<-drop(A%*%w);if(any(v<=0))return(-Inf);sum(shift+log(v))}
 grad<-function(w) colSums(A/drop(A%*%w))
 candidates<-diag(3); labels<-paste0('vertex_',1:3)
 for(a in 1:2)for(b in (a+1):3) {
  makew<-function(t){w<-numeric(3);w[a]<-t;w[b]<-1-t;w}
  z<-optimize(function(t)-obj(makew(t)),c(0,1),tol=1e-12)
  candidates<-rbind(candidates,makew(z$minimum));labels<-c(labels,paste0('edge_',a,b))
 }
 fn<-function(v)-obj(c(v,1-sum(v)))
 gr<-function(v){g<-grad(c(v,1-sum(v)));-(g[1:2]-g[3])}
 z<-constrOptim(c(1/3,1/3),fn,gr,ui=rbind(diag(2),c(-1,-1)),ci=c(0,0,-1),
  control=list(reltol=1e-12,maxit=5000),outer.iterations=100,outer.eps=1e-12)
 candidates<-rbind(candidates,c(z$par,1-sum(z$par)));labels<-c(labels,'interior')
 values<-apply(candidates,1,obj); best<-which.max(values);w<-candidates[best,]
 direct<-sum(apply(sweep(logcontrib,2,log(w),'+'),1,lse))
 g<-grad(w); active<-w>1e-8
 kkt<-max(c(abs(g[active]-n),pmax(g[!active]-n,0),0))
 D<-sweep(A[,1:2,drop=FALSE],1,A[,3],'-')/drop(A%*%w)
 eig<-eigen(crossprod(D),symmetric=TRUE,only.values=TRUE)$values
 globally_strict<-min(eig)>1e-10*max(1,max(eig))
 ai<-which(active)
 face_strict<-if(length(ai)==1) TRUE else {
  F<-sweep(A[,ai[-length(ai)],drop=FALSE],1,A[,ai[length(ai)]],'-')/drop(A%*%w)
  fe<-eigen(crossprod(F),symmetric=TRUE,only.values=TRUE)$values
  min(fe)>1e-10*max(1,max(fe))
 }
 # Strict inactive KKT plus strict concavity on the active face certifies a
 # unique boundary optimum even when global affine curvature is singular.
 identified<-globally_strict||(face_strict&&all(g[!active]<n-1e-7))
 near<-which(values>=max(values)-1e-9)
 separated<-any(vapply(near,function(j)max(abs(candidates[j,]-w))>1e-4,logical(1)))
 status<-if(!identified||separated)'indeterminate_or_nonunique' else if(kkt>1e-5||
  abs(direct-values[best])>1e-9||min(w)<0||abs(sum(w)-1)>1e-12)'numerical_failure' else 'unique_numerical_solution'
 list(status=status,weights=if(status=='unique_numerical_solution')w else rep(NA_real_,3),
  diagnostic_candidate=w,objective=direct,kkt_residual=kkt,eigenvalues=eig,
  simplex_sum=sum(w),simplex_min=min(w),objective_accuracy=abs(direct-values[best]),
  candidates=data.frame(label=labels,candidates,objective=values),solver=z,
  qualification='Numerical identification only, not a proof of posterior or inferential validity')
}

