if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
project_root <- getwd()
for (file in c("00_packages.R", "01_config.R", "02_data.R", "03_model.R"))
  source(file, encoding = "UTF-8")
root <- Sys.getenv("PROJECT_ROOT", normalizePath("..", winslash = "/"))
source(file.path(root, "R", "project.R"), encoding = "UTF-8")
cfg <- read_project_config(root, Sys.getenv("PROJECT_CONFIG_FILE", "local_config.R"))$network
fit <- load_network_model(model_file)
validate_network_model(fit, expected_data = read_network_data())
out <- file.path(output_dir, "predictive_checks")
for (directory in c(out, file.path(out, "cache"), file.path(out, "tables"), file.path(out, "logs")))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
fingerprint <- unname(tools::md5sum(model_file))
nodes <- network_nodes
resp <- gsub("_", "", paste0(nodes, "_t2"))
if (!identical(fit$formula$responses, resp)) stop("Unexpected Gaussian response order.")
write_table <- function(x, name) write_csv(x, file.path(out, "tables", name))
log <- function(...) message(...)

generate_joint_ppc <- function() {
Y <- as.matrix(fit$data[paste0(nodes, "_t2")]); X <- as.matrix(fit$data[paste0(nodes, "_t1")])
colnames(Y) <- colnames(X) <- nodes
n <- nrow(Y); J <- ncol(Y)
lims <- t(vapply(nodes, function(z) if (grepl("^gad", z)) c(1,4) else if (grepl("^mpiq", z)) c(1,7) else c(0,15), numeric(2)))
groups <- list(Overall = seq_len(n), NF = which(fit$data$addiction == "Non-Addiction" & fit$data$gender == "Female"),
               AF = which(fit$data$addiction == "Addiction" & fit$data$gender == "Female"),
               NM = which(fit$data$addiction == "Non-Addiction" & fit$data$gender == "Male"),
               AM = which(fit$data$addiction == "Addiction" & fit$data$gender == "Male"))
G <- length(groups); R <- min(cfg$gaussian_ppc_draws, posterior::ndraws(posterior::as_draws_array(fit)))
set.seed(cfg$ppc_seed)
draw_dimensions <- dim(posterior::as_draws_array(fit))
counts <- rep(R %/% draw_dimensions[2], draw_dimensions[2])
counts[seq_len(R %% draw_dimensions[2])] <- counts[seq_len(R %% draw_dimensions[2])] + 1L
sel <- do.call(rbind, lapply(seq_len(draw_dimensions[2]), function(ch)
  data.frame(chain = ch, iteration = sort(sample.int(draw_dimensions[1], counts[ch])))))
sel$draw <- (sel$chain - 1L) * draw_dimensions[1] + sel$iteration
sel$selection_order <- seq_len(R)
sel$seed <- cfg$ppc_seed
sel$model_md5 <- fingerprint
write_table(sel, "ppc_draw_selection.csv")
cache_manifest_path <- file.path(out,"cache/ppc_cache_manifest.json")
first_pass_path <- file.path(out,"cache/ppc_first_pass.rds")
expected_batches <- vapply(seq(1L,R,by=25L),function(start) file.path(out,"cache",sprintf("ppc_batch_%04d_%04d.rds",start,min(R,start+24L))),character(1))
if (file.exists(first_pass_path) && file.exists(cache_manifest_path) && all(file.exists(expected_batches))) {
  cm <- jsonlite::read_json(cache_manifest_path,simplifyVector=TRUE)
  if (identical(cm$model_md5,fingerprint) && identical(cm$selection_md5,unname(tools::md5sum(file.path(out,"tables/ppc_draw_selection.csv")))) &&
      identical(cm$brms_version,as.character(packageVersion("brms"))) &&
      identical(cm$posterior_version,as.character(packageVersion("posterior"))) && identical(cm$specification,paste0("joint_ppc_v2_n", length(nodes), "_draws", R, "_seed", cfg$ppc_seed))) {
    log("Validated all PPC caches; generation skipped")
    return(invisible(NULL))
  }
}
sigmas <- posterior::as_draws_matrix(fit, variable = paste0("sigma_", resp))[sel$draw, , drop = FALSE]
sigmas <- sigmas[, paste0("sigma_", resp), drop = FALSE]
stopifnot(identical(colnames(sigmas), paste0("sigma_", resp)))
stat_names <- c("mean", "sd", "skewness", "q05", "q25", "q50", "q75", "q95",
                "below_L", "above_U", "outside", "near_L", "near_U", "analytic_outside", "epred_outside")
key <- function(g, j, s) paste(g, nodes[j], s, sep = "|")
node_rep <- new.env(hash = TRUE, parent = emptyenv())
node_obs <- list()
skew <- function(z) { a <- mean(z); v <- mean((z-a)^2); if (v <= 0) return(NA_real_); mean((z-a)^3) / v^(3/2) }
node_stats <- function(z, L, U) {
  q <- unname(quantile(z, c(.05,.25,.5,.75,.95), names = FALSE))
  c(mean=mean(z), sd=sd(z), skewness=skew(z), setNames(q,c("q05","q25","q50","q75","q95")),
    below_L=mean(z<L),above_U=mean(z>U),outside=mean(z<L | z>U),near_L=mean(z<L+.5),near_U=mean(z>U-.5))
}
for (g in names(groups)) for (j in seq_len(J)) {
  st <- node_stats(Y[groups[[g]],j], lims[j,1], lims[j,2])
  for (s in names(st)) node_obs[[key(g,j,s)]] <- st[s]
}
make_pair <- function(M) {
  ij <- which(upper.tri(M), arr.ind = TRUE)
  cbind(source=colnames(M)[ij[,1]], target=colnames(M)[ij[,2]])
}
follow_pairs <- make_pair(cor(Y))
cross_pairs <- as.matrix(expand.grid(source=nodes,target=nodes,stringsAsFactors=FALSE))
follow_rep <- array(NA_real_, c(R,G,nrow(follow_pairs)), dimnames=list(NULL,names(groups),NULL))
cross_rep <- array(NA_real_, c(R,G,nrow(cross_pairs)), dimnames=list(NULL,names(groups),NULL))
follow_obs <- cross_obs <- list()
domains <- list(Support=1:3, Anxiety=4:10, Smartphone_addiction=11:18)
domain_labels <- c("within_support_signed", "within_anxiety_signed", "within_addiction_signed",
                   "support_anxiety_absolute", "support_addiction_absolute", "anxiety_addiction_absolute")
domain_stat <- function(C) c(vapply(domains, function(ids) mean(C[ids,ids][upper.tri(C[ids,ids])]), numeric(1)),
  mean(abs(C[domains[[1]],domains[[2]]])), mean(abs(C[domains[[1]],domains[[3]]])), mean(abs(C[domains[[2]],domains[[3]]])))
domain_rep <- array(NA_real_,c(R,G,6),dimnames=list(NULL,names(groups),domain_labels))
domain_obs <- array(NA_real_,c(G,6),dimnames=list(names(groups),domain_labels))
for (g in names(groups)) {
  ids <- groups[[g]]; fc <- cor(Y[ids,,drop=FALSE]); xc <- cor(X[ids,,drop=FALSE],Y[ids,,drop=FALSE])
  follow_obs[[g]] <- fc[cbind(match(follow_pairs[,1],nodes),match(follow_pairs[,2],nodes))]
  cross_obs[[g]] <- xc[cbind(match(cross_pairs[,1],nodes),match(cross_pairs[,2],nodes))]
  domain_obs[g,] <- domain_stat(fc)
}
resid_sd_obs <- resid_sd_rep <- array(NA_real_,c(R,4,J),dimnames=list(NULL,names(groups)[-1],nodes))
resid_cor_obs <- resid_cor_rep <- array(NA_real_,c(R,4,nrow(follow_pairs)),dimnames=list(NULL,names(groups)[-1],NULL))
mu_sum <- matrix(0,n,J,dimnames=list(NULL,nodes))
cache_batch <- function(start,end,p,m) {
  path <- file.path(out,"cache",sprintf("ppc_batch_%04d_%04d.rds",start,end))
  saveRDS(list(p=p,m=m,draw=sel$draw[start:end],model_md5=fingerprint,seed=cfg$ppc_seed),path,compress="xz")
}
joint_check <- NULL
batch_size <- 25L
for (start in seq(1L,R,by=batch_size)) {
  end <- min(R,start+batch_size-1L); ids <- start:end
  set.seed(cfg$ppc_seed + start)
  p <- brms::posterior_predict(fit,draw_ids=sel$draw[ids],resp=NULL)
  m <- brms::posterior_epred(fit,draw_ids=sel$draw[ids],resp=NULL)
  stopifnot(identical(dim(p),c(length(ids),n,J)),identical(dim(m),c(length(ids),n,J)))
  if (!is.null(dimnames(p)[[3]])) stopifnot(identical(dimnames(p)[[3]],resp))
  if (!is.null(dimnames(m)[[3]])) stopifnot(identical(dimnames(m)[[3]],resp))
  if (start == 1L) {
    resid <- p[1,,]-m[1,,]
    rc <- cor(resid)
    draws <- posterior::as_draws_matrix(fit,variable=grep("^rescor__",posterior::variables(posterior::as_draws_array(fit)),value=TRUE))
    row <- draws[sel$draw[1],]
    names(row) <- colnames(draws)
    target <- diag(J); dimnames(target) <- list(resp,resp)
    for (nm in names(row)) {
      a <- strsplit(sub("^rescor__","",nm),"__",fixed=TRUE)[[1]]
      target[a[1],a[2]] <- target[a[2],a[1]] <- row[nm]
    }
    ij <- which(upper.tri(target), arr.ind=TRUE)
    joint_check <- data.frame(response1=resp[ij[,1]],response2=resp[ij[,2]],
      model=target[ij],empirical=rc[ij],absolute_difference=abs(target[ij]-rc[ij]))
    write_table(joint_check,"joint_residual_simulation_check.csv")
    log("Joint residual correlation check max abs delta: ",max(abs(joint_check$empirical-joint_check$model)))
    stopifnot(nrow(joint_check) == choose(J, 2))
    rm(draws)
  }
  cache_batch(start,end,p,m)
  mu_sum <- mu_sum + apply(m,c(2,3),sum)
  for (b in seq_along(ids)) {
    r <- ids[b]; yy <- p[b,,]; mm <- m[b,,]; eo <- Y-mm; er <- yy-mm
    for (g in names(groups)) {
      ix <- groups[[g]]
      fc <- cor(yy[ix,,drop=FALSE]); xc <- cor(X[ix,,drop=FALSE],yy[ix,,drop=FALSE])
      follow_rep[r,g,] <- fc[cbind(match(follow_pairs[,1],nodes),match(follow_pairs[,2],nodes))]
      cross_rep[r,g,] <- xc[cbind(match(cross_pairs[,1],nodes),match(cross_pairs[,2],nodes))]
      domain_rep[r,g,] <- domain_stat(fc)
      for (j in seq_len(J)) {
        st <- node_stats(yy[ix,j],lims[j,1],lims[j,2])
        st <- c(st,analytic_outside=mean(pnorm((lims[j,1]-mm[ix,j])/sigmas[r,j]) +
                        pnorm((lims[j,2]-mm[ix,j])/sigmas[r,j],lower.tail=FALSE)),
                epred_outside=mean(mm[ix,j]<lims[j,1] | mm[ix,j]>lims[j,2]))
        for (s in names(st)) {
          ky <- key(g,j,s)
          if (!exists(ky,envir=node_rep,inherits=FALSE)) assign(ky,rep(NA_real_,R),envir=node_rep)
          z <- get(ky,envir=node_rep,inherits=FALSE); z[r] <- st[s]; assign(ky,z,envir=node_rep)
        }
      }
      if (g != "Overall") {
        resid_sd_obs[r,g,] <- apply(eo[ix,,drop=FALSE],2,sd)
        resid_sd_rep[r,g,] <- apply(er[ix,,drop=FALSE],2,sd)
        co <- cor(eo[ix,,drop=FALSE]); cr <- cor(er[ix,,drop=FALSE])
        pairids <- cbind(match(follow_pairs[,1],nodes),match(follow_pairs[,2],nodes))
        resid_cor_obs[r,g,] <- co[pairids]; resid_cor_rep[r,g,] <- cr[pairids]
      }
    }
  }
  log("PPC generated draws ",start,"-",end)
  rm(p,m); gc(FALSE)
}
saveRDS(list(node_rep=as.list(node_rep),node_obs=node_obs,follow_rep=follow_rep,follow_obs=follow_obs,
  cross_rep=cross_rep,cross_obs=cross_obs,domain_rep=domain_rep,domain_obs=domain_obs,
  resid_sd_obs=resid_sd_obs,resid_sd_rep=resid_sd_rep,resid_cor_obs=resid_cor_obs,resid_cor_rep=resid_cor_rep,
  mu_bar=mu_sum/R, groups=groups, nodes=nodes, lims=lims, follow_pairs=follow_pairs,cross_pairs=cross_pairs,
  selection=sel, model_md5=fingerprint),file.path(out,"cache/ppc_first_pass.rds"),compress="xz")
jsonlite::write_json(list(model_md5=fingerprint,selection_md5=unname(tools::md5sum(file.path(out,"tables/ppc_draw_selection.csv"))),
  brms_version=as.character(packageVersion("brms")),posterior_version=as.character(packageVersion("posterior")),
  specification=paste0("joint_ppc_v2_n", length(nodes), "_draws", R, "_seed", cfg$ppc_seed),seed=cfg$ppc_seed,draws=R,batches=length(expected_batches)),
  cache_manifest_path,auto_unbox=TRUE,pretty=TRUE)
log("PPC first pass completed")

}
generate_joint_ppc()
a <- readRDS(file.path(out,"cache/ppc_first_pass.rds"))
if (!identical(a$model_md5,fingerprint)) stop("Predictive cache model mismatch.")
nodes <- a$nodes; J <- length(nodes); groups <- a$groups[c("NF","AF","NM","AM","Overall")]; G <- length(groups); R <- nrow(a$selection)
Y <- as.matrix(fit$data[paste0(nodes,"_t2")]); X <- as.matrix(fit$data[paste0(nodes,"_t1")]); colnames(Y) <- colnames(X) <- nodes
paired <- function(obs,rep) {
  if (length(obs)==1L) obs <- base::rep(obs,length(rep))
  valid <- is.finite(obs)&is.finite(rep)
  if (!all(valid)) { obs<-obs[valid]; rep<-rep[valid] }
  if (!length(rep)) return(setNames(base::rep(NA_real_, 18L), c("obs_mean","obs_q025","obs_q50","obs_q975","rep_mean","rep_q025","rep_q50","rep_q975","p_mid","n_greater","n_equal","n_less","n_draws","paired_diff_mean","paired_diff_q025","paired_diff_q50","paired_diff_q975","flag_tail")))
  d <- rep-obs
  qu <- function(z) unname(quantile(z,c(.025,.5,.975),names=FALSE))
  oq<-qu(obs); rq<-qu(rep); dq<-qu(d)
  c(obs_mean=mean(obs),obs_q025=oq[1],obs_q50=oq[2],obs_q975=oq[3],
    rep_mean=mean(rep),rep_q025=rq[1],rep_q50=rq[2],rep_q975=rq[3],
    p_mid=mean(rep>obs)+.5*mean(rep==obs),n_greater=sum(rep>obs),n_equal=sum(rep==obs),n_less=sum(rep<obs),n_draws=length(rep),
    paired_diff_mean=mean(d),paired_diff_q025=dq[1],paired_diff_q50=dq[2],paired_diff_q975=dq[3],
    flag_tail=as.numeric(mean(rep>obs)+.5*mean(rep==obs)<.025 | mean(rep>obs)+.5*mean(rep==obs)>.975))
}
node_rows <- list()
for (g in names(groups)) for (j in seq_len(J)) for (s in c("mean","sd","skewness","q05","q25","q50","q75","q95","below_L","above_U","outside","near_L","near_U","analytic_outside","epred_outside")) {
  ky <- paste(g,nodes[j],s,sep="|")
  obs <- if (s %in% c("analytic_outside","epred_outside")) 0 else a$node_obs[[ky]]
  node_rows[[length(node_rows)+1L]] <- data.frame(group=g,node=nodes[j],n=length(groups[[g]]),statistic=s,
    as.list(paired(obs,a$node_rep[[ky]])),check.names=FALSE)
}
node_table <- do.call(rbind,node_rows)
write_table(subset(node_table,statistic %in% c("mean","sd","skewness","q05","q25","q50","q75","q95","near_L","near_U")),"ppc_node_distribution.csv")
write_table(subset(node_table,statistic %in% c("below_L","above_U","outside","analytic_outside")),"ppc_boundary_rates.csv")
write_table(subset(node_table,statistic=="epred_outside"),"ppc_epred_range.csv")
analytic<-subset(node_table,statistic=="analytic_outside",select=c(group,node,n,rep_mean,rep_q025,rep_q50,rep_q975))
names(analytic)[4:7]<-paste0("analytic_",c("mean","q025","q50","q975"))
simulated<-subset(node_table,statistic=="outside",select=c(group,node,rep_mean,rep_q025,rep_q50,rep_q975))
names(simulated)[3:6]<-paste0("simulated_",c("mean","q025","q50","q975"))
analytic<-merge(analytic,simulated,by=c("group","node"),sort=FALSE)
analytic$simulated_minus_analytic<-analytic$simulated_mean-analytic$analytic_mean
write_table(analytic,"ppc_boundary_analytic_crosscheck.csv")
cor_rows <- list(); cross_rows <- list(); domain_rows <- list()
for (g in names(groups)) {
  for (k in seq_len(nrow(a$follow_pairs))) cor_rows[[length(cor_rows)+1L]] <- data.frame(group=g,source=a$follow_pairs[k,1],target=a$follow_pairs[k,2],
    as.list(paired(a$follow_obs[[g]][k],a$follow_rep[,g,k])),check.names=FALSE)
  for (k in seq_len(nrow(a$cross_pairs))) cross_rows[[length(cross_rows)+1L]] <- data.frame(group=g,source=a$cross_pairs[k,1],target=a$cross_pairs[k,2],
    self_crosswave=a$cross_pairs[k,1]==a$cross_pairs[k,2],as.list(paired(a$cross_obs[[g]][k],a$cross_rep[,g,k])),check.names=FALSE)
  for (s in dimnames(a$domain_rep)[[3]]) domain_rows[[length(domain_rows)+1L]] <- data.frame(group=g,statistic=s,
    as.list(paired(a$domain_obs[g,s],a$domain_rep[,g,s])),check.names=FALSE)
}
write_table(do.call(rbind,cor_rows),"ppc_followup_correlations.csv")
write_table(do.call(rbind,cross_rows),"ppc_crosswave_correlations.csv")
write_table(do.call(rbind,domain_rows),"ppc_domain_summaries.csv")
log("First-pass node and correlation summaries written")
hetero <- list()
for (g in setdiff(names(groups),"Overall")) for (j in seq_len(J)) {
  hetero[[length(hetero)+1L]] <- data.frame(statistic="residual_sd",group=g,node=nodes[j],pair=NA_character_,
    as.list(paired(a$resid_sd_obs[,g,j],a$resid_sd_rep[,g,j])),check.names=FALSE)
}
for (g in setdiff(names(groups),"Overall")) for (k in seq_len(nrow(a$follow_pairs))) {
  hetero[[length(hetero)+1L]] <- data.frame(statistic="residual_correlation",group=g,node=NA_character_,
    pair=paste(a$follow_pairs[k,],collapse="|"),as.list(paired(a$resid_cor_obs[,g,k],a$resid_cor_rep[,g,k])),check.names=FALSE)
}
for (j in seq_len(J)) hetero[[length(hetero)+1L]] <- data.frame(statistic="four_group_sd_max_min_ratio",group="four_groups",node=nodes[j],pair=NA_character_,
  as.list(paired(apply(a$resid_sd_obs[,,j],1,function(z) max(z)/min(z)),
                 apply(a$resid_sd_rep[,,j],1,function(z) max(z)/min(z)))),check.names=FALSE)
pairs <- combn(setdiff(names(groups),"Overall"),2,simplify=FALSE)
frobo <- frobr <- matrix(NA_real_,R,length(pairs))
for (k in seq_along(pairs)) {
  pp<-pairs[[k]]
  frobo[,k] <- sqrt(2*rowSums((a$resid_cor_obs[,pp[1],]-a$resid_cor_obs[,pp[2],])^2))
  frobr[,k] <- sqrt(2*rowSums((a$resid_cor_rep[,pp[1],]-a$resid_cor_rep[,pp[2],])^2))
  hetero[[length(hetero)+1L]] <- data.frame(statistic="pairwise_residual_correlation_frobenius",group="four_groups",node=NA_character_,pair=paste(pp,collapse="|"),
    as.list(paired(frobo[,k],frobr[,k])),check.names=FALSE)
}
hetero[[length(hetero)+1L]] <- data.frame(statistic="max_pairwise_residual_correlation_frobenius",group="four_groups",node=NA_character_,pair=NA_character_,
  as.list(paired(apply(frobo,1,max),apply(frobr,1,max))),check.names=FALSE)
write_table(do.call(rbind,hetero),"ppc_residual_heterogeneity.csv")
rm(hetero);gc(FALSE)
log("Residual covariance summaries written")
bin_records <- list(); bin_indices <- list(); bin_sets <- list()
add_bins <- function(g,j,kind,predictor,raw,continuous=FALSE) {
  ix <- groups[[g]]; z<-raw[ix]
  if (continuous) {
    cuts <- NULL
    for (nb in 5:1) {
      br <- unique(unname(quantile(z,probs=seq(0,1,length.out=nb+1),names=FALSE)))
      b <- if (length(br)==1L) rep(1L,length(z)) else as.integer(cut(z,breaks=br,include.lowest=TRUE,labels=FALSE))
      if (all(tabulate(b, max(1L, min(nb,length(br)-1L)))>=30L)) { cuts<-br; break }
    }
  } else {
    lev <- sort(unique(z)); counts <- tabulate(match(z,lev),nbins=length(lev)); m <- integer(length(lev)); num<-1L; acc<-0L
    for (h in seq_along(lev)) { m[h]<-num; acc<-acc+counts[h]; if (acc>=30L && h<length(lev)) {num<-num+1L;acc<-0L} }
    if (acc<30L && num>1L) m[m==num]<-num-1L
    b <- m[match(z,lev)]; cuts <- NULL
  }
  if (anyNA(b)) stop("Undefined bin")
  set_name <- paste(g,nodes[j],kind,predictor,sep="|")
  bin_sets[[set_name]] <<- character()
  for (h in sort(unique(b))) {
    ii<-ix[b==h]; id<-paste(set_name,h,sep="|")
    bin_indices[[id]] <<- ii
    bin_sets[[set_name]] <<- c(bin_sets[[set_name]],id)
    bin_records[[id]] <<- data.frame(bin_id=id,group=g,node=nodes[j],bin_type=kind,predictor=predictor,bin_number=h,
      n=length(ii),min_predictor=min(raw[ii]),max_predictor=max(raw[ii]),
      lower_cut=if (continuous && length(cuts)>1L) cuts[h] else min(raw[ii]),
      upper_cut=if (continuous && length(cuts)>1L) cuts[h+1L] else max(raw[ii]),
      split_ties=FALSE)
  }
}
for (g in names(groups)) for (j in seq_len(J)) {
  add_bins(g,j,"fitted_mu","mu_bar",a$mu_bar[,j],TRUE)
  add_bins(g,j,"own_baseline",nodes[j],X[,j],FALSE)
  if (j>3L) for (s in 1:3) add_bins(g,j,"baseline_support",nodes[s],X[,s],FALSE)
}
bin_def <- do.call(rbind,bin_records);write_table(bin_def,"ppc_binning_definition.csv")
B<-nrow(bin_def); obs_mean<-rep_mean<-obs_sd<-rep_sd<-matrix(NA_real_,R,B)
colnames(obs_mean)<-colnames(rep_mean)<-colnames(obs_sd)<-colnames(rep_sd)<-bin_def$bin_id
big_rep <- array(NA_real_,c(R,nrow(Y),J))
for (start in seq(1L,R,by=25L)) {
  end<-min(R,start+24L); path<-file.path(out,"cache",sprintf("ppc_batch_%04d_%04d.rds",start,end))
  cc<-readRDS(path)
  stopifnot(identical(cc$model_md5,a$model_md5),identical(cc$draw,a$selection$draw[start:end]),identical(cc$seed,cfg$ppc_seed))
  big_rep[start:end,,] <- cc$p
  for (bb in seq_len(end-start+1L)) {
    r<-start+bb-1L; eo<-Y-cc$m[bb,,]; er<-cc$p[bb,,]-cc$m[bb,,]
    for (j in seq_len(J)) {
      ids<-which(bin_def$node==nodes[j])
      for (h in ids) {
        ix<-bin_indices[[bin_def$bin_id[h]]]
        obs_mean[r,h]<-mean(eo[ix,j]);rep_mean[r,h]<-mean(er[ix,j])
        obs_sd[r,h]<-sd(eo[ix,j]);rep_sd[r,h]<-sd(er[ix,j])
      }
    }
  }
  log("Conditional residuals evaluated for draws ",start,"-",end)
  rm(cc);gc(FALSE)
}
cond_rows <- vector("list",B)
for (h in seq_len(B)) {
  mstat<-paired(obs_mean[,h],rep_mean[,h]); sstat<-paired(obs_sd[,h],rep_sd[,h])
  cond_rows[[h]]<-data.frame(bin_def[h,],mean_obs=unname(mstat["obs_mean"]),mean_rep=unname(mstat["rep_mean"]),
    mean_rep_q025=unname(mstat["rep_q025"]),mean_rep_q50=unname(mstat["rep_q50"]),mean_rep_q975=unname(mstat["rep_q975"]),
    mean_p_mid=unname(mstat["p_mid"]),mean_paired_diff=unname(mstat["paired_diff_mean"]),
    mean_paired_diff_q025=unname(mstat["paired_diff_q025"]),mean_paired_diff_q975=unname(mstat["paired_diff_q975"]),
    sd_obs=unname(sstat["obs_mean"]),sd_rep=unname(sstat["rep_mean"]),sd_rep_q025=unname(sstat["rep_q025"]),
    sd_rep_q50=unname(sstat["rep_q50"]),sd_rep_q975=unname(sstat["rep_q975"]),sd_p_mid=unname(sstat["p_mid"]),
    sd_paired_diff=unname(sstat["paired_diff_mean"]),sd_paired_diff_q025=unname(sstat["paired_diff_q025"]),
    sd_paired_diff_q975=unname(sstat["paired_diff_q975"]),check.names=FALSE)
}
cond_table<-do.call(rbind,cond_rows);write_table(cond_table,"ppc_conditional_residuals.csv")
set_rows <- list()
for (sn in names(bin_sets)) {
  h<-match(bin_sets[[sn]],bin_def$bin_id); bm<-obs_mean[,h,drop=FALSE];br<-rep_mean[,h,drop=FALSE]
  so<-obs_sd[,h,drop=FALSE];sr<-rep_sd[,h,drop=FALSE]
  ratio <- function(z) if (any(!is.finite(z)) || min(z)<=0) NA_real_ else max(z)/min(z)
  meta<-bin_def[h[1],c("group","node","bin_type","predictor")]
  set_rows[[length(set_rows)+1L]]<-data.frame(meta,statistic="max_abs_bin_mean_residual",
    as.list(paired(apply(abs(bm),1,max),apply(abs(br),1,max))),check.names=FALSE)
  set_rows[[length(set_rows)+1L]]<-data.frame(meta,statistic="max_min_bin_residual_sd_ratio",
    as.list(paired(apply(so,1,ratio),apply(sr,1,ratio))),check.names=FALSE)
}
write_table(do.call(rbind,set_rows),"ppc_conditional_set_summaries.csv")
rm(obs_mean,rep_mean,obs_sd,rep_sd);gc(FALSE)
log("All conditional bins summarized")
insample <- list(); ecdf_rows <- list()
for (j in seq_len(J)) {
  z <- big_rep[,,j]
  q <- apply(z,2,quantile,probs=c(.025,.10,.25,.75,.90,.975),names=FALSE)
  for (g in names(groups)) {
    ix<-groups[[g]]
    for (lev in c(50,80,95)) {
      ij<-switch(as.character(lev),"50"=c(3,4),"80"=c(2,5),"95"=c(1,6))
      insample[[length(insample)+1L]]<-data.frame(group=g,node=nodes[j],n=length(ix),interval=lev,
        rmse=sqrt(mean((Y[ix,j]-a$mu_bar[ix,j])^2)),mae=mean(abs(Y[ix,j]-a$mu_bar[ix,j])),
        bias=mean(Y[ix,j]-a$mu_bar[ix,j]),coverage=mean(Y[ix,j]>=q[ij[1],ix] & Y[ix,j]<=q[ij[2],ix]),
        mean_width=mean(q[ij[2],ix]-q[ij[1],ix]))
    }
    grid<-seq(min(c(z[,ix],Y[ix,j])),max(c(z[,ix],Y[ix,j])),length.out=81L)
    observed<-findInterval(grid,sort(Y[ix,j]))/length(ix)
    rr<-matrix(NA_real_,R,length(grid))
    for (r in seq_len(R)) rr[r,]<-findInterval(grid,sort(z[r,ix]))/length(ix)
    qs<-apply(rr,2,quantile,probs=c(.025,.5,.975),names=FALSE)
    ecdf_rows[[length(ecdf_rows)+1L]]<-data.frame(group=g,node=nodes[j],x=grid,observed=observed,
      rep_q025=qs[1,],rep_median=qs[2,],rep_q975=qs[3,])
  }
  log("In-sample coverage and ECDF node ",nodes[j])
}
write_table(do.call(rbind,insample),"in_sample_prediction_summary.csv")
write_table(do.call(rbind,ecdf_rows),"ppc_ecdf_pointwise.csv")
jsonlite::write_json(list(model_md5=a$model_md5,selection_md5=unname(tools::md5sum(file.path(out,"tables/ppc_draw_selection.csv"))),
  brms_version=as.character(packageVersion("brms")),posterior_version=as.character(packageVersion("posterior")),
  specification=paste0("joint_ppc_v2_n", length(nodes), "_draws", R, "_seed", cfg$ppc_seed),seed=cfg$ppc_seed,draws=R,batches=ceiling(R/25L)),
  file.path(out,"cache/ppc_cache_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
log("PPC summarization completed")

writeLines(capture.output(sessionInfo()), file.path(out,"sessionInfo.txt"))
