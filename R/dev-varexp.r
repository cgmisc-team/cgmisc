
#load("~/Dropbox/IgApek_2014/SP_gwas/perm_workspaceSPstatus_mm.Rd")
#data <- data4
#gkin <- data4.gkin
#top.names <- c("BICF2G630275788","BICF2P677609","BICF2G630276012","BICF2G630276743","BICF2S23711428","BICF2S23136666","BICF2G630552985","BICF2S23410286","BICF2P1401994","BICF2S23136652","BICF2P1436700","BICF2P215640","BICF2S2361733","TIGRP2P320304","TIGRP2P404896","BICF2G630772654")
#formula <- "IgA_5PERC_INT ~ age"
#res1 <- varexp(formula, data, gkin, top.names)

varexp <- function(formula, data, gkin, top.names, quiet=F, verbose=F) {
  require(GenABEL)
  require(hglm)
  # Function body
  formula <- as.formula(formula)
  trait <- all.vars(formula)[1]
  y <- data@phdata[,trait]
  # Phenotypic variance
  var.pheno <- var(y)
  
  # Variance due to fixed effects only
  model.fixed <- lm(formula, data=data)
  var.fixed <- var.pheno - var(model.fixed$residuals)
  aov.fixed <- anova(model.fixed)
  aov.sumSq <- aov.fixed$"Sum Sq"
  aov.fixed <- cbind(aov.fixed, PropVarExp=aov.sumSq/sum(aov.sumSq))
  se <- sqrt(deviance(model.fixed)/df.residual(model.fixed))
  # Variance explained by polygenic
  gkin[upper.tri(gkin)] <- t(gkin)[upper.tri(gkin)]
  gkin <- gkin * 2
  s <- svd(gkin)
  L <- s$u %*% diag(sqrt(s$d))
  X <- model.matrix(~ 1, data=data)
  model.random <- hglm(y=model.fixed$residuals, X=X, Z=L, conv=1e-8)
  res.random <- model.fixed$residuals - model.random$fv
  var.random <- var(model.fixed$residuals) - var(res.random)
    
  # Now the top snps
  idx <- which(data@gtdata@snpnames %in% top.names)
  top.gtps <- as.double(data[,idx]@gtdata)
  data.top <- data
  data.top@phdata <- cbind(data@phdata, top.gtps)
  tmp <- all.vars(formula)
  formula.top <- as.formula(paste(tmp[1],"~",paste(tmp[-1], collapse="+"), "+" ,paste(top.names, collapse="+"), sep=" "))
  
  model.fixed.top <- lm(formula.top, data=data.top)
  var.fixed.top <- var.pheno - var(model.fixed.top$residuals)
  prop.var.fixed.top <- var.fixed.top/var.pheno
  aov.fixed.top <- anova(model.fixed.top)
  aov.fixed.top.sumSq <- aov.fixed.top$"Sum Sq"
  aov.fixed.top <- cbind(aov.fixed.top, PropVarExp=aov.fixed.top.sumSq/sum(aov.fixed.top.sumSq))
  
  prop.var.fixed <- var.fixed / var.pheno
  prop.var.random <- var.random / var.pheno # ~= h2
  prop.var.indvar <- var(res.random) / var.pheno # Prop explained by individual variation (no diff. in fixed eff.) -- environment
  prop.var.top <- prop.var.fixed.top - prop.var.fixed

  if (!quiet) {
    cat("=============================================================================================\n")
    cat(sprintf("Proportion of variance explained by fixed effects: %1.4f\n", prop.var.fixed))
    cat(sprintf("Proportion of variance explained by polygenic effects (kinship h2): %1.4f\n", prop.var.random))
    cat(sprintf("Proportion of variance explained by individual variation: %1.4f\n", prop.var.indvar))
    if (verbose){ 
      cat("=============================================================================================\n")
      print(aov.fixed)
    }
    cat("=============================================================================================\n")
    cat(sprintf("Proportion of variance explained by top SNPs: %1.4f\n", prop.var.top))
    if (verbose) {
      print(aov.fixed.top)
    }
    cat("=============================================================================================\n")
  }
  result <- list(propVarFixed=prop.var.fixed, propVarRandom=prop.var.random, propVarInd=prop.var.indvar, propVarTop=prop.var.top)
  return(unlist(result))
}

