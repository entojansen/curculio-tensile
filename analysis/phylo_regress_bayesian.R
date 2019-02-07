# Andrew Jansen
# created iv.15.2018

# Miscellaneous Functions ------------------------------------------------------
#Import the data
GetData = function(file.name) {
  test.data = data.frame(read.csv(file.name, sep = ','))
  return(test.data)
}

#For ANOVA
SqDev = function(x) {var(x) * (length(x) - 1)}

# Analysis 1 -------------------------------------------------------------------
# Compute co-variance matrix from inverse
Analysis1 = function(phylo) {
  inv.phylo = MCMCglmm::inverseA(phylo, nodes="TIPS")
  A = solve(inv.phylo$Ainv)
  rownames(A) = rownames(inv.phylo$Ainv)
  print(A)
  return(A)
}

# Analysis 2 -------------------------------------------------------------------
# fit hierarchical phylogenetic model
Analysis2 = function(raw.data, A, priors) {
  phylo_model = brm(
    phen ~ spec_mean_cf + (1|phylo) + (1|species),
    data = raw.data,
    family = gaussian(),
    cov_ranef = list(phylo = A),
    prior = priors,
    sample_prior = TRUE#,
    #chains = 16,
    #cores = 8,
    #iter = 2000,
    #warmup = 1000
  )
  print(summary(phylo_model))
  plot(phylo_model)
  gg = plot(marginal_effects(phylo_model))[[1]]
  plot(gg + geom_point(inherit.aes = FALSE,
                  data = raw.data,
                  mapping = aes(x = spec_mean_cf, y = phen, color = phylo)))
  
  return(phylo_model)
}

# Analysis 3 -------------------------------------------------------------------
# Hypothesis testing and phylogenetic signal
Analysis3 = function(phylo_model) {
  hyp = paste("sd_phylo__Intercept^2 /",
              "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
  hyp = hypothesis(phylo_model, hyp, class = NULL)
  
  print(hyp)
  plot(hyp)
}

# Analysis 4 -------------------------------------------------------------------
# Update model and re-run phylogenetic signal
Analysis4 = function(raw.data, phylo_model) {
  phylo_updated = update(
    phylo_model,
    formula = ~ . + within_spec_cf,
    newdata = raw.data#,
    #chains = 16,
    #cores = 8,
    #iter = 2000,
    #warmup = 1000
  )
  
  print(summary(phylo_updated))
  #print(coef(phylo_updated))
  
  plot(phylo_updated)
  gg = plot(marginal_effects(phylo_updated))[[1]]
  plot(gg + geom_point(inherit.aes = FALSE,
                  data = raw.data,
                  mapping = aes(x = spec_mean_cf, y = phen, color = phylo)))
  gg = plot(marginal_effects(phylo_updated))[[2]]
  plot(gg + geom_point(inherit.aes = FALSE,
                  data = raw.data,
                  mapping = aes(x = within_spec_cf, y = phen, color = phylo)))
  
  plot(pp_check(phylo_model))
  plot(pp_check(phylo_updated))
  print(loo(phylo_model, phylo_updated, reloo = TRUE))
  
  hyp = paste("sd_phylo__Intercept^2 /",
              "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
  hyp = hypothesis(phylo_updated, hyp, class = NULL)
  
  print(hyp)
  plot(hyp)
  #pairs(phylo_updated)
}

# Main Function ----------------------------------------------------------------
main = function() {
  #Clear workspace
  rm(list = ls())
  setwd(paste("C:/Users/Andrew/Documents/Current Projects/",
              "Curculio Mechanical Testing/Tensile Testing/Analysis", sep=""))
  
  #Load gplots for plotCI
  library(ggplot2)
  library(ggeffects)
  library(ggthemes)
  library(RColorBrewer)
  library(car)
  library(ape)
  library(phylobase)
  library(MCMCglmm)
  library(brms)
  library(phytools)
  library(lmtest)
  library(nlme)
  
  #Create folder for plots
  dir.create("plots/rplots")
  
  #Minimal margins on graphs
  par(mar = c(4, 4, 1.1, 0.6) + 0.1)
  
  #Run analysis functions
  readline(prompt = "reading tensile.csv: press enter to continue.")
  
  #Manual comparisons
  #lnfmax vs lnendo
  raw.data = "tensile.csv"
  phi = GetData(raw.data)
  phi$ratio = phi$exo / phi$endo
  phi$lnlength = log(phi$length)
  phi$lnexo = log(1000000 * phi$exo)
  phi$lnendo = log(1000000 * phi$endo)
  phi$lntotal = log(1000000 * phi$total)
  phi$lnfmax = log(phi$fmax)
  phi$lnuts = log(phi$uts)
  phi$lnE = log(phi$E)
  phi$lnratio = log(phi$ratio)
  phi$phylo = phi$species
  phi$spec_mean_cf = with(phi, sapply(split(lnendo, phylo), mean)[phylo])
  phi$within_spec_cf = phi$lnendo - phi$spec_mean_cf
  phi$phen = phi$lnfmax
  row.names(phi) = phi$spec_ord
  
  print(phi)
  
  #lnfmax vs lnexo
  chi = GetData(raw.data)
  chi$ratio = chi$exo / chi$endo
  chi$lnlength = log(chi$length)
  chi$lnexo = log(1000000 * chi$exo)
  chi$lnendo = log(1000000 * chi$endo)
  chi$lntotal = log(1000000 * chi$total)
  chi$lnfmax = log(chi$fmax)
  chi$lnuts = log(chi$uts)
  chi$lnE = log(chi$E)
  chi$lnratio = log(chi$ratio)
  chi$phylo = chi$species
  chi$spec_mean_cf = with(chi, sapply(split(lnexo, phylo), mean)[phylo])
  chi$within_spec_cf = chi$lnexo - chi$spec_mean_cf
  chi$phen = chi$lnfmax
  row.names(chi) = chi$spec_ord
  
  #lnratio vs lnuts
  xi = GetData(raw.data)
  xi$ratio = xi$exo / xi$endo
  xi$lnlength = log(xi$length)
  xi$lnexo = log(1000000 * xi$exo)
  xi$lnendo = log(1000000 * xi$endo)
  xi$lntotal = log(1000000 * xi$total)
  xi$lnfmax = log(xi$fmax)
  xi$lnuts = log(xi$uts)
  xi$lnE = log(xi$E)
  xi$lnratio = log(xi$ratio)
  xi$phylo = xi$species
  xi$spec_mean_cf = with(xi, sapply(split(lnratio, phylo), mean)[phylo])
  xi$within_spec_cf = xi$lnratio - xi$spec_mean_cf
  xi$phen = xi$lnuts
  row.names(xi) = xi$spec_ord
  
  #lnlength vs lnE
  pi = GetData(raw.data)
  pi$ratio = pi$exo / pi$endo
  pi$lnlength = log(pi$length)
  pi$lnexo = log(1000000 * pi$exo)
  pi$lnendo = log(1000000 * pi$endo)
  pi$lntotal = log(1000000 * pi$total)
  pi$lnfmax = log(pi$fmax)
  pi$lnuts = log(pi$uts)
  pi$lnE = log(pi$E)
  pi$lnratio = log(pi$ratio)
  pi$phylo = pi$species
  pi$spec_mean_cf = with(pi, sapply(split(lnlength, phylo), mean)[phylo])
  pi$within_spec_cf = pi$lnlength - pi$spec_mean_cf
  pi$phen = pi$lnE
  row.names(pi) = pi$spec_ord
  
  tree = ape::read.tree(text = paste("((((aurivestis, uniformis), pardus),",
                                          "((caryae, nasicus), confusor)),",
                                          "((proboscideus, sulcatulus),",
                                          "(((humeralis, pardalis), victoriensis),",
                                          "longidens)));"))
  
  tree.trim = ape::read.tree(text = paste("((uniformis, caryae),",
                                          "((proboscideus, sulcatulus),",
                                          "(humeralis, victoriensis)));"))
  
  tree.star = ape::read.tree(text = paste("(((uniformis1, uniformis2, uniformis3,",
                                          "uniformis4, uniformis5), (caryae1,",
                                          "caryae2, caryae3, caryae4)),",
                                          "(((proboscideus1, proboscideus2,",
                                          "proboscideus3), (sulcatulus1,",
                                          "sulcatulus2, sulcatulus3, sulcatulus4)),",
                                          "((humeralis1, humeralis2), victoriensis1)));"))
  
  tree.trim$edge.length<-rep(1,nrow(tree.trim$edge))
  tree.star$edge.length<-rep(1,nrow(tree.star$edge))
  
  plot(tree.star, type = "phylogram", edge.width = 2)

  priors.phi = c(prior(normal(1,0.5), "b"),
                 prior(normal(6,2), "Intercept"),
                 prior(normal(0.2,0.05), "sd"),
                 prior(normal(0.2,0.05), "sigma"))

  priors.chi = c(prior(normal(0.5,0.25), "b"),
                 prior(normal(3,1), "Intercept"),
                 prior(normal(0.5,0.1), "sd"),
                 prior(normal(0.3,0.1), "sigma"))

  priors.xi = c(prior(normal(-0.5,0.25), "b"),
                prior(normal(4,1), "Intercept"),
                prior(normal(0.3,0.05), "sd"),
                prior(normal(0.3,0.05), "sigma"))

  priors.pi = c(prior(normal(1,0.5), "b"),
                prior(normal(6,1), "Intercept"),
                prior(normal(0.4,0.1), "sd"),
                prior(normal(0.4,0.1), "sigma"))
  
  #PGLS Est. via pgls.ives
  pgls_fit = pgls.Ives(tree.trim,
                       X = setNames(phi$lnexo, phi$species),
                       y = setNames(phi$lnfmax, phi$species))
  
  print(pgls_fit$a)
  
  pgls_null = pgls.Ives(tree.trim,
                       X = setNames(phi$lnexo, phi$species),
                       y = setNames(phi$lnfmax, phi$species),
                       fixed.b1 = 0)
  
  print(pgls_null$a)
  
  print(lrtest(pgls_null,pgls_fit))
  
  lm_mod = lme(lnfmax ~ lnendo,
               random = ~ 1|species,
               data = phi,
               correlation = corBrownian(1, tree.star),
               method = "ML")
  
  print(summary(lm_mod))
  
  lm_mod_alt = lme(lnfmax ~ lnendo,
               random = ~ 1|species,
               data = phi,
               #correlation = corBrownian(1, tree.star),
               method = "ML")
  
  print(summary(lm_mod_alt))
  
  print(lrtest(lm_mod,lm_mod_alt))
  
  #Comparison 1
  readline(prompt = "1st analysis: press enter to continue.")
  cov_mat.phi = Analysis1(tree)

  readline(prompt = "2nd analysis: press enter to continue.")
  mod.phi = Analysis2(phi, cov_mat.phi, priors.phi)

  readline(prompt = "3rd analysis: press enter to continue.")
  Analysis3(mod.phi)

  readline(prompt = "4th analysis: press enter to continue.")
  Analysis4(phi, mod.phi)

  #Comparison 2
  readline(prompt = "5th analysis: press enter to continue.")
  cov_mat.chi = Analysis1(tree)

  readline(prompt = "6th analysis: press enter to continue.")
  mod.chi = Analysis2(chi, cov_mat.chi, priors.chi)

  readline(prompt = "7th analysis: press enter to continue.")
  Analysis3(mod.chi)

  readline(prompt = "8th analysis: press enter to continue.")
  Analysis4(chi, mod.chi)
  
  #Comparison 3
  readline(prompt = "9th analysis: press enter to continue.")
  cov_mat.xi = Analysis1(tree)
  
  readline(prompt = "10th analysis: press enter to continue.")
  mod.xi = Analysis2(xi, cov_mat.xi, priors.xi)
  
  readline(prompt = "11th analysis: press enter to continue.")
  Analysis3(mod.xi)
  
  readline(prompt = "12th analysis: press enter to continue.")
  Analysis4(xi, mod.xi)
  
  #Comparison 4
  readline(prompt = "13th analysis: press enter to continue.")
  cov_mat.pi = Analysis1(tree)
  
  readline(prompt = "14th analysis: press enter to continue.")
  mod.pi = Analysis2(pi, cov_mat.pi, priors.pi)
  
  readline(prompt = "15th analysis: press enter to continue.")
  Analysis3(mod.pi)
  
  readline(prompt = "16th analysis: press enter to continue.")
  Analysis4(pi, mod.pi)
}

# Execute Code -----------------------------------------------------------------
main()