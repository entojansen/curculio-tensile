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
# ADD NAMES FOR LABELS ON GRAPHS TO FUNCTION INPUT
Analysis1 = function(phylo, raw.data, axis.labels) {
  
  lm_mod = lme(phen ~ cf,
               random = ~ 1|species,
               data = raw.data,
               correlation = corBlomberg(6, phylo, fixed = TRUE),
               method = 'ML')
  
  lm_mod_alt = lme(phen ~ cf,
                   random = ~ 1|species,
                   data = raw.data,
                   method = 'ML')
  
  print(summary(lm_mod))
  print(summary(lm_mod_alt))
  print(anova(lm_mod))
  print(anova(lm_mod_alt))
  
  print(lrtest(lm_mod, lm_mod_alt))
  
  #create data.frame with new values for predictors
  #more than one predictor is possible
  cf = raw.data$cf
  new.dat <- data.frame(cf)
  
  #predict response
  new.dat$pred <- predict(lm_mod, newdata=new.dat, level=0)
  
  #create design matrix
  Designmat <- model.matrix(eval(eval(lm_mod$call$fixed)[-2]),
                            new.dat[-ncol(new.dat)])
  
  #compute standard error for predictions
  predvar = diag(Designmat %*% lm_mod$varFix %*% t(Designmat))
  new.dat$SE = sqrt(predvar) 
  new.dat$SE2 = sqrt(predvar+lm_mod$sigma^2)
  
  gg = ggplot(new.dat,aes(x=cf,y=pred)) + 
    geom_ribbon(aes(ymin=pred-2*SE2,ymax=pred+2*SE2),alpha=0.2,fill='black') +
    geom_ribbon(aes(ymin=pred-2*SE,ymax=pred+2*SE),alpha=0.2,fill='black') +
    geom_line(color='black', size=1) +
    geom_point(data=raw.data,aes(x=cf,y=phen,color=species)) +
    scale_y_continuous(axis.labels[2]) +
    scale_x_continuous(axis.labels[1])
  plot(gg)
  
  print(plot(lm_mod, resid(., type = 'n') ~ fitted(.), abline = 0))
  print(plot(lm_mod_alt, resid(., type = 'n') ~ fitted(.), abline = 0))
  
  results = leveneTest(residuals(lm_mod,
                                 type = 'normalized') ~ raw.data$species,
                       data = raw.data)
  print(results)
  
  #residuals grouped/labeled by species or by specimen ID, respectively
  res.mod = residuals(lm_mod, type = 'normalized')
  phy.mod = residuals(lm_mod, type='normalized', level = 0:1)[,2]
  
  qqnorm(res.mod, ylab = 'Residual Quantiles', main = 'normalized')
  qqline(res.mod, col = 'blue')
  print(shapiro.test(res.mod))
  print(ad.test(res.mod))
  
  res.alt = residuals(lm_mod_alt, type = 'normalized')
  phy.alt = residuals(lm_mod, type='normalized', level = 0:1)[,2]
  
  qqnorm(res.alt, ylab = 'Residual Quantiles', main = 'normalized')
  qqline(res.alt, col = 'blue')
  print(shapiro.test(res.alt))
  print(ad.test(res.alt))
  
  print(phylosig(phylo, phy.mod, test = TRUE))
  print(phylosig(phylo, phy.alt, test = TRUE))
  
  r2m = r2beta(lm_mod, method = 'nsj')
  r2sigma = r2beta(lm_mod, method = 'sgv')
  r2m_alt = r2beta(lm_mod_alt, method = 'nsj')
  r2sigma_alt = r2beta(lm_mod_alt, method = 'sgv')
  print(r2m)
  print(r2sigma)
  print(r2m_alt)
  print(r2sigma_alt)
  print(r2dt(x = r2m, y = r2m_alt))
  print(r2dt(x = r2sigma, y = r2sigma_alt))
}

# Analysis 2 -------------------------------------------------------------------
#
Analysis2 = function(phylo, raw.data) {
  
  lm_mod = lme(lnfmax ~ lnendo + lnexo,
               random = ~ 1|species,
               data = raw.data,
               correlation = corBrownian(1, phylo),
               method = 'ML')
  
  lm_mod_alt = lme(lnfmax ~ lnendo,
                   random = ~ 1|species,
                   data = raw.data,
                   correlation = corBrownian(1, phylo),
                   method = 'ML')
  
  print(summary(lm_mod))
  print(summary(lm_mod_alt))
  print(anova(lm_mod))
  print(anova(lm_mod_alt))
  
  print(lrtest(lm_mod, lm_mod_alt))
  
  print(plot(lm_mod))
  print(plot(lm_mod_alt))
  
  
  results = leveneTest(residuals(lm_mod,
                                 type = 'normalized') ~ raw.data$species,
                       data = raw.data)
  print(results)
  
  res.mod = residuals(lm_mod, type = 'normalized')
  phy.mod = residuals(lm_mod, type='normalized', level = 0:1)[,2]
  
  qqnorm(res.mod, ylab = 'Residual Quantiles', main = 'normalized')
  qqline(res.mod, col = 'blue')
  print(shapiro.test(res.mod))
  
  res.alt = residuals(lm_mod_alt, type = 'normalized')
  phy.alt = residuals(lm_mod, type='normalized', level = 0:1)[,2]

  qqnorm(res.alt, ylab = 'Residual Quantiles', main = 'normalized')
  qqline(res.alt, col = 'blue')
  print(shapiro.test(res.alt))

  print(phylosig(phylo, phy.mod, test = TRUE))
  print(phylosig(phylo, phy.alt, test = TRUE))

  r2m = r2beta(lm_mod, method = 'nsj')
  r2sigma = r2beta(lm_mod, method = 'sgv')
  r2m_alt = r2beta(lm_mod_alt, method = 'nsj')
  r2sigma_alt = r2beta(lm_mod_alt, method = 'sgv')
  print(r2m)
  print(r2sigma)
  print(r2m_alt)
  print(r2sigma_alt)
  print(r2dt(x = r2m, y = r2m_alt))
  print(r2dt(x = r2sigma, y = r2sigma_alt))
  print(plot(r2m, r2mthd = 'nsj'))
  print(plot(r2sigma, r2mthd = 'sgv'))
}

# Analysis 3 -------------------------------------------------------------------
#
Analysis3 = function() {

}

# Analysis 4 -------------------------------------------------------------------
#
Analysis4 = function() {

}

# Main Function ----------------------------------------------------------------
main = function() {
  #Clear workspace
  rm(list = ls())
  setwd(paste('C:/Users/Andrew/Documents/Current Projects/',
              'Curculio Mechanical Testing/Tensile Testing/Analysis', sep=''))
  
  #Load gplots for plotCI
  library(ggplot2)
  library(ggeffects)
  library(ggthemes)
  library(RColorBrewer)
  library(car)
  library(ape)
  library(phylobase)
  library(MCMCglmm)
  library(phytools)
  library(lmtest)
  library(nlme)
  library(r2glmm)
  library(nortest)
  
  #Create folder for plots
  dir.create('plots/rplots')
  
  #Minimal margins on graphs
  par(mar = c(4, 4, 1.1, 0.6) + 0.1)
  
  #Manual comparisons
  #lnfmax vs lnendo
  raw.data = 'tensile.csv'
  phi = GetData(raw.data)
  phi$ratio = phi$exo / phi$endo
  phi$lnlength = log(phi$length)
  phi$lnexo = (1000000 * phi$exo)
  phi$lnendo = (1000000 * phi$endo)
  phi$lntotal = (1000000 * phi$total)
  phi$lnfmax = log(phi$fmax)
  phi$lnuts = log(phi$uts)
  phi$lnE = log(phi$E)
  phi$lnratio = log(phi$ratio)
  phi$phylo = phi$species
  phi$phen = phi$lnfmax
  phi$cf = phi$lnendo
  row.names(phi) = phi$specimen
  
  #lnfmax vs lnexo
  chi = GetData(raw.data)
  chi$ratio = chi$exo / chi$endo
  chi$lnlength = log(chi$length)
  chi$lnexo = sqrt(1000000 * chi$exo)
  chi$lnendo = sqrt(1000000 * chi$endo)
  chi$lntotal = sqrt(1000000 * chi$total)
  chi$lnfmax = log(chi$fmax)
  chi$lnuts = (chi$uts)
  chi$lnE = log(chi$E)
  chi$lnratio = log(chi$ratio)
  chi$phylo = chi$species
  chi$spec_mean_cf = with(chi, sapply(split(lnexo, phylo), mean)[phylo])
  chi$within_spec_cf = chi$lnexo - chi$spec_mean_cf
  chi$phen = chi$lnfmax
  chi$cf = chi$lntotal
  row.names(chi) = chi$specimen
  
  #lnratio vs lnuts
  xi = GetData(raw.data)
  xi$ratio = xi$exo / xi$endo
  xi$lnlength = log(xi$length)
  xi$lnexo = sqrt(1000000 * xi$exo)
  xi$lnendo = sqrt(1000000 * xi$endo)
  xi$lntotal = sqrt(1000000 * xi$total)
  xi$lnfmax = log(xi$fmax)
  xi$lnuts = log(xi$uts)
  xi$lnE = log(xi$E)
  xi$lnratio = log(xi$ratio)
  xi$phylo = xi$species
  xi$spec_mean_cf = with(xi, sapply(split(lnratio, phylo), mean)[phylo])
  xi$within_spec_cf = xi$lnratio - xi$spec_mean_cf
  xi$phen = xi$lnuts
  xi$cf = xi$lnratio
  row.names(xi) = xi$specimen
  
  #lnlength vs lnE
  pi = GetData(raw.data)
  pi$ratio = pi$exo / pi$endo
  pi$lnlength = log(pi$length)
  pi$lnexo = sqrt(1000000 * pi$exo)
  pi$lnendo = sqrt(1000000 * pi$endo)
  pi$lntotal = sqrt(1000000 * pi$total)
  pi$lnfmax = log(pi$fmax)
  pi$lnuts = log(pi$uts)
  pi$lnE = log(pi$E)
  pi$lnratio = log(pi$ratio)
  pi$phylo = pi$species
  pi$spec_mean_cf = with(pi, sapply(split(lnlength, phylo), mean)[phylo])
  pi$within_spec_cf = pi$lnlength - pi$spec_mean_cf
  pi$phen = pi$lnE
  pi$cf = pi$lnlength
  row.names(pi) = pi$specimen
  
  tree.text = paste('((((aurivestis, uniformis), pardus),',
                   '((caryae, nasicus), confusor)),',
                   '((proboscideus, sulcatulus),',
                   '(((humeralis, pardalis), victoriensis),',
                   'longidens)));')
  tree = ape::read.tree(text = tree.text)
  
  trim.text = paste('((uniformis, caryae),',
                    '((proboscideus, sulcatulus),',
                    '(humeralis, victoriensis)));')
  tree.trim = ape::read.tree(text = trim.text)
  
  star.text = paste('(((uniformis1, uniformis2, uniformis3,',
                    'uniformis4, uniformis5), (caryae1,',
                    'caryae2, caryae3, caryae4)),',
                    '(((proboscideus1, proboscideus2,',
                    'proboscideus3), (sulcatulus1,',
                    'sulcatulus2, sulcatulus3, sulcatulus4)),',
                    '((humeralis1, humeralis2), victoriensis1)));')
  tree.star = ape::read.tree(text = star.text)
  
  tree.trim$edge.length<-rep(1,nrow(tree.trim$edge))
  tree.star$edge.length<-rep(1,nrow(tree.star$edge))
  
  plot(tree.star, type = 'phylogram', edge.width = 2)
  
  #Comparison 1
  # readline(prompt = '1st analysis: press enter to continue.')
  # Analysis1(tree.star, phi, c('r-endo','lnfmax'))
  # Analysis1(tree.star, chi, c('r-exo','lnfmax'))
  Analysis1(tree.star, xi, c('lnratio','lnuts'))
  # Analysis1(tree.star, pi, c('lnlength','lnE'))

  # Comparison 2
  # readline(prompt = '2nd analysis: press enter to continue.')
  # Analysis2(tree.star, phi)
}

# Execute Code -----------------------------------------------------------------
main()