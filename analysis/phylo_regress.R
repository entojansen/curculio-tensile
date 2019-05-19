# Andrew Jansen
# created x.10.2018

# Miscellaneous Functions ------------------------------------------------------
# Import the data
GetData = function(file.name) {
  test.data = data.frame(read.csv(file.name, sep = ','))
  return(test.data)
}

# Analysis 1 -------------------------------------------------------------------
# Fits models with and without log transformation and phylogenetic covariance
Analysis1 = function(phylogeny, raw.data) {
  
  # Fit models with variable transformation and phylogenetic covariance
  simple.model = lme(phenotype ~ cofactor,
                     random = ~ 1 | species,
                     data = raw.data,
                     method = 'ML')
  
  simple.brownian = lme(phenotype ~ cofactor,
                        random = ~ 1 | species,
                        data = raw.data,
                        correlation = phylogeny,
                        method = 'ML')
  
  y.transformed = lme(log(phenotype) ~ cofactor,
                           random = ~ 1 | species,
                           data = raw.data,
                           method = 'ML')
  
  y.brownian = lme(log(phenotype) ~ cofactor,
                             random = ~ 1 | species,
                             data = raw.data,
                             correlation = phylogeny,
                             method = 'ML')
  
  x.transformed = lme(phenotype ~ log(cofactor),
                      random = ~ 1 | species,
                      data = raw.data,
                      method = 'ML')
  
  x.brownian = lme(phenotype ~ log(cofactor),
                   random = ~ 1 | species,
                   data = raw.data,
                   correlation = phylogeny,
                   method = 'ML')
  
  xy.transformed = lme(log(phenotype) ~ log(cofactor),
                       random = ~ 1 | species,
                       data = raw.data,
                       method = 'ML')
  
  xy.brownian = lme(log(phenotype) ~ log(cofactor),
                    random = ~ 1 | species,
                    data = raw.data,
                    correlation = phylogeny,
                    method = 'ML')
  
  # Summary of model fit and assumption verification, done for each model type
  # No transformation of 
  writeLines(strrep('-',80))
  writeLines('\n\nNo transformation of variables\n\n')
  writeLines(strrep('-',80))
  print(summary(simple.model))
  print(anova(simple.model))
  print(leveneTest(residuals(simple.model,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(simple.model, type = 'normalized')))
  r2sigma = r2beta(simple.model, method = 'sgv')
  print(r2sigma)
  
  print(summary(simple.brownian))
  print(anova(simple.brownian))
  print(leveneTest(residuals(simple.brownian,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(simple.brownian, type = 'normalized')))
  r2sigma.brownian = r2beta(simple.brownian, method = 'sgv')
  print(r2sigma.brownian)
  
  # R2 fit test of covariance structure for each transformation
  print(r2dt(x = r2sigma, y = r2sigma.brownian))
  
  # Natural log transformation of independent variable
  writeLines(strrep('-',80))
  writeLines('\n\nNatural log transformation of independent variable\n\n')
  writeLines(strrep('-',80))
  print(summary(y.transformed))
  print(anova(y.transformed))
  print(leveneTest(residuals(y.transformed,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(y.transformed, type = 'normalized')))
  r2sigma = r2beta(y.transformed, method = 'sgv')
  print(r2sigma)
  
  print(summary(y.brownian))
  print(anova(y.brownian))
  print(leveneTest(residuals(y.brownian,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(y.brownian, type = 'normalized')))
  r2sigma.brownian = r2beta(y.brownian, method = 'sgv')
  print(r2sigma.brownian)
  
  print(r2dt(x = r2sigma, y = r2sigma.brownian))
  
  # Natural log transformation of dependent variable
  writeLines(strrep('-',80))
  writeLines('\n\nNatural log transformation of dependent variable\n\n')
  writeLines(strrep('-',80))
  print(summary(x.transformed))
  print(anova(x.transformed))
  print(leveneTest(residuals(x.transformed,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(x.transformed, type = 'normalized')))
  r2sigma = r2beta(x.transformed, method = 'sgv')
  print(r2sigma)
  
  print(summary(x.brownian))
  print(anova(x.brownian))
  print(leveneTest(residuals(x.brownian,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(x.brownian, type = 'normalized')))
  r2sigma.brownian = r2beta(x.brownian, method = 'sgv')
  print(r2sigma.brownian)
  
  print(r2dt(x = r2sigma, y = r2sigma.brownian))
  
  # Natural log transformation of all variables
  writeLines(strrep('-',80))
  writeLines('\n\nNatural log transformation of both variables\n\n')
  writeLines(strrep('-',80))
  print(summary(xy.transformed))
  print(anova(xy.transformed))
  print(leveneTest(residuals(xy.transformed,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(xy.transformed, type = 'normalized')))
  r2sigma = r2beta(xy.transformed, method = 'sgv')
  print(r2sigma)
  
  print(summary(xy.brownian))
  print(anova(xy.brownian))
  print(leveneTest(residuals(xy.brownian,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(xy.brownian, type = 'normalized')))
  r2sigma.brownian = r2beta(xy.brownian, method = 'sgv')
  print(r2sigma.brownian)
  
  print(r2dt(x = r2sigma, y = r2sigma.brownian))
  
  print(lrtest(simple.brownian, x.brownian, y.brownian, xy.brownian))
  print(AIC(simple.brownian, x.brownian, y.brownian, xy.brownian))
}

# Analysis 2 -------------------------------------------------------------------
# Tests hypothesis that endo is correlated with fmax, and that exo is not
Analysis2 = function(phylogeny, raw.data) {
  
  # Fit models with and without phylogenetic covariance
  full.model = lme(log(fmax) ~ endo * exo,
                        random = ~ 1 | species,
                        data = raw.data,
                        correlation = phylogeny,
                        method = 'ML')
  
  full.simple = lme(log(fmax) ~ endo * exo,
                    random = ~ 1 | species,
                    data = raw.data,
                    method = 'ML')
  
  # Compare with and without phylogenetic covariance
  writeLines(strrep('-',80))
  writeLines('\n\nCompare with and without phylogenetic covariance\n\n')
  writeLines(strrep('-',80))
  print(summary(full.model))
  print(anova(full.model))
  print(leveneTest(residuals(full.model,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(full.model, type = 'normalized')))
  r2sigma = r2beta(full.model, method = 'sgv')
  print(r2sigma)
  
  print(summary(full.simple))
  print(anova(full.simple))
  print(leveneTest(residuals(full.simple,
                             type = 'normalized') ~ raw.data$species,
                   data = raw.data))
  print(shapiro.test(residuals(full.simple, type = 'normalized')))
  r2sigma.simple = r2beta(full.simple, method = 'sgv')
  print(r2sigma.simple)
  
  print(r2dt(x = r2sigma, y = r2sigma.simple))
  
  # Fit reduced models for hypothesis testing
  endo.model = lme(log(fmax) ~ endo,
                   random = ~ 1 | species,
                   data = raw.data,
                   correlation = phylogeny,
                   method = 'ML')
  
  exo.model = lme(log(fmax) ~ exo,
                  random = ~ 1 | species,
                  data = raw.data,
                  correlation = phylogeny,
                  method = 'ML')
  
  # R2 calculation for hypothesis testing
  writeLines(strrep('-',80))
  writeLines('\n\nR2 calculation for hypothesis testing\n\n')
  writeLines(strrep('-',80))
  r2m.full = r2beta(full.model, method = 'nsj')
  r2m.endo = r2beta(endo.model, method = 'nsj')
  r2m.exo = r2beta(exo.model, method = 'nsj')
  
  print(summary(r2m.full))
  print(summary(r2m.endo))
  print(summary(r2m.exo))
  
  print(r2dt(r2m.full, r2m.endo))
  print(r2dt(r2m.full, r2m.exo))
  print(r2dt(r2m.endo, r2m.exo, cor = FALSE))
  
  print(plot(r2m.full, r2mthd = 'nsj'))
  print(plot(r2sigma, r2mthd = 'sgv'))
  
  print(anova(full.model, endo.model))
  print(lrtest(full.model, endo.model))
  print(anova(full.model, exo.model))
  print(lrtest(full.model, exo.model))
  
  # Residual diagnositcs for full model
  print(plot(full.model,
             resid(., type = 'normalized') ~ fitted(.),
             abline = 0))
  
  dev.copy(pdf,
           'plots/rplots/exo_endo_fitted.pdf',
           width = 3,
           height = 3)
  dev.off()
  
  print(plot(full.model,
             resid(., type = 'normalized') ~ fitted(.) | species,
             abline = 0))
  
  dev.copy(pdf,
           'plots/rplots/exo_endo_resid_species.pdf',
           width = 4.5,
           height = 3)
  dev.off()
  
  model.residuals = residuals(full.model, type = 'normalized')
  qqnorm(model.residuals, ylab = 'Residual Quantiles', main = NULL)
  abline(0, 1, col = 'blue')
  
  dev.copy(pdf,
           'plots/rplots/exo_endo_resid_qqplot.pdf',
           width = 3,
           height = 3)
  dev.off()
}

# Analysis 3 -------------------------------------------------------------------
# Final model fitting, hypothesis testing, and phylogenetic signal
Analysis3 = function(phylogeny.vcv, phylogeny.tree, raw.data, axis.labels) {
  
  # Fit preferred model
  final.model = lme(phenotype ~ cofactor,
                   random = ~ 1 | species,
                   data = raw.data,
                   correlation = phylogeny.vcv,
                   method = 'ML')
  
  # REDUNDANT
  writeLines(strrep('-',80))
  writeLines('\n\nRedundant readout of chosen model, see above\n\n')
  writeLines(strrep('-',80))
  print(summary(final.model))
  
  # Phylogenetic signal testing in residuals
  model.residuals = residuals(final.model,
                              type = 'normalized',
                              level = 0:1)[, 2]
  trait = phylo4d(phylogeny.tree, tip.data = model.residuals)
  
  x = axis.labels[2]
  y = axis.labels[1]
  
  # Abouheif's Cmean
  writeLines(strrep('-',80))
  writeLines('\n\nAbouheif\'s Cmean\n\n')
  writeLines(strrep('-',80))
  abouheif.test = abouheif.moran(trait, method = 'oriAbouheif')
  print(abouheif.test)
  plot(abouheif.test)
  
  dev.copy(pdf,
           paste('plots/rplots/', y, x, 'abouheif.pdf', sep = '_'),
           width = 3,
           height = 3)
  dev.off()
  
  # Moran's I
  writeLines(strrep('-',80))
  writeLines('\n\nMoran\'s I\n\n')
  writeLines(strrep('-',80))
  moran.test = abouheif.moran(trait, method='Abouheif')
  print(moran.test)
  plot(abouheif.test)
  
  dev.copy(pdf,
           paste('plots/rplots/', y, x, 'moran.pdf', sep = '_'),
           width = 3,
           height = 3)
  dev.off()
  
  # If an error pops up for tree with no branch lengths, add a column to the
  # csv file, specimen, which has the species name and number.
  
  # Pagel's lambda
  writeLines(strrep('-',80))
  writeLines('\n\nPagel\'s lambda\n\n')
  writeLines(strrep('-',80))
  pagel.test = phylosig(tree = phylogeny.tree,
                        x = model.residuals,
                        method='lambda',
                        test=TRUE)
  print(pagel.test)
  
  # Blomberg's K
  writeLines(strrep('-',80))
  writeLines('\n\nBlomberg\'s kappa\n\n')
  writeLines(strrep('-',80))
  blomberg.test = phylosig(tree = phylogeny.tree,
                           x = model.residuals,
                           method='K',
                           test=TRUE)
  print(blomberg.test)
  
  # R2 hypothesis testing
  writeLines(strrep('-',80))
  writeLines('\n\nR2 hypothesis testing\n\n')
  writeLines(strrep('-',80))
  r2.marginal = r2beta(final.model, method = 'nsj', partial = FALSE)
  r2.sigma = r2beta(final.model, method = 'sgv', partial = FALSE)
  
  print(summary(r2.marginal))
  print(summary(r2.sigma))
  
}

# Analysis 4 -------------------------------------------------------------------
# Data visualization and graphing
Analysis4 = function(phylogeny, raw.data, axis.labels) {
  
  # Fit preferred model
  final.model = lme(phenotype ~ cofactor,
                    random = ~ 1 | species,
                    data = raw.data,
                    correlation = phylogeny,
                    method = 'ML')
  
  # create data.frame with new values for predictors
  cofactor = raw.data$cofactor
  new.dat <- data.frame(cofactor)
  
  # predict normalized
  new.dat$pred <- predict(final.model, newdata = new.dat, level = 0)
  
  # create design matrix
  Designmat <- model.matrix(eval(eval(final.model$call$fixed)[-2]),
                            new.dat[-ncol(new.dat)])
  
  # compute standard error for predictions
  predvar = diag(Designmat %*% final.model$varFix %*% t(Designmat))
  new.dat$SE = sqrt(predvar) 
  new.dat$SE2 = sqrt(predvar + final.model$sigma ^ 2)
  
  # plot marginal effects
  gg = ggplot(new.dat, aes(x = cofactor, y = pred)) + 
    geom_ribbon(aes(ymin = pred - 2 * SE2, ymax = pred + 2 * SE2),
                alpha = 0.2,
                fill = 'black') +
    geom_ribbon(aes(ymin = pred - 2 * SE, ymax = pred + 2 * SE),
                alpha = 0.2,
                fill = 'black') +
    geom_line(color = 'black', size = 1) +
    geom_point(data = raw.data,
               aes(x = cofactor, y = phenotype, color = species)) +
    scale_y_continuous(axis.labels[1]) +
    scale_x_continuous(axis.labels[2])
  plot(gg)
  
  x = axis.labels[2]
  y = axis.labels[1]
  
  dev.copy(pdf,
           paste('plots/rplots/', y, x, 'marginal_effects.pdf', sep = '_'),
           width = 4.5,
           height = 3)
  dev.off()
  
  # plot transformed conditional residuals vs. fitted values
  print(plot(final.model,
             resid(., type = 'normalized') ~ fitted(.),
             abline = 0))
  
  dev.copy(pdf,
           paste('plots/rplots/', y, x, 'fitted.pdf', sep = '_'),
           width = 3,
           height = 3)
  dev.off()
  
  print(plot(final.model,
             resid(., type = 'normalized') ~ fitted(.) | species,
             abline = 0))
  
  dev.copy(pdf,
           paste('plots/rplots/', y, x, 'species.pdf', sep = '_'),
           width = 4.5,
           height = 3)
  dev.off()
  
  print(plot(final.model, species ~ resid(., type = 'normalized')))
  
  dev.copy(pdf,
           paste('plots/rplots/', y, x, 'boxplot.pdf', sep = '_'),
           width = 4.5,
           height = 3)
  dev.off()
  
  model.residuals = residuals(final.model, type = 'normalized')
  qqnorm(model.residuals, ylab = 'Residual Quantiles', main = NULL)
  abline(0, 1, col = 'blue')
  
  dev.copy(pdf,
           paste('plots/rplots/', y, x, 'qqplot.pdf', sep = '_'),
           width = 3,
           height = 3)
  dev.off()
}

# Main Function ----------------------------------------------------------------
main = function() {
  # Clear workspace
  rm(list = ls())
  setwd(paste('C:/Users/Andrew/Documents/GitHub/',
              'curculio-tensile/analysis', sep=''))
  
  # Load gplots for plotCI
  library(ggplot2)
  library(ggeffects)
  library(ggthemes)
  library(RColorBrewer)
  library(car)
  library(ape)
  library(phylobase)
  library(adephylo)
  library(MCMCglmm)
  library(phytools)
  library(lmtest)
  library(nlme)
  library(r2glmm)
  library(nortest)
  
  
  # Create folder for plots
  dir.create('plots/rplots')
  
  # Minimal margins on graphs
  par(mar = c(4, 4, 1.1, 0.6) + 0.1)
  
  # Import raw data and phylogeny
  file.name = 'tensile.csv'
  raw.data = GetData(file.name)
  raw.data$ratio = raw.data$exo / raw.data$endo
  row.names(raw.data) = raw.data$specimen
  
  tree.text = paste('(((uniformis1, uniformis2, uniformis3,',
                    'uniformis4, uniformis5), (caryae1,',
                    'caryae2, caryae3, caryae4)),',
                    '(((proboscideus1, proboscideus2,',
                    'proboscideus3), (sulcatulus1,',
                    'sulcatulus2, sulcatulus3, sulcatulus4)),',
                    '((humeralis1, humeralis2), victoriensis1)));')
  tree.star = ape::read.tree(text = tree.text)
  tree.star$edge.length = rep(1,nrow(tree.star$edge))

  # Generate phylogenetic variance-covariance matrix, Brownian motion
  phylo.vcv = corBrownian(1, tree.star)
  
  plot(tree.star, type = 'phylogram', edge.width = 2)
  
  dev.copy(pdf,
           'plots/rplots/phylogeny.pdf',
           width = 6,
           height = 3)
  dev.off()
  
  # All output to file
  sink(file = 'phylo_regress_output.txt', append =FALSE)
  print('PHYLOGENETIC REGRESSION OUTPUT')
  sink()
  
  sink(file = 'phylo_regress_output.txt', append = TRUE)
  
  # Model exploration and initial testing of phylogenetic covariance
  # Comparison 1a, model covariance and transformation
  readline(prompt = 'fmax vs. endo: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 1a: fmax vs. endo')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$fmax
  raw.data$cofactor = raw.data$endo
  Analysis1(phylo.vcv, raw.data)

  # Comparison 1b, model covariance and transformation
  readline(prompt = 'fmax vs. exo: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 1b: fmax vs. exo')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$fmax
  raw.data$cofactor = raw.data$exo
  Analysis1(phylo.vcv, raw.data)

  # Comparison 1c, model covariance and transformation
  readline(prompt = 'uts vs. ratio: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 1c: uts vs. ratio')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$uts
  raw.data$cofactor = raw.data$ratio
  Analysis1(phylo.vcv, raw.data)
  
  # Comparison 2a, model covariance and transformation
  readline(prompt = 'ufs vs. Etanl: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 2a: ufs vs. Etanl')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$ufs
  raw.data$cofactor = raw.data$Etanl
  Analysis1(phylo.vcv, raw.data)

  # Comparison 2b, model covariance and transformation
  readline(prompt = 'ufs vs. Esec: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 2b: ufs vs. Esec')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$ufs
  raw.data$cofactor = raw.data$uts / raw.data$ufs
  Analysis1(phylo.vcv, raw.data)

  # Comparison 2c, model covariance and transformation
  readline(prompt = 'U vs. Etanl: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 2c: U vs. Etanl')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$U
  raw.data$cofactor = raw.data$Etanl
  Analysis1(phylo.vcv, raw.data)

  # Comparison 2d, model covariance and transformation
  readline(prompt = 'U vs. Esec: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 2d: U vs. Esec')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$U
  raw.data$cofactor = raw.data$uts / raw.data$ufs
  Analysis1(phylo.vcv, raw.data)

  # Comparison 3a, model covariance and transformation
  readline(prompt = 'Etan low vs. length: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 3a: Etan low vs. length')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$Etanl
  raw.data$cofactor = raw.data$length
  Analysis1(phylo.vcv, raw.data)

  # Comparison 3b, model covariance and transformation
  readline(prompt = 'Esec vs. length: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 3b: Esec vs. length')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = raw.data$uts / raw.data$ufs
  raw.data$cofactor = raw.data$length
  Analysis1(phylo.vcv, raw.data)

  # Comparison 1 Test, compound cross sectional area
  readline(prompt = 'log(fmax) ~ endo + exo: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Comparison 5: log(fmax) ~ endo + exo')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  Analysis2(phylo.vcv, raw.data)

  # Final model fitting, hypothesis testing, and phylogenetic signal
  # Model 1a: log(fmax) ~ endo
  readline(prompt = 'fmax vs. endo: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 1a: log(fmax) ~ log(endo)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$fmax)
  raw.data$cofactor = log(raw.data$endo)
  comparison = c('ln(Fmax [N])', 'ln(Endocuticle CSA [mm2])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 1b: log(fmax) ~ exo
  readline(prompt = 'fmax vs. exo: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 1b: log(fmax) ~ exo')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$fmax)
  raw.data$cofactor = log(raw.data$exo)
  comparison = c('ln(Fmax [N])', 'ln(Exocuticle CSA [mm2])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 1c:log(uts) ~ log(ratio)
  readline(prompt = 'uts vs. ratio: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 1c:log(uts) ~ log(ratio)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$uts)
  raw.data$cofactor = log(raw.data$ratio)
  comparison = c('ln(Ultimate Strength [MPa])',
                 'ln(Exocuticle-Endocuticle Ratio)')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 2a: log(ufs) ~ log(Etanl)
  readline(prompt = 'ufs vs. Etanl: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 2a: log(ufs) ~ log(Etanl)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$ufs)
  raw.data$cofactor = log(raw.data$Etanl)
  comparison = c('ln(Ultimate Strain [mm2])', 'ln(Elastic Modulus [MPa])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 2b: log(ufs) ~ log(Esec)
  readline(prompt = 'ufs vs. Esec: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 2b: log(ufs) ~ log(Esec)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$ufs)
  raw.data$cofactor = log(raw.data$Esec)
  comparison = c('ln(Ultimate Strain [mm2])', 'ln(Secant Modulus [MPa])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 2c: log(U) ~ log(Etanl)
  readline(prompt = 'U vs. Etanl: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 2a: log(U) ~ log(Etanl)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$U)
  raw.data$cofactor = log(raw.data$Etanl)
  comparison = c('ln(Work of Fracture [Jm-3])', 'ln(Elastic Modulus [MPa])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 2d: log(U) ~ log(Esec)
  readline(prompt = 'U vs. Esec: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 2b: log(U) ~ log(Esec)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$U)
  raw.data$cofactor = log(raw.data$Esec)
  comparison = c('ln(Work of Fracture [Jm-3])', 'ln(Secant Modulus [MPa])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 3a: log(Etanl) ~ log(length)
  readline(prompt = 'Etan low vs. length: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 3a: log(Etanl) ~ log(length)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$Etanl)
  raw.data$cofactor = log(raw.data$length)
  comparison = c('ln(Elastic Modulus [MPa])', 'ln(Rostrum Length [mm2])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  # Model 3b: log(Esec) ~ log(length)
  readline(prompt = 'Esec vs. length: press enter to continue.')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  writeLines('Model 3b: log(Esec) ~ log(length)')
  writeLines(strrep('=',80))
  writeLines(strrep('=',80))
  raw.data$phenotype = log(raw.data$Esec)
  raw.data$cofactor = log(raw.data$length)
  comparison = c('ln(Secant Modulus [MPa])', 'ln(Rostrum Length [mm2])')
  Analysis3(phylo.vcv, tree.star, raw.data, comparison)
  Analysis4(phylo.vcv, raw.data, comparison)

  #Close output file
  sink()
}

# Execute Code -----------------------------------------------------------------
main()