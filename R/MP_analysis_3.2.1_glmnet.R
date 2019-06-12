library(glmnet)
library(FLCore)
library(doParallel)

### ------------------------------------------------------------------------ ###
### glmnet with extended parameters ####
### ------------------------------------------------------------------------ ###

### load data
glm_data <- readRDS("output/stats_3.2.1.rds")
# glm_data <- glm_data[order(glm_data$scenario), ]
# ### load new stats
# glm_data_stats <- readRDS("output/stats_corrected_scn.RDS")
# ### subset
# glm_data_stats <- glm_data_stats[glm_data_stats$scenario %in% glm_data$scenario, ]

library(doParallel)
library(ggplot2)
library(tidyr)
library(dplyr)
library(glmnet)
cl <- makeCluster(4)
clusterEvalQ(cl = cl, expr = {library(glmnet)})
registerDoParallel(cl)

### select predictor and response variables
var_pred <- c("a", "b", "linf", "k", "t0", "a50", "alpha", "beta", "spr0",
              "lopt", "r", "rc", "m", "mk", "fmsym", "bmsyb0")
var_resp <- c("f_rel", "ssb_rel", "collapse_total", "ssb_below_blim_total",
              "yield_rel_MSY", "iav")
#var_resp2 <- c("collapse_total", "yield_rel_MSY")

### first: test lasso only
ow <- cv.glmnet(
  x = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred]),
  y = as.matrix(glm_data[glm_data$fhist == "one-way", var_resp]),
  foldid = 1:29, nlambda = 1000, family = "mgaussian", alpha = 1)
coef(ow, s = "lambda.min")
coef(ow, s = "lambda.1se")
### mean squared error
min(ow$cvm)

### fit lm separately to each of the response variables
lms_ow <- lapply(var_resp, function(x) {
  lm(glm_data[glm_data$fhist == "one-way", x] ~
              glm_data[glm_data$fhist == "one-way", "k"])
})
names(lms_ow) <- var_resp
lapply(lms_ow, summary)
lms <- lapply(var_resp, function(x) {
  tmp <- lm(glm_data[glm_data$fhist == "one-way", x] ~
              glm_data[glm_data$fhist == "one-way", "k"])
  tmp <- as.data.frame(t(data.frame(coef(tmp))))
  colnames(tmp) <- c("Intercept", "k")
  tmp$var <- x
  tmp
})
names(lms) <- var_resp
lms <- do.call(rbind, lms)
res_glmnet <- coef(ow, s = "lambda.min")

### save data
saveRDS(object = list(data = glm_data, glmnet = res_glmnet, lm = lms_ow, 
                      var_resp = var_resp, var_pred = var_pred), 
        file = "output/glmnet_results_ow.rds")

### plot
plot_df <- glm_data[glm_data$fhist == "one-way", c(var_resp, "k")] %>%
  gather(key = "response", value = "value", 1:5)
ggplot(data = plot_df, aes(x = k, y = value)) +
  geom_point() + theme_bw() + facet_wrap(~ response, scales = "free_y")
  
  
p1 <- ggplot(data = glm_data[glm_data$fhist == "one-way", ],
               aes(x = k, y = f_rel)) +
  geom_point() + theme_bw() + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(intercept = res_glmnet$f_rel["(Intercept)",],
              slope = res_glmnet$f_rel["k",], col = "red") +
  geom_abline(data = lms[lms$var == "f_rel", ], 
              aes(intercept = Intercept, slope = k)) + 
  labs(y = "F/Fmsy")
p2 <- ggplot(data = glm_data[glm_data$fhist == "one-way", ],
             aes(x = k, y = ssb_rel)) +
  geom_point() + theme_bw() + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(intercept = res_glmnet$ssb_rel["(Intercept)",],
              slope = res_glmnet$ssb_rel["k",], col = "red") +
  geom_abline(data = lms[lms$var == "ssb_rel", ], 
              aes(intercept = Intercept, slope = k)) + 
  labs(y = "B/Bmsy")
p3 <- ggplot(data = glm_data[glm_data$fhist == "one-way", ],
             aes(x = k, y = collapse_total)) +
  geom_point() + theme_bw() + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(intercept = -res_glmnet$collapse_total["(Intercept)",],
              slope = -res_glmnet$collapse_total["k",], col = "red") +
  geom_abline(data = lms[lms$var == "collapse_total", ], 
              aes(intercept = -Intercept, slope = -k)) + 
  labs(y = "collapse risk") +
  scale_y_reverse()
p4 <- ggplot(data = glm_data[glm_data$fhist == "one-way", ],
             aes(x = k, y = ssb_below_blim_total)) +
  geom_point() + theme_bw() + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(intercept = -res_glmnet$ssb_below_blim_total["(Intercept)",],
              slope = -res_glmnet$ssb_below_blim_total["k",], col = "red") +
  geom_abline(data = lms[lms$var == "ssb_below_blim_total", ], 
              aes(intercept = -Intercept, slope = -k)) + 
  labs(y = "risk Blim") +
  scale_y_reverse()
p5 <- ggplot(data = glm_data[glm_data$fhist == "one-way", ],
             aes(x = k, y = yield_rel_MSY)) +
  geom_point() + theme_bw() + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(intercept = res_glmnet$yield_rel_MSY["(Intercept)",],
              slope = res_glmnet$yield_rel_MSY["k",], col = "red") +
  geom_abline(data = lms[lms$var == "yield_rel_MSY", ], 
              aes(intercept = Intercept, slope = k)) + 
  labs(y = "yield/MSYyield")
p6 <- ggplot(data = glm_data[glm_data$fhist == "one-way", ],
             aes(x = k, y = iav)) +
  geom_point() + theme_bw() + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(intercept = -res_glmnet$iav["(Intercept)",],
              slope = -res_glmnet$iav["k",], col = "red") +
  geom_abline(data = lms[lms$var == "iav", ], 
              aes(intercept = -Intercept, slope = -k)) + 
  labs(y = "catch iav") +
  scale_y_reverse()
par(mfrow = c(2, 3))
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
### save plot
png("output/perfect_knowledge/plots/3.2.1/glmnet/lm.png", 
    width = 15, height = 10, units = "cm", res = 300, type = "cairo")
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
dev.off()


### test: remove penalty for k
penalities <- sapply(var_pred, function(x) ifelse(x == "k", 0, 1))
test_k <- cv.glmnet(
  x = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred]),
  y = as.matrix(glm_data[glm_data$fhist == "one-way", var_resp]),
  foldid = 1:29, nlambda = 1000, family = "mgaussian", alpha = 1, 
  penalty.factor = penalities)
coef(test_k, s = "lambda.min")
coef(test_k, s = "lambda.1se")
min(test_k$cvm)
### same as when lm applied individually
### -> difference between lm and glmnet is because of penalty

### test alpha (lasso vs. ridge)
foldid <- 1:29
alphas <- seq(from = 0, to = 1, 0.01)
names(alphas) <- alphas
fits <- foreach(x = alphas, .errorhandling = "pass") %dopar% {
  cv.glmnet(
    foldid = 1:29, nlambda = 1000, alpha = x, family = "mgaussian",
    x = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred]),
    y = as.matrix(glm_data[glm_data$fhist == "one-way", var_resp]))
}
names(fits) <- alphas
### extract mean-sqared error
mses <- lapply(alphas, function(x) {
  data.frame(lambda = fits[[as.character(x)]]$lambda, 
             mse = fits[[as.character(x)]]$cvm, alpha = x)
})
mses <- do.call(rbind, mses)
ggplot(data = mses[], 
       aes(x = log(lambda), y = mse, colour = alpha,
           group = alpha)) +
  geom_line() + theme_bw() +
  labs(x = "log(Lambda)", y = "Mean-Squared Error")
ggsave(filename = "output/perfect_knowledge/plots/3.2.1/glmnet/alphas_mgauss.png",
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo")
### find alpha with lowest error
a <- mses$alpha[which.min(mses$mse)]
### coefficients
coef(fits[[as.character(a)]])
### mse
min(fits[[as.character(a)]]$cvm)

#best_fit <- fits[[as.character(a)]]
best_fit <- fits[[length(fits)]]
pred <- predict(best_fit, 
                newx = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred]),
                s = "lambda.1se")
### predict for many linf & k
pred_data <- glm_data[glm_data$fhist == "one-way", var_pred][1, ]
pred_data[] <- 0
pred_data[2:100, ] <- 0
pred_data$k <- seq(0, 1, length.out = 100)
pred_values <- predict(best_fit, 
                       newx = as.matrix(pred_data),
                       s = "lambda.1se")
pred_data <- cbind(pred_data, pred_values[,, 1])
### predict for stocks

pred_stocks_vals <- predict(best_fit, 
                            newx = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred]),
                            s = "lambda.1se")
pred_stocks <- cbind(glm_data[glm_data$fhist == "one-way", var_pred], 
                     pred_stocks_vals[,, 1])
### plot
ggplot(data = pred_data[, ],
       aes(x = k, y = yield_rel_MSY)) +
  geom_line(size = 1) +
  theme_bw()
### data points from simulation
sim_data <- scns_res[scns_res$fhist == "one-way",
                     c("linf", "k", "a", "b", "t0", "a50", "yield_rel_MSY",
                       "ssb_below_blim_total")]
df_plot <- rbind(cbind(pred_data, source = "model"), 
                 cbind(sim_data, source = "data"),
                 cbind(pred_stocks, source = "fitted"))
df_plot <- gather(df_plot, key = "stat", value = "value", yield_rel_MSY,
                  ssb_below_blim_total)
ggplot(data = df_plot,
       aes(x = k, y = value, colour = linf, group = linf, fill = linf)) +
  geom_line(data = df_plot[df_plot$source == "model", ], size = 1) +
  theme_bw() +
  facet_grid(~ stat, scales = "free_y") +
  ylim(0, NA) +
  geom_point(data = df_plot[df_plot$source == "fitted", ],
             colour = "black", shape = 24, size = 2,
             aes(x = k, y = value, fill = linf)) +
  geom_point(data = df_plot[df_plot$source == "data", ],
             colour = "black", shape = 21, size = 3,
             aes(x = k, y = value, fill = linf)) +
  geom_line(data = df_plot[df_plot$source %in% c("data", "fitted"), ],
            colour = "black")
### inverse
ggplot(data = df_plot,
       aes(x = linf, y = value, colour = k, group = k, fill = k)) +
  geom_line(data = df_plot[df_plot$source == "model", ], size = 1) +
  theme_bw() +
  facet_grid(~ stat, scales = "free_y") +
  ylim(0, NA) +
  geom_point(data = df_plot[df_plot$source == "fitted", ],
             colour = "black", shape = 24, size = 2,
             aes(x = linf, y = value, fill = k)) +
  geom_point(data = df_plot[df_plot$source == "data", ],
             colour = "black", shape = 21, size = 3,
             aes(x = linf, y = value, fill = k)) +
  geom_line(data = df_plot[df_plot$source %in% c("data", "fitted"), ],
            colour = "black", aes(group = linf))

### ------------------------------------------------------------------------ ###
### extended parameters & iterations ####
### ------------------------------------------------------------------------ ###

stats_iter <- readRDS(file = "output/1237_6804_stats_iter.rds")
### default 3.2.1 scenarios
glm_data_sub <- glm_data[glm_data$fhist == "one-way", ]
### get stats
stats_iter <- stats_iter[as.character(glm_data_sub$scenario)]
stats_iter_df <- lapply(names(stats_iter), function(x) {
  data.frame(scenario = as.numeric(x),
             iter_iav = c(stats_iter[[x]]$iav),
             iter_f_rel = c(stats_iter[[x]]$f_rel),
             iter_ssb_rel = c(stats_iter[[x]]$ssb_rel),
             iter_catch_MSY_prop = c(stats_iter[[x]]$catch_MSY_prop),
             iter_collapse_risk = c(stats_iter[[x]]$collapse_risk))
})
stats_iter_df <- do.call(rbind, stats_iter_df)
### merge with stock parameters
stats_iter_df <- merge(stats_iter_df, glm_data_sub, all = TRUE)

### glmnet
### first: test lasso only
glmnet_iter <- cv.glmnet(
  x = as.matrix(stats_iter_df[, var_pred]),
  y = as.matrix(stats_iter_df[, var_resp]),
  foldid = 1:29, nlambda = 1000, 
  family = "mgaussian", alpha = 1)
coef(test, s = "lambda.min")
coef(test, s = "lambda.1se")

