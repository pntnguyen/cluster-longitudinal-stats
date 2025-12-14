library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
# install.packages("mice")

library(mice)

potthoffroy |> 
  pivot_longer(cols = -c(id,sex)) |> 
  mutate(age = str_extract(name,'\\d+') |> as.numeric()) |> 
  ggplot(aes(x = age,y = value))+
  geom_line(aes(group = id))+
  facet_wrap(~sex)+
  theme_bw()

# install.packages("ggcorrplot")
library(ggcorrplot)

## correlation matrix
potthoffroy |> 
  filter(sex == "F") |> 
  select(-c(id,sex)) |> 
  cor(use="pairwise.complete.obs")|> 
  ggcorrplot(lab = TRUE) +
  labs(tag = "Female")

## correlation matrix
potthoffroy |> 
  filter(sex == "M") |>
  select(-c(id,sex)) |> 
  cor(use="pairwise.complete.obs") |> 
  ggcorrplot(lab = TRUE) +
  labs(tag = "Male")

data <- potthoffroy |> 
  pivot_longer(cols = -c(id,sex)) |> 
  mutate(age = str_extract(name,'\\d+') |> as.numeric(),
         sex = factor(sex,labels = c(0,1))) 

lm_model <- lm(value ~ sex + age, data = data)

summary(lm_model)
library(gtsummary)
lm_model |> 
tbl_regression(
  intercept = TRUE,
  estimate_fun = purrr::partial(style_ratio, digits = 3)
) |> 
  remove_row_type(
    type = "reference"
  )

### The estimation of beta (sex,age,intercept) are unbiased, 
### but the standard errors are uncorrect if the covariance is not independent 
library(lmtest)
library(sandwich)

model <- lm(value ~ sex + age, data = data)
summary(model)
# cluster-robust SE by subj
cov_cl <- vcovCL(model, cluster = data$id)

coeftest(model, vcov = cov_cl) |>  
  tbl_regression(
    intercept = TRUE,
    estimate_fun = purrr::partial(style_ratio, digits = 3)
  ) 

cov_cl |> 
  cor(use="pairwise.complete.obs")|> 
  ggcorrplot(lab = TRUE)

lm_model

sqrt(diag(vcov(lm_model)))[2]

se_robust_model[, "Std. Error"][2]

sqrt(diag(vcov(model)))[2]

## step-by-step
X  <- model.matrix(model)          # design matrix (n × p)
e  <- resid(model)                 # residuals (n × 1)
id <- data$id             

bread <- solve(t(X) %*% X)

clusters <- unique(id)

G <- length(clusters)

meat <- matrix(0, ncol = ncol(X), nrow = ncol(X))

for (g in clusters) {
  g = 1
  idx <- which(id == g)
  Xg  <- X[idx, , drop = FALSE]
  eg  <- e[idx]
  
  # cluster score = sum over cluster of x_i * e_i
  score_g <- t(Xg) %*% eg
  
  # add outer product to meat
  meat <- meat + score_g %*% t(score_g)
}

bread %*% meat %*% bread

appendix <- appendix |> 
  mutate(lognchanges = log(nchanges),
         group = factor(group,labels = c("Lap","Open")))  

  ggplot(aes(lognchanges))+
  geom_histogram()

appendix |> 
  ggplot(aes(x = day, y = lognchanges))+
  geom_line(aes(group = patid))+
  geom_point()+
  facet_wrap(~group)+
  theme_bw()

## lmer model

### random intercept model
library(lme4)

data <- potthoffroy |> 
  pivot_longer(cols = -c(id,sex)) |> 
  mutate(age = str_extract(name,'\\d+') |> as.numeric(),
         sex = factor(sex,labels = c(0,1)),
         time = age-8,
         y = value) 


random_intecept_model <- lmer(y ~ time + (1 | id),data = data, subset = (sex == "1"))

summary(random_intecept_model)

data |> 
  filter(sex == 1) |> 
  mutate(fitted = predict(random_intecept_model)) |> 
  ggplot(aes(x = time, y = fitted))+
  geom_line(aes(group = id))+
  geom_point(aes(y = value))


## random slope and intecept

random_slope_intecept_model <- lmer(y ~ time + (1 + time | id),data = data, 
                                    subset = (sex == "Male"), na.action=na.exclude,REML=T)

summary(random_slope_intecept_model)

data |> 
  filter(sex == 1) |> 
  mutate(fitted = predict(random_slope_intecept_model)) |> 
  ggplot(aes(x = time, y = fitted))+
  geom_line(aes(group = id))+
  geom_point(aes(y = value))

library(gtsummary)

lmer_gt_tbl <- function(mod){
  
 var_s <<- as.data.frame(VarCorr(mod))$vcov[1]
 var_re <<- as.data.frame(VarCorr(mod))$vcov[length(as.data.frame(VarCorr(mod))$vcov)]
 
  mod |> 
    tbl_regression(
      intercept = TRUE,
      tidy_fun = broom.mixed::tidy,
      estimate_fun = purrr::partial(style_ratio, digits = 3)
    ) |>
    as_gt() |>
    gt::grand_summary_rows(
      columns = everything(),
      fns = list(
        "Variance (subject)"  = ~ sprintf("%.3f", var_s),
        "Variance (residual)" = ~ sprintf("%.3f", var_re)
      )
    )  
}

lmer_gt_tbl(random_intecept_model)

lmer_gt_tbl(random_slope_intecept_model)


as.data.frame(VarCorr(random_intecept_model))$vcov

as.data.frame(VarCorr(random_slope_intecept_model))$vcov |> length()


random_slope_intecept_model2 <- lmer(y ~ sex*time + (1 + time | id),data = data)

summary(random_slope_intecept_model2)
  
data |> 
  mutate(fitted = predict(random_slope_intecept_model2),
         sex2 = factor(sex) |> as.numeric(),
         pred.pop = predict(
           random_slope_intecept_model2,
           re.form = NA   # fixed effects only
         )) |> 
  ggplot(aes(x = time))+
  geom_line(aes(y = pred.pop,group = id),size = 1)+
  # geom_line(aes(y = fitted,group = id),alpha = .5)+
  geom_line(aes(y = value,group = id),alpha = .5)+
  geom_point(aes(y = value))+
  facet_wrap(~sex)




