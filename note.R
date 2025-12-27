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

###

library(haven)
smoke_prev <- read_dta("data/smoke_prev.dta")

## smoking change in each gender 

### change with age

smoke_prev |> 
  group_by(wave,sex) |> 
  summarise(mean_regsmoke = mean(regsmoke),.groups = "drop") |>
  ggplot(aes(x = wave,y = mean_regsmoke))+
  geom_point()+
  ylim(c(0,1))+
  facet_wrap(~sex)

smoke_prev <- transform(smoke_prev, age_c = age - 16.2,
                        wave_c = (wave - 3.5)/2)

smoke_prev |> 
  ggplot(aes(x = wave_c, y = age_c))+
  geom_line(aes(group = id))

smoke_prev |> 
  ggplot(aes(x = wave, y = age))+
  geom_line(aes(group = id))

mean(smoke_prev$wave)

## centering make the estimation more intepretable, (mean = 0, the negative estimation => below the mean)
## to create a meaningful zero point => make the intercept more meaningful
## run the regression with interaction terms
## important in multilevel modeling 

library(Hmisc)
library(geepack); library(gee) # two alternative packages for fitting gee models
library(doBy)  


glm.fit1 <- glm(regsmoke ~ sex+age_c, data = smoke_prev, family =binomial)
# summary(glm.fit1)

smoke_prev |> 
  mutate(pred = exp(predict(glm.fit1))) |> 
  group_by(age,sex) |> 
  summarise(meanreg = mean(pred),.groups = "drop") |>
  pivot_wider(names_from = sex,values_from = meanreg)
  ggplot(aes(x = age, y = meanreg))+
  geom_line()+
  facet_wrap(~sex)+
  scale_y_continuous(limits = c(0,1))

predict(glm.fit1)

glm.fit1.1 <- glm(regsmoke ~ sex+I(age_c*(sex=="male"))+I(age_c*(sex=="female")), 
                  data = smoke_prev, family =binomial)
summary(glm.fit1.1)


smoke_prev |> 
  group_by(wave,age) |> 
  summarise(mean_regsmoke = mean(regsmoke),.groups = "drop") |>
  ggplot(aes(x = wave,y = mean_regsmoke))+
  geom_point()+
  ylim(c(0,1))+
  facet_wrap(~sex)


### respiratory data
library(geepack)

library(janitor)

respiratory |> 
  group_by(treat,visit) |> 
  summarise(out = mean(outcome),
            base = mean(baseline),.groups = "drop") |> 
  pivot_wider(names_from = visit,values_from = out) |> 
  magrittr::set_colnames(c("treat","0","1","2","3","4")) |> 
  pivot_longer(cols = -treat,names_to = "visit") |> 
  ggplot(aes(x = visit, y = value,color = treat))+
  geom_point()+
  ylim(c(0,1))

respiratory |> 
  ggplot(aes(x = visit, y = outcome,color = treat))+
  # geom_jitter(aes(group = id))
  geom_point(position = position_jitterdodge(0.1, dodge.width = .1))+
  stat_summary(fun = "mean", size = 2, 
               position = position_dodge(0.1), geom = "point")+
  stat_summary(fun.data = "mean_se" ,width = 0, geom = "errorbar", 
               position = position_dodge(0.1)) 

respiratory |> 
  pivot_wider(id_cols = c(baseline,visit,outcome))

?respiratory

fit <- geeglm(outcome ~ treat*visit, data=respiratory, id=id, 
       family=binomial(), corstr="unstructured") 

summary(fit)

factor(smoke_prev$sex)

library(dplyr)
library(magrittr)
smoke_prev %<>% 
  mutate(sex = factor(sex,labels = c("female","male"),levels = c(1,0))) 

smoke_prev |> 
  mutate(sex = factor)

attr(smoke_prev$sex, "label")

library(gtsummary)

glm.fit1 |> 
  tbl_regression(exponentiate = TRUE)
