library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
install.packages("mice")

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
