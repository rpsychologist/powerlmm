# This recreates TABLE 2 in Raudenbush and Liu (2001)


library(dplyr)
library(tidyr)
p <- study_parameters(n1 = 6,
                      n2 = seq(100,800, by=100)/2,
                      T_end = c(5),
                      sigma_subject_intercept = 4,
                      sigma_subject_slope = sqrt(0.0030),
                      sigma_error = sqrt(0.0262),
                      effect_size = sqrt(0.0030)*0.4 * 5)
x <- get_power(p)
round(x$power, 2)
x <- as.data.frame(x)

x %>%
    group_by(T_end) %>%
    select(n1, n2, power) %>%
    spread(n2, power) %>%
    as.data.frame %>%
    print(digits = 5)

get_slope_diff(p)


time <- c(-2, -1, 0, 1, 2)

sum((time - mean(time))^2)

time2 <- poly(time)
sum( time2^2)

n1 <- 5
(n1+1)*n1*(n1-1)/12


time <- seq(0, 4, length.out = 10)
n1 <- length(time)
sx <- sum( (time - mean(time))^2)/n1


error <- 1.2
u1 <- 0.3
n2 <- 30

(error^2 + n1*u1^2*sx) / (n1*n2*sx)

sx2 <- sum( (time - mean(time))^2)
(error^2 + u1^2*sx2) / (n2*sx2)



sum( (poly(time) - mean(poly(time)))^2)


## TAble 1
# Column 2
p <- study_parameters(n1 = 3,
                      n2 = 238/2,
                      T_end = c(2),
                      sigma_subject_intercept = 4,
                      sigma_subject_slope = sqrt(0.0030),
                      sigma_error = sqrt(0.0262),
                      effect_size = sqrt(0.0030)*0.4 * 2)
x <- get_power(p)
round(x$power, 2)


p <- study_parameters(n1 = 4,
                      n2 = 238/2,
                      T_end = c(3),
                      sigma_subject_intercept = 4,
                      sigma_subject_slope = sqrt(0.0030),
                      sigma_error = sqrt(0.0262),
                      effect_size = sqrt(0.0030) * 0.4 * 3)
x <- get_power(p)
round(x$power, 2)

p <- study_parameters(n1 = 5,
                      n2 = 238/2,
                      T_end = c(4),
                      sigma_subject_intercept = 4,
                      sigma_subject_slope = sqrt(0.0030),
                      sigma_error = sqrt(0.0262),
                      effect_size = sqrt(0.0030) * 0.4 * 4)
x <- get_power(p)
round(x$power, 2)


## Column 3
f <- 2
D <- 2
n1 <- f * D +  1
p <- study_parameters(n1 = n1,
                      n2 = 238/2,
                      T_end = D,
                      sigma_subject_intercept = 4,
                      sigma_subject_slope = sqrt(0.0030),
                      sigma_error = sqrt(0.0262),
                      effect_size = sqrt(0.0030) * 0.4 * D)
x <- get_power(p)
round(x$power, 2)

## For this ES, don't make it a post-test ES.
