chapmanRichards <- function(y, k, t, p) {
  y <- max(y)*(1 - exp(-k*t))^p
}

plot(1:200, chapmanRichards(1:200*15, 0.025, 1:200, 3))

a$totMerch

nlsout <- nlrob(totMerch ~ A * (1 - exp(-k * ageNumeric))^p,
                data = a,
                start = list(A = 83, k = 0.03, p = 4),
                trace = TRUE)
summary(nlsout)
summary(nlsout)$coefficients[1 : 3]

plot(fitted(nlsout))
