#!/usr/bin/env Rscript

library("ipumsr")
sessionInfo()

dir.create("data-clean", showWarnings = FALSE)
setwd("data-raw")
raw <- read_ipums_micro(read_ipums_ddi("cps_00001.xml"))

`%notin%` <- Negate(`%in%`)
tails <- function(x) {
  p <- c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.25)
  quantile(x, probs = c(p, 0.5, 1 - rev(p)))
}

###################
## Data cleaning ##
###################
d <- raw
d <- d[d$YEAR >= 1994, ] # According to the codebook, UHRSWORKORG is messed up until 1994.
dim(d)

# Location
head(d$STATECENSUS)
tb <- table(d$STATECENSUS)
length(tb) # 50 states, plus D.C.
sort(tb)
rm(tb)
d$DIVISION <- d$STATECENSUS %/% 10
sort(table(d$DIVISION))
d$REGION <- numeric(length(d$DIVISION))
d$REGION[d$DIVISION %in% c(1, 2)] <- 1
d$REGION[d$DIVISION %in% c(3, 4)] <- 2
d$REGION[d$DIVISION %in% c(5, 6, 7)] <- 3
d$REGION[d$DIVISION %in% c(8, 9)] <- 4
table(d$REGION)
cd <- attr(d$STATECENSUS, "labels")
d$STATENAME <- names(cd)[match(d$STATECENSUS, cd)]
rm(cd)
sort(table(d$STATENAME))

# Time
table(d$YEAR)

# Other Covariates
head(d$AGE)
summary(d$AGE)

head(raw$SEX)
table(raw$SEX)
d$MALE <- 1 * (d$SEX == 1)
table(d$MALE)
d$SEX <- NULL

head(d$RACE)
table(d$RACE)
d$WHITE <- 1 * (d$RACE == 100)
table(d$WHITE)
mean(d$WHITE) # GBD: Seems too high.
d$RACE <- NULL

head(d$HISPAN)
table(d$HISPAN)
d$HISPANIC <- 1 * (d$HISPAN != 0)
table(d$HISPANIC)
d$HISPAN <- NULL

head(d$EDUC)
table(d$EDUC)
d$HSDROPOUT <- 1 * (d$EDUC < 72)
d$HSGRAD <- 1 * (d$EDUC %in% c(72, 73))
d$SOMECOL <- 1 * (d$EDUC > 73 & d$EDUC <= 109)
d$BACHELORS <- 1 * (d$EDUC %in% c(110, 111))
d$ADVDEGREE <- 1 * (d$EDUC > 111)
table(d$HSDROPOUT)
table(d$HSGRAD)
table(d$SOMECOL)
table(d$BACHELORS)
table(d$ADVDEGREE)
sum(d$HSDROPOUT + d$HSGRAD + d$SOMECOL + d$BACHELORS + d$ADVDEGREE)
nrow(d)
d$EDUC <- NULL

# Skip METRO because it doesn't seem to have a clear interpretation.
head(d$METRO)
table(d$METRO)
d$METRO <- NULL

head(d$UNION)
table(d$UNION)
d$UNION <- 1 * (d$UNION %in% c(2, 3))
table(d$UNION)

head(d$MARST)
table(d$MARST)
d$MARRIED <- 1 * (d$MARST %in% c(1, 2))
table(d$MARRIED)
d$MARST <- NULL

# Response 1: Employed (0/1 variable)
head(d$EMPSTAT)
table(d$EMPSTAT)
d$EMPLOYED <- 1 * (d$EMPSTAT %in% c(1, 10, 12))
table(d$EMPLOYED)
d$EMPSTAT <- NULL

# Fork the datasets before filtering that requires employment.
e <- d
e$EARNWEEK <- NULL
e$UHRSWORKORG <- NULL
names(e)
d$EMPLOYED <- NULL

# Response 2: Usual hours worked per week (continuous)
head(d$UHRSWORKORG)
d <- d[d$UHRSWORKORG %notin% c(998, 999), ]
dim(d)
tails(d$UHRSWORKORG)
d <- d[d$UHRSWORKORG >= quantile(d$UHRSWORKORG, 0.01) & d$UHRSWORKORG <= quantile(d$UHRSWORKORG, 0.99), ]
dim(d)
tails(d$UHRSWORKORG)
d$WEEKLYHRS <- d$UHRSWORKORG
d$UHRSWORKORG <- NULL

# Response 3: Average hourly earnings (continuous)
head(d$EARNWEEK)
d <- d[d$EARNWEEK != 9999.99, ]
dim(d)
years <- seq(min(d$YEAR), max(d$YEAR))
earns_week_top_code <- cbind(YEAR = years, EARNWEEKTOP = c(1923.00, 2884.61)[1 + (years >= 1998)])
d <- merge(d, earns_week_top_code)
rm(years, earns_week_top_code)
d$EARNWEEK.is.top.coded <- d$EARNWEEK == d$EARNWEEKTOP
x <- aggregate(EARNWEEK.is.top.coded ~ YEAR, data = d, mean)
x
max_percentage_top_coded <- ceiling(max(x$EARNWEEK.is.top.coded) * 100) / 100
max_percentage_top_coded
d <- d[d$EARNWEEK < quantile(d$EARNWEEK, 1 - max_percentage_top_coded), ]
dim(d)
sum(d$EARNWEEK.is.top.coded)
d$EARNWEEK.is.top.coded <- NULL
d$EARNWEEKTOP <- NULL
rm(x, max_percentage_top_coded)
cpi <- read.csv("inflation.csv")
cpi$INFADJ <- cpi$CPI[cpi$YEAR == 2020] / cpi$CPI
cpi$CPI <- NULL
d <- merge(d, cpi)
d$EARNWEEK2020 <- d$EARNWEEK * d$INFADJ
rm(cpi)
d$AHE <- d$EARNWEEK2020 / d$WEEKLYHRS
tails(d$AHE)
d <- d[d$AHE >= quantile(d$AHE, 0.005) & d$AHE <= quantile(d$AHE, 0.995), ]
tails(d$AHE)
d$EARNWEEK2020 <- NULL
d$EARNWEEK <- NULL
d$INFADJ <- NULL

# Final variables

names(e)
names(d)
setdiff(names(e), names(d))
setdiff(names(d), names(e))
intersect(names(e), names(d))

# Sanity checks

dim(d)

summary(lm(
  data = d,
  log(AHE) ~ AGE + I(AGE^2) + MALE + WHITE + HISPANIC + UNION + HSDROPOUT + SOMECOL + BACHELORS + ADVDEGREE
))

table(d$YEAR)

sapply(1994:2020, function(i) {
  coef(lm(
    data = d[d$YEAR == i, ],
    log(AHE) ~ AGE + I(AGE^2) + MALE + WHITE + HISPANIC + UNION + HSDROPOUT + SOMECOL + BACHELORS + ADVDEGREE
  ))
})

table(d$DIVISION)

sapply(1:9, function(i) {
  coef(lm(
    data = d[d$DIVISION == i, ],
    log(AHE) ~ AGE + I(AGE^2) + MALE + WHITE + HISPANIC + UNION + HSDROPOUT + SOMECOL + BACHELORS + ADVDEGREE
  ))
})

table(d$REGION)

sapply(1:4, function(i) {
  coef(lm(
    data = d[d$REGION == i, ],
    log(AHE) ~ AGE + I(AGE^2) + MALE + WHITE + HISPANIC + UNION + HSDROPOUT + SOMECOL + BACHELORS + ADVDEGREE
  ))
})


# Save the data

saveRDS(e, "../data-clean/d1.rds", compress = "xz")
saveRDS(d, "../data-clean/d2.rds", compress = "xz")

sessionInfo()
