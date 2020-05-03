# dhs_namibia_2013 vaccination coverage file downloaded from the website
data("dhs_namibia_2013_namibia_2013")

# For all observations together (population mean)
tmp <- get_beta_params_age(age=dhs_namibia_2013$age.in.months,
                           vacc=dhs_namibia_2013$measles.y,
                           breaks=NULL)


# Two age groups < 12 months and > 12 months
tmp <- get_beta_params_age(age=dhs_namibia_2013$age.in.months,
                           vacc=dhs_namibia_2013$measles.y,
                           breaks=12)
par(mfrow=c(1,2))
for(i in 1:2) curve(dbeta(x, tmp$shape1[i], tmp$shape2[i]), 0, 1,
                    xlab='Proportion vaccinated', ylab='Density',
                    main=tmp$age[i], lwd=2)


# 6-month age groups
tmp <- get_beta_params_age(age=dhs_namibia_2013$age.in.months,
                           vacc=dhs_namibia_2013$measles.y,
                           breaks=seq(0, 60, 6))
par(mfrow=c(2,5))
for(i in 1:10) curve(dbeta(x, tmp$shape1[i], tmp$shape2[i]), 0, 1,
                     xlab='Proportion vaccinated', ylab='Density',
                     main=tmp$age[i], lwd=2)


# Each unique age in months
tmp <- get_beta_params_age(age=dhs_namibia_2013$age.in.months,
                           vacc=dhs_namibia_2013$measles.y,
                           breaks=seq(0, 60, 1))

par(mfrow=c(1,1))
plot(tmp$mu, type='l')


# 3-month age groups for each region
require(foreach)
tmp <- foreach(i=unique(dhs_namibia_2013$region.residence), .combine='rbind') %do% {

  sel <- dhs_namibia_2013$region.residence == i

  cbind(
    data.frame(region=i),
    get_beta_params_age(age=dhs_namibia_2013$age.in.months[sel],
                        vacc=dhs_namibia_2013$measles.y[sel],
                        breaks=seq(0, 60, 3))
  )
}
