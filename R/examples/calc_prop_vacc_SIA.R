# Two dose measles vaccine: routine with SIA
calc_prop_vacc_SIA(V=c(0.9, 0.7),
                   S=0.95,
                   effectiveness=c(0.84, 0.941, 0.99),
                   independent=FALSE)

calc_prop_vacc_SIA(V=c(0.9, 0.7),
                   S=0.95,
                   effectiveness=c(0.84, 0.941, 0.99),
                   independent=TRUE)

# Three dose vaccine: routine with SIA
calc_prop_vacc_SIA(V=c(0.9, 0.8, 0.7),
                   S=0.95,
                   effectiveness=c(0.85, 0.9, 0.95, 0.99),
                   independent=FALSE)

calc_prop_vacc_SIA(V=c(0.9, 0.8, 0.7),
                   S=0.95,
                   effectiveness=c(0.85, 0.9, 0.95, 0.99),
                   independent=TRUE)
