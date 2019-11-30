# Two dose measles vaccine: routine without SIA
calc_prop_vacc(V=c(0.9, 0.8),
               effectiveness=c(0.84, 0.941),
               independent=FALSE)

calc_prop_vacc(V=c(0.9, 0.8),
               effectiveness=c(0.84, 0.941),
               independent=TRUE)

# Three dose vaccine: routine without SIA
calc_prop_vacc(V=c(0.9, 0.8, 0.7),
               effectiveness=c(0.85, 0.9, 0.95))

