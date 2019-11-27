calc_prop_vacc_SIA(v1=0.9, v2=0.8, S=0.97) # High coverage two-dose routine and campaign

calc_prop_vacc_SIA(v1=0.9, v2=0.85, v3=0.8, S=0.95) # High coverage two-dose routine and campaign

calc_prop_vacc_SIA(v1=0.5, v2=0.5, v3=0.5, S=0.95) # Low routine coverage and high campaign

# Boundary conditions
calc_prop_vacc_SIA(v1=1, v2=1, v3=1, S=1)
calc_prop_vacc_SIA(v1=1, v2=0, v3=0, S=0)

# Complex models reduce to simpler models
calc_prop_vacc_SIA(v1=0.9, v2=0.8, v3=0.7, S=0)
calc_prop_vacc(v1=0.9, v2=0.8, v3=0.7)

calc_prop_vacc_SIA(v1=0.9, v2=0, v3=0, S=0)
calc_prop_vacc(v1=0.9, v2=0, v3=0)
calc_prop_vacc(v1=0.9, v2=0)
