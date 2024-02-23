# Repro-R

## function conventions

#### For convenience, I'll list some of the data structures and their names in this file.
'''p_value(lower_bd, upper_bd, t_init, seeds, G, s_obs, T_stat)'''

-'''lower_bd''': a vector containing the lower bounds for parameter theta (same for upper_bd).
-'''t_init''': an initial theta (k-dim) for the optim function to start with.
-'''seeds''': an R by d matrix whose rows are d-dim seeds.
-'''G''': the data generating function $G(u_i, \theta)$.
-'''s_obs''': the observed data (d-dim).
-'''T_stat''': the function that calculates the statistics, its inputs are (in order): an l by d matrix containing l vectors whose statistics we want to compute; an R+1 by d matrix containing R+1 vectors with respect to which the statistics is computed; parameter theta (may or may not be explicited used in the function)