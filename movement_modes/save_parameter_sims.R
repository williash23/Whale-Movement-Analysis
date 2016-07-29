# Sara Williams
# 2/22/2016; 
# Save parameter simulation iterations from models formulated after code from Morales et al. 2004
#   and run in JAGS.
################################################################################



#  Save "single" model ouput


sim_reps_v1_doub <- out_double_dcat$BUGS$sims.list$v[,1]
sim_reps_v2_doub <- out_double_dcat$BUGS$sims.list$v[,2]
sim_reps_scale1_doub <- out_double_dcat$BUGS$sims.list$scale[,1]
sim_reps_scale2_doub <- out_double_dcat$BUGS$sims.list$scale[,2]
sim_reps_mu1_doub <- out_double_dcat$BUGS$sims.list$mu[,1]
sim_reps_mu2_doub <- out_double_dcat$BUGS$sims.list$mu[,2]
sim_reps_mean.q1_doub <- out_double_dcat$BUGS$sims.list$mean.q[,1]
sim_reps_mean.q2_doub <- out_double_dcat$BUGS$sims.list$mean.q[,2]  

