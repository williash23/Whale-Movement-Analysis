# Sara Williams
# 2/22/2016; 
# Save parameter simulation reps from models formulated after code from Morales et al. 2004
#   and run in JAGS.
################################################################################

#  Save "single" model ouput
save(out_single, file="out_single.RData")

sim_reps_b0_sing <- out_single$BUGS$sims.list$b0
sim_reps_a0_sing <- out_single$BUGS$sims.list$a0
sim_reps_mu0_sing <- out_single$BUGS$sims.list$mu0
    
save(sim_reps_b0_sing, file="sim_reps_b0_sing.RData")
save(sim_reps_a0_sing, file="sim_reps_a0_sing.RData")
save(sim_reps_mu0_sing, file="sim_reps_mu0_sing.RData")

#  Save "double" model ouput
save(out_double_dcat, file="out_double_dcat.RData")

sim_reps_v1_doub <- out_double_dcat$BUGS$sims.list$v[,1]
sim_reps_v2_doub <- out_double_dcat$BUGS$sims.list$v[,2]
sim_reps_scale1_doub <- out_double_dcat$BUGS$sims.list$scale[,1]
sim_reps_scale2_doub <- out_double_dcat$BUGS$sims.list$scale[,2]
sim_reps_mu1_doub <- out_double_dcat$BUGS$sims.list$mu[,1]
sim_reps_mu2_doub <- out_double_dcat$BUGS$sims.list$mu[,2]
sim_reps_mean.q1_doub <- out_double_dcat$BUGS$sims.list$mean.q[,1]
sim_reps_mean.q2_doub <- out_double_dcat$BUGS$sims.list$mean.q[,2]  

save(sim_reps_v1_doub, file="sim_reps_v1_doub.RData")
save(sim_reps_v2_doub, file="sim_reps_v2_doub.RData")
save(sim_reps_scale1_doub, file="sim_reps_scale1_doub.RData")
save(sim_reps_scale2_doub, file="sim_reps_scale2_doub.RData")
save(sim_reps_mu1_doub, file="sim_reps_mu1_doub.RData")
save(sim_reps_mu2_doub, file="sim_reps_mu2_doub.RData")
save(sim_reps_mean.q1_doub, file="sim_reps_mean.q1_doub.RData")
save(sim_reps_mean.q2_doub, file="sim_reps_mean.q2_doub.RData")

#  Save "double switch" model ouput
save(out_double_sw_dcat, file="out_doub_sw_dcat.RData")

sim_reps_v1_doub_sw <- out_double_sw_dcat$BUGS$sims.list$v[,1]
sim_reps_v2_doub_sw <- out_double_sw_dcat$BUGS$sims.list$v[,2]
sim_reps_scale1_doub_sw <- out_double_sw_dcat$BUGS$sims.list$scale[,1]
sim_reps_scale2_doub_sw <- out_double_sw_dcat$BUGS$sims.list$scale[,2]
sim_reps_mu1_doub_sw <- out_double_sw_dcat$BUGS$sims.list$mu[,1]
sim_reps_mu2_doub_sw <-out_double_sw_dcat$BUGS$sims.list$mu[,1]
sim_reps_mean.q1_doub_sw <- out_double_sw_dcat$BUGS$sims.list$mean.q[,1]
sim_reps_mean.q2_doub_sw <- out_double_sw_dcat$BUGS$sims.list$mean.q[,2]
sim_reps_beta01_doub_sw <- out_double_sw_dcat$BUGS$sims.list$beta0[,1]
sim_reps_beta02_doub_sw <- out_double_sw_dcat$BUGS$sims.list$beta0[,2]
    
save(sim_reps_v1_doub_sw, file="sim_reps_v1_doub_sw.RData")
save(sim_reps_v2_doub_sw, file="sim_reps_v2_doub_sw.RData")
save(sim_reps_scale1_doub_sw, file="sim_reps_scale1_doub_sw.RData")
save(sim_reps_scale2_doub_sw, file="sim_reps_scale2_doub_sw.RData")
save(sim_reps_mu1_doub_sw, file="sim_reps_mu1_doub_sw.RData")
save(sim_reps_mu2_doub_sw, file="sim_reps_mu2_doub_sw.RData")
save(sim_reps_mean.q1_doub_sw, file="sim_reps_mean.q1_doub_sw.RData")
save(sim_reps_beta01_doub_sw, file="sim_reps_beta01_doub_sw.RData")
save(sim_reps_beta02_doub_sw, file="sim_reps_beta02_doub_sw.RData")

#  Save "double with covariate " model ouput
save(out_doub_cov_dcat, file="out_doub_cov_dcat.RData")

sim_reps_v1_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$v[,1]
sim_reps_v2_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$v[,2]
sim_reps_scale1_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$scale[,1]
sim_reps_scale2_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$scale[,2]
sim_reps_mu1_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$mu[,1]
sim_reps_mu2_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$mu[,2]
sim_reps_beta0_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$beta0
sim_reps_beta1_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$beta1
sim_reps_mean.q1_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$mean.q[,1]
sim_reps_mean.q2_doub_cov <- out_doub_cov_dcat$BUGS$sims.list$mean.q[,2]

save(sim_reps_v1_doub_cov, file="sim_reps_v1_doub_cov.RData")  
save(sim_reps_v2_doub_cov, file="sim_reps_v2_doub_cov.RData")    
save(sim_reps_scale1_doub_cov, file="sim_reps_scale1_doub_cov.RData")
save(sim_reps_scale2_doub_cov, file="sim_reps_scale2_doub_cov.RData")
save(sim_reps_mu1_doub_cov, file="sim_reps_mu1_doub_cov.RData")
save(sim_reps_mu2_doub_cov, file="sim_reps_mu2_doub_cov.RData")
save(sim_reps_beta0_doub_cov, file="sim_reps_beta0_doub_cov.RData")
save(sim_reps_beta1_doub_cov, file="sim_reps_beta1_doub_cov.RData")
save(sim_reps_mean.q1_doub_cov, file="sim_reps_mean.q1_doub_cov.RData")
save(sim_reps_mean.q2_doub_cov, file="sim_reps_mean.q2_doub_cov.RData")

#  Save "double switch with covariate" model ouput
save(out_double_sw_cov_dcat, file="out_double_sw_cov_dcat.RData")

sim_reps_v1_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$v[,1]
sim_reps_v2_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$v[,2]
sim_reps_scale1_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$scale[,1]
sim_reps_scale2_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$scale[,2]
sim_reps_mu1_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$mu[,1]
sim_reps_mu2_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$mu[,2]
sim_reps_mean.q1_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$mean.q[,1]
sim_reps_mean.q2_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$mean.q[,2]
sim_reps_beta01_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$beta0[,1]
sim_reps_beta02_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$beta0[,2]
sim_reps_beta11_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$beta1[,1]
sim_reps_beta12_doub_sw_cov <- out_double_sw_cov_dcat$BUGS$sims.list$beta1[,2]
    
save(sim_reps_v1_doub_sw_cov, file="sim_reps_v1_doub_sw_cov.RData")
save(sim_reps_v2_doub_sw_cov, file="sim_reps_v2_doub_sw_cov.RData")
save(sim_reps_scale1_doub_sw_cov, file="sim_reps_scale1_doub_sw_cov.RData")
save(sim_reps_scale2_doub_sw_cov, file="sim_reps_scale2_doub_sw_cov.RData")
save(sim_reps_mu1_doub_sw_cov, file="sim_reps_mu1_doub_sw_cov.RData")
save(sim_reps_mu2_doub_sw_cov, file="sim_reps_mu2_doub_sw_cov.RData")
save(sim_reps_beta01_doub_sw_cov, file="sim_reps_beta01_doub_sw_cov.RData")
save(sim_reps_beta02_doub_sw_cov, file="sim_reps_beta02_doub_sw_cov.RData")
save(sim_reps_beta11_doub_sw_cov, file="sim_reps_beta11_doub_sw_cov.RData")
save(sim_reps_beta12_doub_sw_cov, file="sim_reps_beta12_doub_sw_cov.RData")
save(sim_reps_mean.q1_doub_sw_cov, file="sim_reps_mean.q1_doub_sw_cov.RData")
save(sim_reps_mean.q2_doub_sw_cov, file="sim_reps_mean.q2_doub_sw_cov.RData")

#  Save "triple switch" model ouput
save(out_triple_sw_dcat, file="out_triple_sw_dcat.RData")

sim_reps_v1_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$v[,1]
sim_reps_v2_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$v[,2]
sim_reps_v3_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$v[,3]
sim_reps_scale1_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$scale[,1]
sim_reps_scale2_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$scale[,2]
sim_reps_scale3_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$scale[,3]
sim_reps_mu1_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$mu[,1]
sim_reps_mu2_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$mu[,2]
sim_reps_mu3_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$mu[,3]
sim_reps_mean.q1_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$mean.q[,1]
sim_reps_mean.q2_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$mean.q[,2]
sim_reps_mean.q3_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$mean.q[,3]
sim_reps_beta01_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$beta0[,1]
sim_reps_beta02_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$beta0[,2]
sim_reps_beta03_triple_sw_dcat <- out_triple_sw_dcat$BUGS$sims.list$beta0[,3]

save(sim_reps_v1_triple_sw_dcat, file="sim_reps_v1_triple_sw_dcat.RData")
save(sim_reps_v2_triple_sw_dcat, file="sim_reps_v2_triple_sw_dcat.RData")
save(sim_reps_v3_triple_sw_dcat, file="sim_reps_v3_triple_sw_dcat.RData")
save(sim_reps_scale1_triple_sw_dcat, file="sim_reps_scale1_triple_sw_dcat.RData")
save(sim_reps_scale2_triple_sw_dcat, file="sim_reps_scale2_triple_sw_dcat.RData")
save(sim_reps_scale3_triple_sw_dcat, file="sim_reps_scale3_triple_sw_dcat.RData")
save(sim_reps_mu1_triple_sw_dcat, file="sim_reps_mu1_triple_sw_dcat.RData")
save(sim_reps_mu2_triple_sw_dcat, file="sim_reps_mu2_triple_sw_dcat.RData")
save(sim_reps_mu3_triple_sw_dcat, file="sim_reps_mu3_triple_sw_dcat.RData")
save(sim_reps_beta01_triple_sw_dcat, file="sim_reps_beta01_triple_sw_dcat.RData")
save(sim_reps_beta02_triple_sw_dcat, file="sim_reps_beta02_triple_sw_dcat.RData")
save(sim_reps_beta03_triple_sw_dcat, file="sim_reps_beta03_triple_sw_dcat.RData")
save(sim_reps_mean.q1_triple_sw_dcat, file="sim_reps_mean.q1_triple_sw_dcat.RData")
save(sim_reps_mean.q2_triple_sw_dcat, file="sim_reps_mean.q2_triple_sw_dcat.RData")
save(sim_reps_mean.q3_triple_sw_dcat, file="sim_reps_mean.q3_triple_sw_dcat.RData")