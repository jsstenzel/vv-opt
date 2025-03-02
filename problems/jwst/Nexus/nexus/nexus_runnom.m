%nexus_runnom.m
clear all
close all
tstart=cputime;
nexus_init
nexus_perf
disp(['Total Runtime: ' num2str(cputime-tstart) ' [sec]'])

% Lyapunov simulation results
disp('Lyapunov simulation results')
RMMS_WFE_lyap=z1 % RMMS WFE (root-mean-mean-square wavefront error)
RSS_LOS_lyap=z2 %RSS LOS (root-sum-square line-of-sight)

% time domain simulation results
disp('Time-domain simulation results')
RMMS_WFE_time=std(WFE)
RSS_LOS_time=sqrt(sum(std(LOS).^2))
