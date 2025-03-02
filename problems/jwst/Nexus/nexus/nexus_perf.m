% nexus_perf.m
% NEXUS closed loop performance evaluation
% Comparison of Lyapunov and Time Domain Simulation
% original June 10,2001 , dWo
% updated June 1, 2005, dWo

%-------------------------------------------------------------
% Solve Lyapunov Equation
%-------------------------------------------------------------

warning off
if diagnostics
   warning on
   tstartnewlyap=cputime;
end
syszd=ss(Azd,Bzd,Czd,Dzd); scf=0;
[syszdm,T]=ss2mod7(syszd,scf);
blocksize=20;
[Azdm,Bzdm,Czdm,Dzdm]=ssdata(syszdm);
Sqm = newlyap(Azdm,Bzdm,blocksize);
Szm=Czdm*Sqm*Czdm';
z1=sqrt((1/nray)*sum(diag(Szm(1:nray,1:nray)))); %RMMS WFE
z2=sqrt(sum(diag((Szm(nray+1:nray+2,nray+1:nray+2))))); %RSS LOS
if diagnostics
   disp(['Performance Evaluation (newlyap) took: ' num2str(cputime-tstartnewlyap) ' seconds'])
   disp(['Performance Results: z1= ' num2str(z1) ' [nm], z2= ' num2str(z2) ' [\mum]'])
end
warning on

%-------------------------------------------------------------
% Automatically Execute Time-Domain Simulation
%-------------------------------------------------------------

% Time simulation parameters

dt=0.001;   % 2000 Hz simulation sample frequency = 5*~411.3 Hz highest freq in model
Seed=round(1e5*rand(10,1));  % random initial seeds for white noise sources
Tmax=5;    % Time simulation length
nray=size(dwdu,1);  % number of rays in dynamics model
ndof=size(Dp,1)-6; dwdu=dwdu(1:nray,1:ndof);

nexus_block  % calls SIMULINK block diagram
sim('nexus_block')
