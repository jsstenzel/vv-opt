% nexus_plant.m
% A MATLAB m-file, which runs the nexus_fem.m file
% and assembles m and k, then solves the eigenproblem
% and builds a state space system of the plant in modal
% form.
% Author: Olivier L. de Weck
% (c) Massachusetts Institute of Technology
% Date: 1 November 2000
% Modified: 3 June 2001
if diagnostics
   tic
end
%tic
nexus_fem
   %disp(['Assembled NEXUS FEM in ' num2str(toc) ' [sec]'])

if diagnostics
   disp(['Assembled NEXUS FEM in ' num2str(toc) ' [sec]'])
end
if diagnostics
   tic
end
%tic
[phi,omeg]=eigfem(k,m);
   %disp(['Solved Eigenproblem in ' num2str(toc) ' [sec]'])

omeg=abs(omeg);
phi=real(phi);
%tic
phi=normphi(phi,m); % mass normalize
   %disp(['onormphi in ' num2str(toc) ' [sec]'])

if diagnostics
   disp(['Solved Eigenproblem in ' num2str(toc) ' [sec]'])
   end

% expand phi back out to full set for future use;
% phifull=zeros(ndof,size(phi,2));
% phifull(nset,:)=phi;
% phifull(mapdof(:,1),:)=phi(mapdof(:,2),:);
% phi=phifull;
% clear phifull
 phi=mce_dis(nset,mset,rg,phi);
% alternative but much slower procedure below
% phig=gset(phi,bc);
% phig(mset,:)=gm*phig(nset,:);
% phi=phig;
if diagnostics
   disp(['Expanded phi matrix out to full set of dof in ' num2str(toc) ' [sec]']);
   end

% replace rigid body modes with geometrically computed ones
nrbm=6;
phi(:,1:nrbm)=phirb;
omeg(1:nrbm,1)=zeros(nrbm,1);
% Fix SM actuator problem
phi(reshape(bci(ids([14 15 16 17 18 19 20 21 23 24 25 26 39 40 41 42 37 38 33 34 35 36 31 32 117 118 116 121 120 119]),:)',180,1),:)=phi(repmat(bci(ids(13),:)',30,1),:);

% select modes for analysis
% nm=119;
nm=79; % number of modes in FEM: results in ~500 Hz as highest flexible mode
omeg=omeg(1:nm+3);
phi=phi(:,1:nm+3);
if diagnostics
   tic
   end
%  (1) build state space model using mode2ss2
%  define input dofs
%  input nodes:  
%  79,80,81,82  Individual RWA locations
%  83  center of RWA assembly for ACS torques
%  207 instrument node for cryocooler (Cryo z-axis colinear with +y basic)
ig=[...
       79 1 1 1 1 1 1
       80 1 1 1 1 1 1 
       81 1 1 1 1 1 1
       82 1 1 1 1 1 1
      207 1 1 1 0 0 0
       83 0 0 0 1 1 1
       ];
       % ig= input vector
       % INPUT VECTOR
       %  1 RW1-Fx
       %  2 RW1-Fy
       %  3 RW1-Fz
       %  4 RW1-Mx
       %  5 RW1-My
       %  6 RW1-Mz
       %  7 RW2-Fx
       %  8 RW2-Fy
       %  9 RW2-Fz
       % 10 RW2-Mx
       % 11 RW2-My
       % 12 RW2-Mz
       % 13 RW3-Fx
       % 14 RW3-Fy
       % 15 RW3-Fz
       % 16 RW3-Mx
       % 17 RW3-My
       % 18 RW3-Mz
       % 19 RW4-Fx
       % 20 RW4-Fy
       % 21 RW4-Fz
       % 22 RW4-Mx
       % 23 RW4-My
       % 24 RW4-Mz
       % 25 CYO-Fx
       % 26 CYO-Fy
       % 27 CYO-Fz
       % 28 ACS-Mx
       % 29 ACS-My
       % 30 ACS-Mz
str_in= str2mat('1 RW1-Fx','2 RW1-Fy','3 RW1-Fz','4 RW1-Mx','5 RW1-My',...
   '6 RW1-Mz','7 RW2-Fx','8 RW2-Fy','9 RW2-Fz','10 RW2-Mx','11 RW2-My',...
   '12 RW2-Mz','13 RW3-Fx','14 RW3-Fy','15 RW3-Fz','16 RW3-Mx','17 RW3-My',...
   '18 RW3-Mz','19 RW4-Fx','20 RW4-Fy','21 RW4-Fz','22 RW4-Mx','23 RW4-My',...
   '24 RW4-Mz','25 CYO-Fx','26 CYO-Fy','27 CYO-Fz','28 ACS-Mx','29 ACS-My',...
   '30 ACS-Mz');

      
 %  define displacement output dofs
 %  output nodes: 
 %   129,149,169  PM petal vertex nodes
 %   202  SM vertex
 %   207  Detector (node where back end optics and detector are located)
 %   84           ST Star Tracker (S/C node) angles

 dg=[...
    129  1 1 1 1 1 1
    149  1 1 1 1 1 1
    169  1 1 1 1 1 1
    202  1 1 1 1 1 1 
    207  1 1 1 1 1 1
     84  0 0 0 1 1 1 
 ];
 %  define velocity output dofs (three angular velocities at S/C bus node)
 vg = [ 84 0 0 0 1 1 1 ];
       % OUTPUT VECTOR
       %  1 PMS1-x
       %  2 PMS1-y
       %  3 PMS1-z
       %  4 PMS1-rx
       %  5 PMS1-ry
       %  6 PMS1-rz
       %  7 PMS2-x
       %  8 PMS2-y
       %  9 PMS2-z
       % 10 PMS2-rx
       % 11 PMS2-ry
       % 12 PMS2-rz
       % 13 PMS3-x
       % 14 PMS3-y
       % 15 PMS3-z
       % 16 PMS3-rx
       % 17 PMS3-ry
       % 18 PMS3-rz
       % 19 SM-x
       % 20 SM-y
       % 21 SM-z
       % 22 SM-rx
       % 23 SM-ry
       % 24 SM-rz
       % 25 IM-x
       % 26 IM-y
       % 27 IM-z
       % 28 IM-rx
       % 29 IM-ry
       % 30 IM-rz
       % 31 SC-rx
       % 32 SC-ry
       % 33 SC-rz
       % 34 SC-rrx rate
       % 35 SC-rry rate
       % 36 SC-rrz rate
    str_out=str2mat('1 PMS1-x','2 PMS1-y','3 PMS1-z','4 PMS1-rx','5 PMS1-ry',...
   '6 PMS1-rz','7 PMS2-x','8 PMS2-y','9 PMS2-z','10 PMS2-rx','11 PMS2-ry',...
   '12 PMS2-rz','13 PMS3-x','14 PMS3-y','15 PMS3-z','16 PMS3-rx','17 PMS3-ry',...
   '18 PMS3-rz','19 SM-x','20 SM-y','21 SM-z','22 SM-rx','23 SM-ry','24 SM-rz',...
   '25 IM-x','26 IM-y','27 IM-z','28 IM-rx','29 IM-ry','30 IM-rz','31 SC-rx',...
   '32 SC-ry','33 SC-rz','34 SC-rrx','35 SC-rry','36 SC-rrz');

%za=0.005; % Enter global modal damping coefficient (refine later)
% assume baseline zeta=0.001 for elastic structural modes
damping=ones(1,length(omeg)-3);
damping([4:8 16:18])=5*ones(1,8); % sunshield damping 5%
damping([9:11 14:15 19:25])=20*ones(1,12); % isolator damping 10%
damping([12:13 26:27])=20*ones(1,4); % solar panel damping 5%
Damping=diag(damping);
za=zeta*damping;
% remove translational rigid body modes since unobservable/uncontrollable
phi=phi(:,4:nm+3);
omeg=omeg(4:nm+3);
nrbm=3;
%  compute state space model
[Ap,Bp,Cp,Dp,lb,lc] = mode2ss2(xyz,bc,nm,ig,dg,vg,nrbm,za,phi,omeg);
sysp=ss(Ap,Bp,Cp,Dp);
O=diag(omeg);

if verification
%  (2) build state space model manually to doublecheck results
%  inputs
n_in=30;
dof_in=[ bci(78,:) bci(79,:) bci(80,:) bci(81,:)  bci(202,[1 2 3])  bci(82,4:6)]';
betau=zeros(ndof,n_in);
for ind=1:n_in
   betau(dof_in(ind),ind)=1;
end
n_out=36;
dof_out=[ bci(128,:) bci(148,:) bci(168,:) bci(197,:) bci(202,:) bci(83,4:6)]';
betay=zeros(length(dof_out),ndof);
for ind=1:length(dof_out)
   betay(ind,dof_out(ind))=1;
end
dof_outrate=[bci(83,4:6)]';
betaydot=zeros(length(dof_outrate),ndof);
for ind=1:length(dof_outrate)
   betaydot(ind,dof_outrate(ind))=1;
end
% assemble matrices
O=diag(omeg);
Z=zeta*Damping;
Apm=[zeros(nm,nm) eye(nm); -O.^2 -2*Z.*O];
Bpm=[zeros(nm,n_in); phi'*betau];
Cpm=[betay*phi zeros(length(dof_out),nm); 
   zeros(length(dof_outrate),nm) betaydot*phi];
Dpm=zeros(n_out,n_in);
syspm=ss(Apm,Bpm,Cpm,Dpm);
if diagnostics
   disp(['Finished Manual SS assembly'])
   end
end
% conclusion: both manual method and mode2ss2 give identical results ... GOOD !
if 0
nexus2xview
if diagnostics
   disp(['Created NEXUS file in XVIEW format'])
end
end
