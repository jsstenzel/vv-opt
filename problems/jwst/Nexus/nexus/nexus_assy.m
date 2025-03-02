% nexus_assy.m

% conversion factors
r2a = 3600*180/pi;	% radians to arc-sec
a2r = 1/r2a;			% arc-sec to radians
nray= 134;
%tic
nexus_rwa

nexus_cryo

nexus_acsnoise

nexus_gs
%disp(['Assy' num2str(toc)])
%tic
nexus_plant
%disp(['plant' num2str(toc)])
%tic
nexus_optics

nexus_acsloop

nexus_fsm
%disp(['optcont' num2str(toc)])

%-------------------------------------------------------------
%RWA disturbance matrix
Rrwa=[R1 zeros(6,18); zeros(6,6) R2 zeros(6,12);
      zeros(6,12) R3 zeros(6,6); zeros(6,18) R4];

ndw=size(Adw,1);
Adw4=[Adw zeros(ndw,3*ndw); zeros(ndw,ndw) Adw zeros(ndw,2*ndw);
   zeros(ndw,2*ndw) Adw zeros(ndw,ndw); zeros(ndw,3*ndw) Adw];
Bdw4=[Bdw zeros(ndw,3); zeros(ndw,1) Bdw zeros(ndw,2);
   zeros(ndw,2) Bdw zeros(ndw,1); zeros(ndw,3) Bdw];
Cdw4=Rrwa*[Cdw zeros(6,3*ndw); zeros(6,ndw) Cdw zeros(6,2*ndw);
   zeros(6,2*ndw) Cdw zeros(6,ndw); zeros(6,3*ndw) Cdw];
Ddw4=Rrwa*zeros(24,4);
%-------------------------------------------------------------
%   Assemble overall ss system
%   d1-10  10 white noise sources
%   d1-4    RWA
%   d5      CRYO
%   d6-8    ACS
%   d9-10   GS
% state count
ndw=size(Adw4,1);
ndc=size(Adc,1);
nds=size(Ads,1);
ndg=size(Adg,1);
np=size(Ap,1);
nca=size(Aca,1);
ncf=size(Acf,1);

% partial matrices
Bp1=Bp(:,1:24); Bp2=Bp(:,25:27); Bp3=Bp(:,28:30);
Cp1=Cp(1:30,:); Cp2=Cp(31:33,:); Cp3=Cp(34:36,:);
Bca1=Bca(:,1:3); Bca2=Bca(:,4:6); Bca3=Bca(:,7:8);
dwdu1=dwdu(:,1:30); dwdu2=dwdu(:,31:32);

% closed loop integrated model
Azd=[Adw4 zeros(ndw,ndc+nds+ndg+np+nca+ncf);
   zeros(ndc,ndw) Adc zeros(ndc,nds+ndg+np+nca+ncf);
   zeros(nds,ndw+ndc) Ads zeros(nds,ndg+np+nca+ncf);
   zeros(ndg,ndw+ndc+nds) Adg zeros(ndg,np+nca+ncf);
   Bp1*Cdw4  Bp2*Cdc zeros(np,nds+ndg) Ap  Bp3*Cca zeros(np,ncf);
   zeros(nca,ndw+ndc) Bca2*Cds zeros(nca,ndg) Bca1*Cp3+Bca2*Cp2 Aca Bca3*(-1)*Ccf;
   zeros(ncf,ndw+ndc+nds) Bcf*psc*Cdg Bcf*dcdu*Cp1 zeros(ncf,nca) Acf-Bcf*Kfsm*(-1)*Ccf];

Bzd=[1*Bdw4 zeros(ndw,6);
   zeros(ndc,4) 1*Bdc zeros(ndc,5);
   zeros(nds,5) 1*Bds zeros(nds,2);
   zeros(ndg,8) 1*Bdg;
   zeros(np+nca+ncf,10)];

Czd=[zeros(size(dwdu,1),ndw+ndc+nds+ndg) m2nm*dwdu1*Cp1 zeros(size(dwdu,1),nca) m2nm*dwdu2*(-1)*Ccf;
   zeros(size(dcdu,1),ndw+ndc+nds+ndg) m2micron*dcdu*Cp1 zeros(size(dcdu,1),nca) -m2micron*Kfsm*(-1)*Ccf];

Dzd=zeros(size(dwdu,1)+size(dcdu,1),10);
