% nexus_acsloop.m 
% NEXUS Attitude Control Model -  Version 1.0
% 
% Purpose
% ===========================================================
% This function closes a simple ACS controller
% around the NEXUS model in
% order to stabilize the rotational rigid body modes.  The ACS 
% controller consists of parallel PID-controller loops with 
% 2nd-order rolloff LPF. The MIMO controller is then rotated 
% from the decoupled loops in principal axes of inertia 
% to the spacecraft axes.
%
% Inputs
% ===========================================================
% mass: Rigid body mass matrix (6x6)
% fc: Target ACS bandwidth[Hz]
% plotflag: set to 1 to compare closed loop and open loop FRF's
% (Caution: setting plotflag to 1 for the full order model is slow)
% Author: Olivier L. de Weck
% Copyright: MIT Space Systems Laboratory
% Date: 25 May 2000
% Modified for use on NEXUS: 9 June 2001
%=================================
% 1. Define ACS Control parameters
%=================================
% fc=0.01;                % ACS control bandwidth [Hz]
Kcc=Kc*ones(1,3);         % FSM-to-ACS cross coupling gain
wc=2*pi*fca;
zetac=sqrt(2)/2;   % Butterworth pattern - critical damping
% Extract Inertia's from  6x6 rigid body mass matrix
% Find principal axes of inertia 
load mass
[P,I]=eig(mass(4:6,4:6));
I=diag(I);
% Make sure the principle axes are a right hand rule coordinate system
if ((cross(P(:,1),P(:,2))'*P(:,3)) <0)  %P1 x P2 = -P3
   P(:,1)=-P(:,1);
end
%====================================
%  2. Design ACS controller for NEXUS
%====================================
for ind=1:3
   Kp(ind)=I(ind)*wc^2;        % proportional gain
   Ki(ind)=0.05;               % integral gain fixed
   Kd(ind)=2*zetac*wc*I(ind);  % derivative gain
   wl=7*wc;                    % corner frequency for LPF
   zetal=sqrt(2)/2;            % LPF critical damping
   
   %numPD=[Kd(ind) Kp(ind) Ki(ind)];
   %denPD=[1 0];
   %numLPF=[wl^2];
   %denLPF=[1 2*zetal*wl wl^2];
   %numACS=conv(numPD,numLPF);
   %denACS=conv(denPD,denLPF);
   %[Atemp,Btemp,Ctemp,Dtemp]=tf2ss(numACS,denACS);
   Atemp=[-2*zetal*wl 1 0; -wl^2 0 1;  0 0 0];
   %Btemp=[Kd(ind)*wl^2 Kp(ind)*wl^2 Ki(ind)*wl^2]';
   Btemp=[0 0 0; Kp(ind)*wl^2 Kd(ind)*wl^2 Kcc(ind)*wl^2; Ki(ind)*wl^2 0 0];
   Ctemp=[1 0 0]; Dtemp=[0 0 0];
   sysTEMP=ss(Atemp,Btemp,Ctemp,Dtemp);
   if ind==1
      sysACS=sysTEMP;
   else
      sysACS=parallel(sysACS,sysTEMP,[],[],[],[]);
   end
end  
[Aca,Bca,Cca,Dca]=ssdata(sysACS);
% Transform to non-principle coordinates
Tsel=[0 0 0 1 0 0 0 0 0;
      1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0
      0 0 0 0 1 0 0 0 0;
      0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 1 0;
      0 0 0 0 0 1 0 0 0;
      0 0 1 0 0 0 0 0 0
      0 0 0 0 0 0 0 0 1];

Bca=Bca*Tsel*[P' zeros(3,5); zeros(3,3) P' zeros(3,2); zeros(3,6) P'*H];
Cca=(-1)*P*Cca;
Dca=(-1)*P*Dca*Tsel*[P' zeros(3,5); zeros(3,3) P' zeros(3,2); zeros(3,6) P'*H];
%========================================
   

