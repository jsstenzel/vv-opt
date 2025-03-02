% nexus_optics.m
% Build optical linear sensitivity matrices
% from raw MACOS output (prescription: xci.in)
% dWo, 7 June 2001

mm2m=0.001;   % conversion from mm to m
m2micron=1e6; % conversion from m to micrometers (microns)
m2nm=1e9;     % conversion from m to nanometer
R=[0 1 0; 0 0 1; 1 0 0];  % rotation matrix from MACOS frame to FEM frame
% Note in FEM frame the +Y axis is optical boresight, in MACOS it is +Z

load nexopt   % previously computed with MACOS
   
% 1. source rotation -> plate scale
Kr=(R*dcdu(:,99:101)')';
H=Kr'*inv(Kr*Kr')*mm2m*dcdu(:,97:98);
sourcerot=mm2m*dcdu(:,99:101);
psc=mean(abs([sourcerot(1,1) sourcerot(2,2)])); %plate scale meters/rad
% psc corrsponds to 0.4894 microns/milli-arcsecond

%2. FSM actuation matrix
Kfsm=mm2m*dcdu(:,97:98);  % meters centroid pos XY change/rad FSM tilt

% 3. centroid sensitivities
% reorganize by optical element
% rotate from MACOS fram to FEM frame
% put translations X Y Z first, then rotations RX RZ RZ
dcdupmseg1=m2micron*[(R*dcdu(:,4:6)')' mm2m*(R*dcdu(:,1:3)')'];
dcdupmseg2=m2micron*[(R*dcdu(:,16:18)')' mm2m*(R*dcdu(:,13:15)')']; 
% note PMSeg 2 and 3 are switched in FEM
dcdupmseg3=m2micron*[(R*dcdu(:,10:12)')' mm2m*(R*dcdu(:,7:9)')'];
dcdusm=m2micron*[(R*dcdu(:,22:24)')' mm2m*(R*dcdu(:,19:21)')'];
% add all sensitivities for the back optics together
% Fold, TM, DM, FSM etc.... since they are moving together
% as a rigid body in the FEM, recall RBE2 from grid 84 (S/C) to 207 (IM)
dcduimrot=[sum(dcdu(:,25:6:91)')' sum(dcdu(:,26:6:92)')' sum(dcdu(:,27:6:93)')'];
dcduimtrans=[sum(dcdu(:,28:6:94)')' sum(dcdu(:,29:6:95)')' sum(dcdu(:,30:6:96)')'];
dcduim=m2micron*[(R*dcduimtrans')' mm2m*(R*dcduimrot')'];
% overwrite dcdu
dcdu=inv(m2micron)*[dcdupmseg1 dcdupmseg2 dcdupmseg3 dcdusm dcduim]; 
%back to m due to numerical conditioning

% wavefront sensitivities
% subsample only every 10-th ray since we only care about RMS WFE
% dwdu=dwdu(1:10:end,:);
% reorganize by optical element
% rotate from MACOS fram to FEM frame
% put translations X Y Z first, then rotations RX RZ RZ
dwdupmseg1=m2nm*[(R*dwdu(:,4:6)')' mm2m*(R*dwdu(:,1:3)')'];
dwdupmseg2=m2nm*[(R*dwdu(:,16:18)')' mm2m*(R*dwdu(:,13:15)')']; 
% note PMSeg 2 and 3 are switched in FEM
dwdupmseg3=m2nm*[(R*dwdu(:,10:12)')' mm2m*(R*dwdu(:,7:9)')'];
dwdusm=m2nm*[(R*dwdu(:,22:24)')' mm2m*(R*dwdu(:,19:21)')'];
% add all sensitivities for the back optics together
% Fold, TM, DM, FSM etc.... since they are moving together
% as a rigid body in the FEM, recall RBE2 from grid 84 (S/C) to 207 (IM)
dwduimrot=[sum(dwdu(:,25:6:91)')' sum(dwdu(:,26:6:92)')' sum(dwdu(:,27:6:93)')'];
dwduimtrans=[sum(dwdu(:,28:6:94)')' sum(dwdu(:,29:6:95)')' sum(dwdu(:,30:6:96)')'];
dwduim=m2nm*[(R*dwduimtrans')' mm2m*(R*dwduimrot')'];
dwdufsm=m2nm*mm2m*dwdu(:,97:98);
% overwrite dcdu
dwdu=[dwdupmseg1 dwdupmseg2 dwdupmseg3 dwdusm dwduim dwdufsm];
dwdu=dwdu(1:10:end,:); % sample every 10th ray
dwdu=inv(m2nm)*dwdu; % back to m

% save nexus_optics psc Kfsm dcdu dwdu 



