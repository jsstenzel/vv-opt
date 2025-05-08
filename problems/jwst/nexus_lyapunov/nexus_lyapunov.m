%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RMMS_WFE_lyap, RSS_LOS_lyap] = nexus_lyapunov(x)
tstart=cputime;

%%%Unwrap the input struct x
Ru=x.Ru;      % Ru:  upper operational wheel speed  [RPM]
Us=x.Us;       % Us:  static wheel imbalance         [gcm]    (%0.7160 test)
Ud=x.Ud;        % Ud:  dynamic wheel imbalance        [gcm^2]  (%29.536 test)
fc=x.fc;        % fc:  cryocooler drive frequency     [Hz]
Qc=x.Qc;     % Qc:  cryocooler attenuation factor  [-] 
Tst=x.Tst;	     % Tst: Star Tracker update rate       [sec]
Srg=x.Srg;    % Srg: rate gyro noise intensity      [(rad/s)^2/Hz] 
Sst=x.Sst;        % Sst: Star Tracker one sigma         [arcsec]
Tgs=x.Tgs;     % Tgs: Guider integration time         [sec]
m_SM=x.m_SM;          % m_SM: mass of secondary mirror (SM) [kg]
m_SMhub=x.m_SMhub;  %SM hub mass
I_SMhubt=x.I_SMhubt; %SM hub inertia transverse
I_SMhuba=x.I_SMhuba; %SM hub inertia axial
m_RW=x.m_RW; % RW mass (ITHACO E-wheel)
K_yPM=x.K_yPM; % K_yPM: PM bipod act trans stiffness [N/m]
I_xRWA=x.I_xRWA;  % Ix moment of inertia of RWA chassis
I_yRWA=x.I_yRWA;  % Iy moment of inertia of RWA chassis
I_zRWA=x.I_zRWA; % Iz moment of inertia of RWA chassis
I_RWt=x.I_RWt; % Transverse Inertia of ITHACO E reaction wheel
I_RWa=x.I_RWa; % Axial Inertia of ITHACO E reaction wheel
m_ISO=x.m_ISO;  % mass of RWA isolator strut actuator
I_ISOa=x.I_ISOa;  % RWA isolator strut axial inertia
I_ISOt=x.I_ISOt;  % RWA isolator strut transverse inertia
K_yISO=x.K_yISO;  % Hexapod (i.e. secondary mirror support structure) isolator strut axial stiffness
K_xISO=x.K_xISO;  % Hexapod (i.e. secondary mirror support structure) isolator strut transverse stiffness
m_RWAchx=x.m_RWAchx; % mass RWA pyramid chassis
I_bus=x.I_bus;    % NEXUS spacecraft bus principal inertia
m_bus=x.m_bus;    % m_bus:  NEXUS spacecraft bus mass   [kg]
m_prop=x.m_prop;   % propulsion system mass (incl. propellant)
I_propt=x.I_propt;  % propulsion system transverse inertia % I_propt: propulsion sys trv inertia [kgm^2]
I_propa=x.I_propa;  % propulsion system axial inertia
m_instr=x.m_instr;  % instrument mass (incl. CCD detector)
I_i1=x.I_i1;  % Instrument moment of inertia 1-axis
I_i2=x.I_i2;  % Instrument moment of inertia 2-axis
I_i3=x.I_i3;  % Instrument moment of inertia 3-axis
% SM supporting structure (spider): Tangential Bipod design
A_sptop=x.A_sptop;  % Cross sectional area of SM spider (top)
D_sp=x.D_sp;           % SM spider (bottom) diameter
t_sp=x.t_sp;         % SM spider (bottom) wall thickness  [m]
A_spbot=x.A_spbot; % Cross sectional area of SM spider (bottom)
J_sp=x.J_sp;   % SM Spider torsional moment of inertia
I_sp=x.I_sp;   % SM Spider bending moment of inertia
I_ss=x.I_ss;     % Sunshield out-of-plane bending moment of inertia, % I_ss: Sunshield bending mof inertia [m^4]     
K_rad1=x.K_rad1;   % rotational stiffness of SM and PM actuators in 1 direction
K_rad2=x.K_rad2;   % rotational stiffness of SM and PM actuators in 2 direction 
K_rISO=x.K_rISO;  % K_rISO: rot stiff RWA iso struts [Nm/rad]
K_aISO=x.K_aISO; % axial rot RWA iso struts 
K_act1=x.K_act1;  % actuator axial stiffness 
K_act2=x.K_act2;  % actuator radial stiffness +12
I_iso=x.I_iso;    % RW wheel isolator bending m.o.I. [m^4]
K_zpet=x.K_zpet;  % K_zpet: Deployable petal hinge stiffness [N/m] (nom:0.90000E+08)
zeta=x.zeta;         % zeta: Global modal damping ratio    [-]
lambda=x.lambda;	     % lambda: center optical wavelength   [m]
Ro=x.Ro;		        % Ro:  optical surface transmissivity [-]
QE=x.QE ;		     % QE:  CCD quantum efficiency    [-]
Mgs=x.Mgs;             % Mgs: magnitude of faint guide star  [mag]
fca=x.fca;           % fca: Target ACS control bandwidth    [Hz]
Kc=x.Kc;             % Kc:  FSM-to-ACS cross coupling gain  [0-1]
Kcf=x.Kcf;    % Kcf: FSM controller gain             [-]
nray=x.nray;
caxi=x.caxi;
crad=x.crad;
ctor=x.ctor;
haxi=x.haxi;
hrad=x.hrad;
htor=x.htor;
zeta1=x.zeta1;     % empirical number
a=x.a;         % HPF corner fraction of top wheel speed
wheel_locs=x.wheel_locs;
n=x.n;       % number of harmonics to include (including fundamental)
h=x.h;
C=x.C;
%ACS_GYRO_nw=x.ACS_GYRO_nw;	% gyro rate random walk (rad^2/s^3) or (rad^2/s^4)/Hz
%ACS_GYRO_nf=x.ACS_GYRO_nf;	   % gyro rate flicker noise (rad^2/s^2) 
Nsurf=x.Nsurf;	         			% number of optical surfaces before detector
D=x.D;		            		% aperture diameter (m) NEXUS
BP=x.BP;		         		% bandpass (microns)
PH=x.PH;		      		% photons/m^2/micron/sec
R0=x.R0;						   % detector readout noise, 60 electrons per pixel, 4 pixels
mass=x.mass;
FgsNom=x.FgsNom;			      % nominal FSM sample rate (Hz)
% STUFF I ADDED
K_pm1=x.K_pm1;
K_pm3=x.K_pm3;
K_pm4=x.K_pm4;
K_pm5=x.K_pm5;
K_pm6=x.K_pm6;
K_act_pm2=x.K_act_pm2;
K_act_pm3=x.K_act_pm3;
K_act_pm4=x.K_act_pm4;
K_act_pm5=x.K_act_pm5;
K_act_pm6=x.K_act_pm6;
K_xpet=x.K_xpet;
K_cryo=x.K_cryo;
K_IA=x.K_IA;
c_cryo=x.c_cryo;
c_IA=x.c_IA;
c_RWA=x.c_RWA;
c_RWAI=x.c_RWAI;
c_SM_act=x.c_SM_act;
c_PM=x.c_PM;
c_PM_act=x.c_PM_act;
c_petal=x.c_petal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Set up flags and conversions
diagnostics=1;       % print diagnostic messages on command line
plotflag=0;          % plots generated in subroutines
r2a = 3600*180/pi;	% radians to arc-sec
a2r = 1/r2a;			% arc-sec to radians


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_rwa: Reaction wheels
%----------------------------------------------------------
dR=Ru/2;              % +/- wheel speed deviation
scale(1)=Us/0.7160;  % scale factor for Us %TODO
scale(2)=Ud/29.536;  % scale factor for Ud %TODO
f=logspace(-2,3,500);
type=0;             % uniform wheel speed pdf
model=[];

Rl=0;
Crad1=((2*pi)^2/(1000*60^2*100))*Us;
Caxi1=caxi(1);  %note axial force not affected by imbalances
Ctor1=((2*pi)^2/(1000*60^2*100^2))*Ud;
Sradm=pi*Crad1^2/(2*(Ru-Rl)*(pi/30)^5)*(2*pi*(Ru/60))^4; %PSD value
Saxim=pi*Caxi1^2/(2*(Ru-Rl)*(pi/30)^5)*(2*pi*(Ru/60))^4; %PSD value
Storm=pi*Ctor1^2/(2*(Ru-Rl)*(pi/30)^5)*(2*pi*(Ru/60))^4; %PSD value
f1max=Ru/60;
f2r=max(hrad)*f1max;
f2a=max(haxi)*f1max;
f2t=max(htor)*f1max;
krad1=4*pi*zeta1*f1max*sqrt(Sradm);
kaxi1=4*pi*zeta1*f1max*sqrt(Saxim);
ktor1=4*pi*zeta1*f1max*sqrt(Storm);

wl=f1max*2*pi;
wh=2*pi*[f2r f2a f2t];
kr=[krad1 kaxi1 ktor1];
% radial force ss approximation
Adr_r=[0 1 0 0; -a*wl*wh(1) -(a*wl+wh(1)) 0 0;
   0 0 0 1; 0 wh(1) -wl^2 -2*zeta1*wl];
Bdr_r=[0 1 0 0]';
Cdr_r=[0 0 0 kr(1)];
Ddr_r=[0]; sysr=ss(Adr_r,Bdr_r,Cdr_r,Ddr_r);
% axial force ss approximation
Adr_a=[0 1 0 0; -a*wl*wh(2) -(a*wl+wh(2)) 0 0;
   0 0 0 1; 0 wh(2) -wl^2 -2*zeta1*wl];
Bdr_a=[0 1 0 0]';
Cdr_a=[0 0 0 kr(2)];
Ddr_a=[0]; sysa=ss(Adr_a,Bdr_a,Cdr_a,Ddr_a);
% radial torque ss approximation
Adr_t=[0 1 0 0; -a*wl*wh(3) -(a*wl+wh(3)) 0 0;
   0 0 0 1; 0 wh(3) -wl^2 -2*zeta1*wl];
Bdr_t=[0 1 0 0]';
Cdr_t=[0 0 0 kr(3)];
Ddr_t=[0]; syst=ss(Adr_t,Bdr_t,Cdr_t,Ddr_t);
% place ss in parallel and duplicate radial force and radial torque channels
Adw=[Adr_r,zeros(4,8); zeros(4,4) Adr_a zeros(4,4); zeros(4,8) Adr_t];
Bdw=[Bdr_r ; Bdr_a ; Bdr_t];
Cdw=[Cdr_r zeros(1,8); Cdr_r zeros(1,8); zeros(1,4) Cdr_a zeros(1,4);
   zeros(1,8) Cdr_t; zeros(1,8) Cdr_t; zeros(1,12)];
Ddw=[zeros(6,1)];

for ind=1:4
   temp=wheel_locs(ind,2:4)-wheel_locs(5,2:4);
   r(:,ind)=temp'/norm(temp,2); %normalized wheel axis orientation
   beta(ind)=0;
   gamm(ind)=asin(-r(2,ind));
   theta(ind)=acos(r(3,ind)/cos(gamm(ind)));
end
% compute rotation matrices Ri , i=1,2,3,4
for ind =1:4
  cb = cos(beta(ind)) ;
  sb = sin(beta(ind)) ;
  ct = cos(theta(ind)) ;
  st = sin(theta(ind)) ;
  cg = cos(gamm(ind)) ;
  sg = sin(gamm(ind)) ;
  
  R = [cb*ct  cb*st*sg-sb*cg  cb*st*cg+sb*sg ;
       sb*ct  sb*st*sg+cb*cg  sb*st*cg-cb*sg ;
       -st        ct*sg           ct*cg     ] ;        % 3x3 rotation matrix
    
    eval(['R' num2str(ind) '=[R zeros(3,3); zeros(3,3) R];'])
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_cryo: Cryo cooler noise
wc=2*pi*fc; % drive frequency
Q=Qc*ones(1,n);  
%h just goes 1 to 6, not sure what it represents
%C has two rows of values in [N], apparently from Cassini -- these are magnitudes of the disturbance?

 for ind1=1:2
   for ind2=1:n
         zij(ind1,ind2)=inv(2*wc*h(ind2));
         kij(ind1,ind2)=2*zij(ind1,ind2)*wc*h(ind2)*Q(ind2)*C(ind1,ind2);
   end
 end
Adc=[0 1 zeros(1,10); -wc^2 -2*zij(1,1)*wc zeros(1,10);
   zeros(1,2) 0 1 zeros(1,8); zeros(1,2) -2^2*wc^2 -2*zij(1,2)*2*wc zeros(1,8);
   zeros(1,4) 0 1 zeros(1,6); zeros(1,4) -3^2*wc^2 -2*zij(1,3)*3*wc zeros(1,6);
   zeros(1,6) 0 1 zeros(1,4); zeros(1,6) -1^2*wc^2 -2*zij(2,1)*1*wc zeros(1,4);
   zeros(1,8) 0 1 zeros(1,2); zeros(1,8) -2^2*wc^2 -2*zij(2,2)*2*wc zeros(1,2);
   zeros(1,10) 0 1 ; zeros(1,10) -3^2*wc^2 -2*zij(2,3)*3*wc ];
Bdc=[0 1 0 1 0 1 0 1 0 1 0 1]';
Cdc=[...
      zeros(1,6) 0 kij(2,1) 0 kij(2,2) 0 kij(2,3) ; 
      zeros(1,6) 0 kij(2,1) 0 kij(2,2) 0 kij(2,3) ;
      0       kij(1,1) 0 kij(1,2) 0 kij(1,3) zeros(1,6)];
Ddc=[zeros(3,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_acsnoise: ACS noise
ACS_ST_period=Tst; % Star Tracker update rate (sec) nom: Cassini 100 sec
ACS_GYRO_nr=Srg;       	% gyro rate noise (rad^2/s) or (rad/s)^2/Hz nom: 3e-14ST
ACS_ST_sigma=Sst;             % Star Tracker one sigma resolution (arcsec)
Sstrad=Sst/60/60*pi/180;	%convert to radians
ACS_ST_sigma=ACS_ST_sigma/60/60*pi/180;	%convert to radians

%  frequency base
FreqBase=logspace(-5,3,1000);
rad2arcsec=(3600*180/pi);
f1=find(FreqBase <= 2/ACS_ST_period);

Kn=sqrt((1/2)*(Sstrad^2+Srg./(2/Tst).^2));  % ACS noise gain

% state space filter approximation to composite
% ACS sensor noise model.
wn=2*pi*(2/Tst);
Ads=[-wn]; 
Bds=[1]';
Cds=Kn*[wn];
Dds=[0]; 
SYSds=ss(Ads,Bds,Cds,Dds); %send to Bode plot, see magnitude

[magds,phsds]=bode(SYSds,FreqBase*2*pi); %magnitude plot & phase plot
magds=abs(squeeze(magds));
Snn=magds.^2; %this gets us rms
%TODO what does intrms do
%rms is a number under a frequency range; intrms integrates that in frequency bands
%bode plots how rms error changes over frequency space
rmsssrad=intrms(2*Snn,FreqBase);
rmsssarcsec=rad2arcsec*rmsssrad; %this is likely the control error
if diagnostics
  disp(['RMS Approximation (arcsec) = ',num2str(rmsssarcsec)]);
  end
  if plotflag
	 hold on
	 loglog(FreqBase,rad2arcsec*(sqrt(Snn)),'c-')
	 h=gcf; line_handle=findobj(h,'Type','line'); set(line_handle,'LineWidth',2);
	 legend('Star Tracker',...
	 'Gyro Rate Random Walk',...
	 'Gyro Rate Noise',...
	 'Composite',...
	 'State Space Approx',3);
	 text(1e-5,1e-3,['RMS Model [asec] = ',num2str(ACS_PSD_noise_rms_arcsec)])
	 text(1e-5,0.2e-3,['RMS Approximation [asec] = ',num2str(rmsssarcsec)])
  end

%---------------------------------------
% State space for 3 axes
%---------------------------------------
Systmp=parallel(SYSds,SYSds,[],[],[],[]);
SYSds=parallel(Systmp,SYSds,[],[],[],[]);
[Ads,Bds,Cds,Dds]=ssdata(SYSds);
% Note state space outputs of SYSds are in rad !

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_gs: Optics

% overall transmissibility 
Reff = Ro^Nsurf;
% area of primary mirror
Area = pi * D^2 / 4;
% slope of centroiding transfer function
ks = (16/3/pi) * (D/lambda) / r2a;
alpha=Reff * QE * Area * BP * 10.^(-0.4.*Mgs) * PH;
% total number of detected photons
N = alpha * Tgs;
% raw NEA
NEA = sqrt(1+R0./N)/ks./sqrt(N); % compare this with rms output from SS

if diagnostics
   disp('Compute ss model of FGS noise')
end
zeta_fgs=0.707;
w_fgs=2*pi/Tgs;
k_NEA=2*NEA*sqrt(zeta_fgs/w_fgs);
Adg=[ 0 1 ; -w_fgs^2 -2*zeta_fgs*w_fgs ]; 
Adg=[Adg zeros(2,2); zeros(2,2) Adg];   % 2 channels in parallel
Bdg=[0 1]';Bdg=[Bdg [0 0]';[0 0]' Bdg];
Cdg=[ k_NEA*w_fgs^2 0]; 
Cdg=a2r*[Cdg [0 0]; [0 0] Cdg]; % output in radians NEA
Ddg=a2r*zeros(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_plant: Plant & FEM
if diagnostics
   tic
end

% --------------------------------------------------------------
% -------------------------   plant stuff  ---------------------
% --------------------------------------------------------------
% ======================================
% 1. Create grid location matrix
% ======================================
%assigns numbers to points on an xyz grid
%TODO categorize everything at high level: Optical Telescope element (OTE), Integrated Science Instrument Model (ISIM), Spacecraft Element (SCE)
%TODO identify or modify to specifically include: 
	%IEC (contains electronics in ISIM)
	%IEC harness??
	%sunshield
	%spacecraft bus
	%magnetic tuned mass dampener (isolator for secondary mirror support structure, the little mirrors)
	%Cryocooler CJAA (on OTE/ISIM payload, purpose unclear to me)
	%Cryocooler CCA (isolates compressors from rest of SCE)
	%isolator assembly (connects SCE to DTA)
	%PT cryocooler compressor
	%JT cryocooler compressor
	%deployable tower assembly (connects SCE to OTE)
	%Reaction Wheel Isolator Assembly (RWIA)
%TODO note that many of these points are actually repeated - namely all the ones used for nicelas2!
%not sure what it means... but, correlate them in the notation here
xyz=[ ...
1  0.00000E+00  3.41300E+00  2.72000E-01 %spider
2  0.00000E+00  2.87800E+00  7.20000E-01
3 -2.35560E-01  3.41300E+00 -1.36000E-01
4 -6.23540E-01  2.87800E+00 -3.60000E-01 %spider
5  2.35560E-01  3.41300E+00 -1.36000E-01 %spider
6  6.23540E-01  2.87800E+00 -3.60000E-01 %spider
7 -6.50000E-01  0.00000E+00  1.10000E+00 %isolator?
8  6.50000E-01  0.00000E+00  1.10000E+00 %isolator?
9 -1.13800E+00  0.00000E+00 -6.30000E-02 %isolator?
10 -6.50000E-01  0.00000E+00 -9.07000E-01 %isolator?
11  6.50000E-01  0.00000E+00 -9.07000E-01 %isolator?
12  1.13800E+00  0.00000E+00 -6.30000E-02 %isolator?
13  0.00000E+00  3.41300E+00  0.00000E+00 %secondary mirror hub?
14 -1.89470E-01  3.45930E+00 -1.23830E-01 %secondary mirror actuator
15 -2.01970E-01  3.45930E+00 -1.02170E-01 %secondary mirror actuator
16  1.89470E-01  3.45930E+00 -1.23830E-01 %secondary mirror actuator
17  2.01970E-01  3.45930E+00 -1.02170E-01 %secondary mirror actuator
18 -1.25000E-02  3.45930E+00  2.26000E-01 %secondary mirror actuator
19  1.25000E-02  3.45930E+00  2.26000E-01 %secondary mirror actuator
20 -9.88200E-02  3.32485E+00  4.23890E-02 %secondary mirror actuator
21 -8.61200E-02  3.32485E+00  6.43860E-02 %secondary mirror actuator
22  0.00000E+00  3.29369E+00  0.00000E+00 %secondary mirror RBE2 %isolator?
23 -1.27000E-02  3.32485E+00 -1.06780E-01 %secondary mirror actuator
24  1.27000E-02  3.32485E+00 -1.06780E-01 %secondary mirror actuator
25  8.61200E-02  3.32485E+00  6.43860E-02 %secondary mirror actuator
26  9.88200E-02  3.32485E+00  4.23890E-02 %secondary mirror actuator
27  6.15960E-01 -3.08310E-01  1.72528E+00 %S/C top corner (4pts) %solar panel?
28 -6.15960E-01 -3.08310E-01  1.72528E+00 %S/C top corner (4pts)
29  0.00000E+00  1.24610E-01 -2.37940E-01 %primary mirror
30 -6.15960E-01 -1.46593E+00  1.30394E+00 %S/C top corner (4pts)
31 -9.88200E-02  3.32485E+00  4.23890E-02 %=20
32 -8.61200E-02  3.32485E+00  6.43860E-02 %=21
33  8.61200E-02  3.32485E+00  6.43860E-02 %=25
34  9.88200E-02  3.32485E+00  4.23890E-02 %=26
35  1.27000E-02  3.32485E+00 -1.06780E-01 %=24
36 -1.27000E-02  3.32485E+00 -1.06780E-01 %=23
37 -1.89470E-01  3.45930E+00 -1.23830E-01 %=14
38 -2.01970E-01  3.45930E+00 -1.02170E-01 %=15
39 -1.25000E-02  3.45930E+00  2.26000E-01 %=18
40  1.25000E-02  3.45930E+00  2.26000E-01 %=19
41  2.01970E-01  3.45930E+00 -1.02170E-01 %=17
42  1.89470E-01  3.45930E+00 -1.23830E-01 %=16
43  6.15960E-01 -8.87120E-01  1.51461E+00 %sunshield
44  0.00000E+00 -3.08310E-01  1.72528E+00 %sunshield
45 -6.15960E-01 -8.87120E-01  1.51461E+00 %sunshield
46  0.00000E+00 -1.46593E+00  1.30394E+00 %sunshield
47  0.00000E+00 -1.79290E+00  2.20228E+00 %mid point, S/C bottom edge (pitch)
48  6.15960E-01 -1.21409E+00  2.41295E+00 %solar array attach
49  0.00000E+00 -6.35280E-01  2.62362E+00 %mid point, S/C bottom edge (pitch)
50 -6.15960E-01 -1.21409E+00  2.41295E+00 %solar array attach
51  3.61596E+00  9.56000E-01  0.00000E+00
52  0.00000E+00  9.56000E-01 -3.61596E+00
53 -3.61596E+00  9.56000E-01  0.00000E+00
55  0.00000E+00  9.56000E-01  6.21596E+00
56 -3.81596E+00  0.00000E+00  0.00000E+00
57  3.81596E+00  0.00000E+00  0.00000E+00
58  0.00000E+00  0.00000E+00  0.00000E+00 %origin - not CM
59  0.00000E+00 -2.17070E+00  1.04742E+00
60  0.00000E+00 -2.87547E+00  7.90910E-01
61  0.00000E+00 -3.58024E+00  5.34390E-01
62  1.36596E+00 -8.87120E-01  1.51461E+00
63  2.11596E+00 -8.87120E-01  1.51461E+00
64  6.15960E-01 -1.46593E+00  1.30394E+00 %S/C top corner (4pts)
65  2.86596E+00 -8.87120E-01  1.51461E+00
66 -1.36596E+00 -8.87120E-01  1.51461E+00
67 -6.90000E-01  0.00000E+00 -1.50000E-01
68  6.90000E-01  0.00000E+00 -1.50000E-01
69  6.15960E-01 -6.35280E-01  2.62362E+00 %bottom 4 attach pts to launch vehicle
70 -6.15960E-01 -6.35280E-01  2.62362E+00 %bottom 4 attach pts to launch vehicle
71 -6.15960E-01 -1.79290E+00  2.20228E+00 %bottom 4 attach pts to launch vehicle
72  6.15960E-01 -1.79290E+00  2.20228E+00 %bottom 4 attach pts to launch vehicle
73 -6.54380E-01 -6.65870E-01  1.63769E+00
74 -4.04380E-01 -5.79050E-01  2.13009E+00
75 -6.54380E-01 -9.62070E-01  2.45148E+00
76 -4.04380E-01 -1.43192E+00  2.28047E+00
77 -6.54380E-01 -1.51874E+00  1.78807E+00
78 -4.04380E-01 -1.13572E+00  1.46668E+00
79 -4.74130E-01 -8.13750E-01  2.04467E+00 %RWA location
80 -4.74130E-01 -9.63310E-01  1.72393E+00 %RWA location
81 -4.74130E-01 -1.28405E+00  1.87349E+00 %RWA location
82 -4.74130E-01 -1.13448E+00  2.19423E+00 %RWA location
83 -5.67280E-01 -1.04890E+00  1.95908E+00 %center of RWA assembly (for torques)
84  0.00000E+00 -1.05061E+00  1.96378E+00 %spacecraft bus CM, star tracker?
85 -6.54380E-01 -6.65870E-01  1.63769E+00
86 -4.04380E-01 -5.79050E-01  2.13009E+00 %86 = 112 RWA attach point
87 -6.54380E-01 -9.62070E-01  2.45148E+00
88 -4.04380E-01 -1.43192E+00  2.28047E+00 %88 = 113 RWA attach point
89 -6.54380E-01 -1.51874E+00  1.78807E+00
90 -4.04380E-01 -1.13572E+00  1.46668E+00 %90 = 114 RWA attach point
91  0.00000E+00  0.00000E+00  0.00000E+00
92 -2.11596E+00 -8.87120E-01  1.51461E+00
93 -2.86596E+00 -8.87120E-01  1.51461E+00
94  0.00000E+00  1.00726E+00  2.20410E+00
95  0.00000E+00  2.32283E+00  2.68293E+00
96  0.00000E+00  3.63840E+00  3.16176E+00
97 -5.29380E-01 -9.00800E-01  1.55218E+00
98 -5.29380E-01 -6.22460E-01  1.88389E+00
99 -5.29380E-01 -1.32723E+00  1.62737E+00
100 -5.29380E-01 -1.47533E+00  2.03427E+00
101 -5.29380E-01 -7.70560E-01  2.29079E+00
102 -5.29380E-01 -1.19700E+00  2.36598E+00
103 -6.54380E-01 -6.65870E-01  1.63769E+00
104 -6.54380E-01 -6.65870E-01  1.63769E+00
105 -6.54380E-01 -9.62070E-01  2.45148E+00
106 -6.54380E-01 -9.62070E-01  2.45148E+00
107 -6.54380E-01 -1.51874E+00  1.78807E+00
108 -6.54380E-01 -1.51874E+00  1.78807E+00
109 -4.04380E-01 -5.79050E-01  2.13009E+00
110 -4.04380E-01 -1.43192E+00  2.28047E+00
111 -4.04380E-01 -1.13572E+00  1.46668E+00
112 -4.04380E-01 -5.79050E-01  2.13009E+00 %86 = 112 RWA attach point
113 -4.04380E-01 -1.43192E+00  2.28047E+00 %88 = 113 RWA attach point
114 -4.04380E-01 -1.13572E+00  1.46668E+00 %90 = 114 RWA attach point
115  4.34740E-01  1.54070E-01 -4.87220E-01
116  1.01090E-01  3.39207E+00 -1.15300E-01
117 -1.01090E-01  3.39207E+00 -1.15300E-01
118 -1.50400E-01  3.39207E+00 -2.99000E-02
119 -4.93100E-02  3.39208E+00  1.45190E-01
120  4.93100E-02  3.39207E+00  1.45190E-01
121  1.50400E-01  3.39207E+00 -2.99000E-02
122  4.34740E-01  2.07420E-01 -9.86380E-01
123  0.00000E+00  2.32800E-01 -1.23610E+00
124 -4.34740E-01  2.07420E-01 -9.86380E-01
125 -4.34740E-01  1.54070E-01 -4.87220E-01
126  0.00000E+00  1.51070E-01 -6.23780E-01
127  1.00030E-01  1.70140E-01 -7.95980E-01
128 -1.00030E-01  1.70140E-01 -7.95980E-01
129  0.00000E+00  2.24540E-01 -7.32210E-01 %deployable PM petal vertex
130  0.00000E+00  2.06200E-01 -7.28870E-01
131 -4.34740E-01  2.07420E-01 -9.86380E-01
132  0.00000E+00  2.32800E-01 -1.23610E+00
133  4.34740E-01  2.07420E-01 -9.86380E-01
134  4.34740E-01  1.54070E-01 -4.87220E-01
135  0.00000E+00  1.24610E-01 -2.37940E-01
136 -4.34740E-01  1.54070E-01 -4.87220E-01
137  0.00000E+00  1.51070E-01 -6.23780E-01
138 -1.00030E-01  1.70140E-01 -7.95980E-01
139  1.00030E-01  1.70140E-01 -7.95980E-01
140  2.06060E-01  1.24610E-01  1.18970E-01
141  2.04570E-01  1.54070E-01  6.20110E-01
142  6.36860E-01  2.07420E-01  8.69690E-01
143  1.07049E+00  2.32800E-01  6.18050E-01
144  1.07160E+00  2.07420E-01  1.16690E-01
145  6.39320E-01  1.54070E-01 -1.32890E-01
146  5.40210E-01  1.51070E-01  3.11890E-01
147  6.39330E-01  1.70140E-01  4.84610E-01
148  7.39350E-01  1.70140E-01  3.11360E-01
149  6.34110E-01  2.24540E-01  3.66100E-01 %fixed PM petals vertex
150  6.31220E-01  2.06200E-01  3.64430E-01
151  1.07160E+00  2.07420E-01  1.16690E-01
152  1.07049E+00  2.32800E-01  6.18050E-01
153  6.36860E-01  2.07420E-01  8.69690E-01
154  2.04570E-01  1.54070E-01  6.20110E-01
155  2.06060E-01  1.24610E-01  1.18970E-01
156  6.39320E-01  1.54070E-01 -1.32890E-01
157  5.40210E-01  1.51070E-01  3.11890E-01
158  7.39350E-01  1.70140E-01  3.11360E-01
159  6.39330E-01  1.70140E-01  4.84610E-01
160 -2.06060E-01  1.24610E-01  1.18970E-01
161 -6.39320E-01  1.54070E-01 -1.32890E-01
162 -1.07160E+00  2.07420E-01  1.16690E-01
163 -1.07049E+00  2.32800E-01  6.18050E-01
164 -6.36860E-01  2.07420E-01  8.69690E-01
165 -2.04570E-01  1.54070E-01  6.20110E-01
166 -5.40210E-01  1.51070E-01  3.11890E-01
167 -7.39350E-01  1.70140E-01  3.11360E-01
168 -6.39330E-01  1.70140E-01  4.84610E-01
169 -6.34110E-01  2.24540E-01  3.66110E-01 %fixed PM petals vertex
170 -6.31220E-01  2.06200E-01  3.64440E-01
171 -6.36860E-01  2.07420E-01  8.69690E-01
172 -1.07049E+00  2.32800E-01  6.18050E-01
173 -1.07160E+00  2.07420E-01  1.16690E-01
174 -6.39320E-01  1.54070E-01 -1.32890E-01
175 -2.06060E-01  1.24610E-01  1.18970E-01
176 -2.04570E-01  1.54070E-01  6.20110E-01
177 -5.40210E-01  1.51070E-01  3.11890E-01
178 -6.39330E-01  1.70140E-01  4.84610E-01
179 -7.39350E-01  1.70140E-01  3.11360E-01
180 -2.17370E-01  2.20110E-01 -1.11124E+00
181  4.34740E-01  1.80750E-01 -7.36800E-01
182 -2.17370E-01  1.39340E-01 -3.62580E-01
183 -8.55460E-01  1.80750E-01 -8.10000E-03
184 -2.05320E-01  1.39340E-01  3.69540E-01
185 -8.53670E-01  2.20110E-01  7.43870E-01
190 -4.34740E-01  5.10000E-02 -1.62940E-01
191 -4.34740E-01 -5.10000E-02 -1.62940E-01
192  4.34740E-01  5.10000E-02 -1.62940E-01
193  4.34740E-01 -5.10000E-02 -1.62940E-01
194 -1.31000E-02  3.24000E-02 -6.41000E-01
195 -4.34740E-01 -5.10000E-02 -1.62940E-01
196 -4.34740E-01  5.10000E-02 -1.62940E-01
197  4.34740E-01  5.10000E-02 -1.62940E-01
198  4.34740E-01 -5.10000E-02 -1.62940E-01
199  8.53670E-01  2.20110E-01  7.43870E-01
200  8.55460E-01  1.80750E-01 -8.10000E-03
201  2.05320E-01  1.39340E-01  3.69540E-01
202  0.00000E+00  3.25485E+00  0.00000E+00 %SM secondary mirror
203  3.81596E+00  0.00000E+00 -6.25000E-01
204  3.81596E+00  0.00000E+00  6.25000E-01
205 -3.81596E+00  0.00000E+00  6.25000E-01
206 -3.81596E+00  0.00000E+00 -6.25000E-01
207  0.00000E+00 -5.25000E-01 -3.50000E-01 %instrument (where back end optics and detector are) (where cryocooler is)
208  1.14929E+00 -1.21409E+00  2.41295E+00
209  1.68262E+00 -1.21409E+00  2.41295E+00
210  2.21596E+00 -1.21409E+00  2.41295E+00
211  2.74929E+00 -1.21409E+00  2.41295E+00
212  3.28262E+00 -1.21409E+00  2.41295E+00
213 -1.14929E+00 -1.21409E+00  2.41295E+00
214 -1.68262E+00 -1.21409E+00  2.41295E+00
215 -2.21596E+00 -1.21409E+00  2.41295E+00
216 -2.74929E+00 -1.21409E+00  2.41295E+00
217 -3.28262E+00 -1.21409E+00  2.41295E+00
218  0.00000E+00 -1.81832E+00  1.17568E+00
219  9.90960E-01 -8.87120E-01  1.51461E+00
220 -9.90960E-01 -8.87120E-01  1.51461E+00
221 -4.29550E-01  3.14550E+00 -2.48000E-01
222  0.00000E+00  3.14550E+00  4.96000E-01
223  4.29550E-01  3.14550E+00 -2.48000E-01
224  0.00000E+00  3.49470E-01  1.96469E+00
225  0.00000E+00 -2.52309E+00  9.19160E-01
226  0.00000E+00 -3.22785E+00  6.62650E-01
227  0.00000E+00 -3.93262E+00  4.06130E-01
228  1.74096E+00 -8.87120E-01  1.51461E+00
229  2.49096E+00 -8.87120E-01  1.51461E+00
230  3.24096E+00 -8.87120E-01  1.51461E+00
231 -1.74096E+00 -8.87120E-01  1.51461E+00
232 -2.49096E+00 -8.87120E-01  1.51461E+00
233 -3.24096E+00 -8.87120E-01  1.51461E+00
234  0.00000E+00  1.66504E+00  2.44352E+00
235  0.00000E+00  2.98061E+00  2.92235E+00
236  0.00000E+00  4.29618E+00  3.40117E+00
237 -4.87500E-01  7.19500E-01  1.00500E+00
238 -3.25000E-01  1.43900E+00  9.10000E-01
239 -1.62500E-01  2.15850E+00  8.15000E-01
240  4.87500E-01  7.19500E-01  1.00500E+00
241  3.25000E-01  1.43900E+00  9.10000E-01
242  1.62500E-01  2.15850E+00  8.15000E-01
243  1.00938E+00  7.19500E-01 -1.37250E-01
244  8.80770E-01  1.43900E+00 -2.11500E-01
245  7.52150E-01  2.15850E+00 -2.85750E-01
246  6.43380E-01  7.19500E-01 -7.70250E-01
247  6.36770E-01  1.43900E+00 -6.33500E-01
248  6.30150E-01  2.15850E+00 -4.96750E-01
249 -1.00938E+00  7.19500E-01 -1.37250E-01
250 -8.80770E-01  1.43900E+00 -2.11500E-01
251 -7.52150E-01  2.15850E+00 -2.85750E-01
252 -6.43380E-01  7.19500E-01 -7.70250E-01
253 -6.36770E-01  1.43900E+00 -6.33500E-01
254 -6.30150E-01  2.15850E+00 -4.96750E-01
255  6.24470E-01 -2.31230E-01  1.56896E+00
256  6.32980E-01 -1.54160E-01  1.41264E+00
257  6.41490E-01 -7.71000E-02  1.25632E+00
258 -6.24470E-01 -2.31230E-01  1.56896E+00
259 -6.32980E-01 -1.54160E-01  1.41264E+00
260 -6.41490E-01 -7.71000E-02  1.25632E+00
261  2.99470E-01 -2.31230E-01  1.56896E+00
262 -1.70000E-02 -1.54160E-01  1.41264E+00
263 -3.33510E-01 -7.71000E-02  1.25632E+00
264 -2.99470E-01 -2.31230E-01  1.56896E+00
265  1.70210E-02 -1.54160E-01  1.41264E+00
266  3.33510E-01 -7.71000E-02  1.25632E+00
267 -6.34470E-01 -1.09945E+00  9.40450E-01
268 -6.52980E-01 -7.32970E-01  5.76970E-01
269 -6.71490E-01 -3.66480E-01  2.13480E-01
270  6.34470E-01 -1.09945E+00  9.40450E-01
271  6.52980E-01 -7.32970E-01  5.76970E-01
272  6.71490E-01 -3.66480E-01  2.13480E-01
273 -6.24470E-01 -1.09945E+00  1.25295E+00
274 -6.32980E-01 -7.32970E-01  1.20197E+00
275 -6.41490E-01 -3.66480E-01  1.15098E+00
276  6.24470E-01 -1.09945E+00  1.25295E+00
277  6.32980E-01 -7.32970E-01  1.20197E+00
278  6.41490E-01 -3.66480E-01  1.15098E+00 
301  0.00000E+00 -5.25000E-01 -3.50000E-01 %new pt - instrument nicelas
302  6.15960E-01 -1.46593E+00  1.30394E+00 %new pt - 64 S/C top corner (4pts)
303  6.15960E-01 -3.08310E-01  1.72528E+00 %new pt - 27 S/C top corner (4pts)
304 -6.15960E-01 -3.08310E-01  1.72528E+00 %new pt - 28 S/C top corner (4pts)
305 -6.15960E-01 -1.46593E+00  1.30394E+00 %new pt - 30 S/C top corner (4pts)
];
	 
%====================================================
% 2. Local to Basic coordinate system transformations
%====================================================
% ci and cf are for coord_in

% Create coord system index matrix: 
% grid id, these correspond to the nodes in xyz
% We're identifying a set of coordinate system transformations, which are defined by cf
% each transformation goes from the input csid to the output csid
ci=[ ...
14       0      79
15       0      78
16       0      80
17       0      75
18       0      77
19       0      76
20       0      78
21       0      77
23       0      79
24       0      80
25       0      76
26       0      75
29       0      16
31       0      78
32       0      77
33       0      76
34       0      75
35       0      80
36       0      79
37       0      79
38       0      78
39       0      77
40       0      76
41       0      75
42       0      80
51       3       0
52       3       0
53       3       0
55       3       0
56       3       0
57       3       0
58       3       0
73       0      10
74       0      10
75       0      12
76       0      12
77       0       8
78       0       8
85       0      10
86       0      10
87       0      12
88       0      12
89       0       8
90       0       8
103       0       9
104       0       9
105       0      11
106       0      11
107       0      13
108       0      13
109       0      11
110       0      13
111       0       9
112       0      11
113       0      13
114       0       9
115       0      15
122       0      15
123       0      14
124       0      14
125       0      16
126       0      17
127       0      18
128       0      19
131       0      14
132       0      14
133       0      15
134       0      15
135       0      16
136       0      16
137       0      17
138       0      19
139       0      18
140       0      22
141       0      22
142       0      25
143       0      25
144       0      24
145       0      24
146       0      29
147       0      31
148       0      30
151       0      24
152       0      25
153       0      25
154       0      22
155       0      22
156       0      24
157       0      29
158       0      30
159       0      31
160       0      26
161       0      28
162       0      28
163       0      27
164       0      27
165       0      26
166       0      34
167       0      33
168       0      32
171       0      27
172       0      27
173       0      28
174       0      28
175       0      26
176       0      26
177       0      34
178       0      32
179       0      33
203       3       0
204       3       0
205       3       0
206       3       0 ];
	 
	 
% Create coordinate system definition matrix
%TODO figure out what exactly cf is
%   cf(:,1) is the coordinate system id.
%   cf(:,2) is the coordinate system type.
%1= rectangular (x,y,z) 2= cylindrical (r,theta,z) 3= spherical (r,theta,phi) (theta from z-axis and phi from x-axis in xy plane)
%   cf(:,3) is the reference coordinate system id used for a,b,c point locations (0=basic, meaning the one that xyz is defined in).
%          a= the origin of the new coordinate system.
%          b= a point on the z-axis of the new coordinate system.
%          c= a point in the +zx quadrant.
%   cf(:,4) cf(:,5) cf(:,6) is the 3 coordinates of point a in reference coord system.
%   cf(:,7) cf(:,8) cf(:,9) is the 3 coordinates of point b in reference coord system.
%   cf(:,10) cf(:,11) cf(:,12) is the 3 coordinates of point c in reference coord system.
cf=[ ...
1       2       0  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.10000E+01  0.10000E+01  0.00000E+00  0.10000E+01
2       3       0  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.10000E+01  0.10000E+01  0.00000E+00  0.10000E+01
3       1       0  0.00000E+00 -0.12141E+01  0.24129E+01  0.00000E+00 -0.27440E+00  0.27550E+01  0.10000E+01 -0.12141E+01  0.24129E+01
4       1       0 -0.47413E+00 -0.81375E+00  0.20447E+01 -0.14113E+01 -0.48596E+00  0.21640E+01 -0.47413E+00 -0.47173E+00  0.11050E+01
5       1       0 -0.47413E+00 -0.96331E+00  0.17239E+01 -0.14113E+01 -0.84400E+00  0.13961E+01 -0.47413E+00 -0.19030E+01  0.13819E+01
6       1       0 -0.47413E+00 -0.12840E+01  0.18735E+01 -0.14113E+01 -0.16118E+01  0.17542E+01 -0.47413E+00 -0.16261E+01  0.28132E+01
7       1       0 -0.47413E+00 -0.11345E+01  0.21942E+01 -0.14113E+01 -0.12538E+01  0.25220E+01 -0.47413E+00 -0.19479E+00  0.25362E+01
8       1       0 -0.40438E+00 -0.11357E+01  0.14667E+01 -0.95210E+00 -0.41776E+00  0.18962E+01 -0.11115E+01 -0.12585E+01  0.77031E+00
9       1       0 -0.40438E+00 -0.11357E+01  0.14667E+01 -0.95210E+00 -0.16635E+01  0.21159E+01  0.30273E+00 -0.10129E+01  0.21630E+01
10       1       0 -0.40438E+00 -0.57905E+00  0.21301E+01 -0.95210E+00 -0.13100E+01  0.25371E+01 -0.11115E+01  0.85413E-01  0.23719E+01
11       1       0 -0.40438E+00 -0.57905E+00  0.21301E+01 -0.95210E+00 -0.87743E+00  0.13484E+01  0.30273E+00 -0.12435E+01  0.18882E+01
12       1       0 -0.40438E+00 -0.14319E+01  0.22805E+01 -0.95210E+00 -0.14189E+01  0.14439E+01 -0.11115E+01 -0.19736E+01  0.27350E+01
13       1       0 -0.40438E+00 -0.14319E+01  0.22805E+01 -0.95210E+00 -0.60580E+00  0.24129E+01  0.30273E+00 -0.89024E+00  0.18259E+01
14       1       0 -0.43474E+00  0.20742E+00 -0.98638E+00  0.65255E-01  0.11713E+00 -0.12507E+00  0.43128E+00  0.25955E+00 -0.14837E+01
15       1       0  0.43474E+00  0.20742E+00 -0.98638E+00 -0.56526E+00  0.20742E+00 -0.98638E+00  0.43474E+00  0.10316E+00  0.81700E-02
16       1       0  0.00000E+00  0.12461E+00 -0.23794E+00  0.50000E+00  0.21491E+00 -0.10992E+01 -0.86603E+00  0.17674E+00 -0.73521E+00
17       1       0  0.00000E+00  0.15107E+00 -0.62378E+00  0.00000E+00  0.46806E-01  0.37077E+00  0.10000E+01  0.15107E+00 -0.62378E+00
18       1       0  0.10003E+00  0.17014E+00 -0.79598E+00  0.96605E+00  0.22227E+00 -0.12933E+01 -0.39997E+00  0.26043E+00 -0.16573E+01
19       1       0 -0.10003E+00  0.17014E+00 -0.79598E+00 -0.96605E+00  0.22227E+00 -0.12933E+01 -0.60003E+00  0.79849E-01  0.65326E-01
22       1       0  0.20606E+00  0.12461E+00  0.11897E+00  0.12020E+01  0.21490E+00  0.11661E+00  0.20842E+00  0.72476E-01 -0.87967E+00
24       1       0  0.10716E+01  0.20742E+00  0.11669E+00  0.57160E+00  0.20742E+00  0.98272E+00  0.19329E+01  0.31168E+00  0.61397E+00
25       1       0  0.63686E+00  0.20742E+00  0.86969E+00  0.14095E+00  0.11713E+00  0.60250E-02 -0.22681E+00  0.15529E+00  0.13711E+01
26       1       0 -0.20606E+00  0.12461E+00  0.11897E+00 -0.12020E+01  0.21490E+00  0.11661E+00 -0.20370E+00  0.17674E+00  0.11176E+01
27       1       0 -0.63686E+00  0.20742E+00  0.86969E+00 -0.14095E+00  0.11713E+00  0.60250E-02 -0.15005E+01  0.25955E+00  0.36833E+00
28       1       0 -0.10716E+01  0.20742E+00  0.11669E+00 -0.57160E+00  0.20742E+00  0.98272E+00 -0.21030E+00  0.10316E+00 -0.38058E+00
29       1       0  0.54021E+00  0.15107E+00  0.31189E+00 -0.32110E+00  0.46806E-01 -0.18539E+00  0.40209E-01  0.15107E+00  0.11779E+01
30       1       0  0.73935E+00  0.17014E+00  0.31136E+00  0.16030E+01  0.22227E+00 -0.19000E+00  0.24344E+00  0.79849E-01 -0.55230E+00
31       1       0  0.63933E+00  0.17014E+00  0.48461E+00  0.63697E+00  0.22227E+00  0.14832E+01  0.16352E+01  0.26043E+00  0.48225E+00
32       1       0 -0.63933E+00  0.17014E+00  0.48461E+00 -0.63697E+00  0.22227E+00  0.14832E+01  0.35659E+00  0.79849E-01  0.48697E+00
33       1       0 -0.73935E+00  0.17014E+00  0.31136E+00 -0.16030E+01  0.22227E+00 -0.19000E+00 -0.12353E+01  0.26043E+00  0.11750E+01
34       1       0 -0.54021E+00  0.15107E+00  0.31189E+00  0.32110E+00  0.46806E-01 -0.18539E+00 -0.10402E+01  0.15107E+00 -0.55414E+00
75       1       0  0.20197E+00  0.34593E+01 -0.10217E+00 -0.65717E+00  0.35851E+01 -0.59820E+00 -0.26112E+00  0.28557E+01  0.54683E+00
76       1       0  0.12500E-01  0.34593E+01  0.22600E+00  0.87164E+00  0.33335E+01  0.72203E+00  0.34301E+00  0.28557E+01 -0.49955E+00
77       1       0 -0.12500E-01  0.34593E+01  0.22600E+00  0.84664E+00  0.35851E+01 -0.27003E+00 -0.34301E+00  0.28557E+01 -0.49955E+00
78       1       0 -0.20197E+00  0.34593E+01 -0.10217E+00 -0.10611E+01  0.33335E+01  0.39385E+00  0.26112E+00  0.28557E+01  0.54683E+00
79       1       0 -0.18947E+00  0.34593E+01 -0.12383E+00 -0.18947E+00  0.35851E+01  0.86823E+00  0.60413E+00  0.28557E+01 -0.47300E-01
80       1       0  0.18947E+00  0.34593E+01 -0.12383E+00  0.18947E+00  0.33335E+01 -0.11159E+01 -0.60413E+00  0.28557E+01 -0.47300E-01
];
  
%==============================================================
% 3. Create FEM Element connectivity (nodal incidence matrices)
%==============================================================
% Create bar element connectivity matrix
%nibar is the (:,5), (:,7), or (:,10) nodal connectivity array, each row describes an element:
%   ni(:,1) is the element number
%   ni(:,2) is the node number of one end of the element
%   ni(:,3) is the node number of the other end of the element
%   ni(:,4) is the node number of a third node defining the plane for in-plane and out-of-plane bending
%   ni(:,5) is the property number, or row index into prop
%   ni(:,6) is the optional pin flag for the first end, ni(:,2).
%   ni(:,7) is the optional pin flag for the second end, ni(:,3).
%   ni(:,8) is the optional first component of the orient vector
%   ni(:,9) is the optional second component of the orient vector 
%   ni(:,10) is the optional third component of the orient vector 
nibar=[ ...
1      46     218       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
2      43     219       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
3      45     220       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
4       7     237       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
5       8     240       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
6      12     243       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
7      11     246       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
8       9     249       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
9      10     252       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
14      27     255       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
26      28     258       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
94       3     221       0       1       0       0  0.00000E+00  1.00000E+00  0.00000E+00
95       1     222       0       1       0       0  0.00000E+00  1.00000E+00  0.00000E+00
96       5     223       0       1       0       0  0.00000E+00  1.00000E+00  0.00000E+00
97      44     224       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
98      48     208       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
99      50     213       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
108      59     225       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
109      60     226       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
110      61     227       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
111      62     228       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
112      63     229       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
113      65     230       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
114      66     231       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
115      27     261       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
116      28     264       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
117      92     232       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
118      93     233       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
119      94     234       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
120      95     235       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
121      30     267       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
122      64     270       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
131      30     273       0      14       0       0  0.00000E+00  1.00000E+00  0.00000E+00
132      64     276       0      14       0       0  0.00000E+00  1.00000E+00  0.00000E+00
134      96     236       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
415     208     209       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
416     209     210       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
417     210     211       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
418     211     212       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
419     212      57       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
420     213     214       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
421     214     215       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
422     215     216       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
423     216     217       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
424     217      56       0      22       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
425     218      59       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
426     219      62       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
427     220      66       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
428     221       4       0       1       0       0  0.00000E+00  1.00000E+00  0.00000E+00
429     222       2       0       1       0       0  0.00000E+00  1.00000E+00  0.00000E+00
430     223       6       0       1       0       0  0.00000E+00  1.00000E+00  0.00000E+00
431     224      94       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
432     225      60       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
433     226      61       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
434     227      52       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
435     228      63       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
436     229      65       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
437     230      51       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
438     231      92       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
439     232      93       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
440     233      53       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
441     234      95       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
442     235      96       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
443     236      55       0      25       0       0  0.00000E+00  3.42020E-01 -9.39690E-01
444     237     238       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
445     238     239       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
446     239       2       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
447     240     241       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
448     241     242       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
449     242       2       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
450     243     244       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
451     244     245       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
452     245       6       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
453     246     247       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
454     247     248       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
455     248       6       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
456     249     250       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
457     250     251       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
458     251       4       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
459     252     253       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
460     253     254       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
461     254       4       0       2       0       0  0.00000E+00  1.00000E+00  0.00000E+00
462     255     256       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
463     256     257       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
464     257       8       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
465     258     259       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
466     259     260       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
467     260       7       0      11       0       0  0.00000E+00  1.00000E+00  0.00000E+00
468     261     262       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
469     262     263       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
470     263       7       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
471     264     265       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
472     265     266       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
473     266       8       0      12       0       0  0.00000E+00  1.00000E+00  0.00000E+00
474     267     268       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
475     268     269       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
476     269      67       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
477     270     271       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
478     271     272       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
479     272      68       0      13       0       0  0.00000E+00  1.00000E+00  0.00000E+00
480     273     274       0      14       0       0  0.00000E+00  1.00000E+00  0.00000E+00
481     274     275       0      14       0       0  0.00000E+00  1.00000E+00  0.00000E+00
482     275       7       0      14       0       0  0.00000E+00  1.00000E+00  0.00000E+00
483     276     277       0      14       0       0  0.00000E+00  1.00000E+00  0.00000E+00
484     277     278       0      14       0       0  0.00000E+00  1.00000E+00  0.00000E+00  
%501      83      79       0      99       0       0  0.00000E+00  1.00000E+00  0.00000E+00
%502      83      80       0      99       0       0  0.00000E+00  1.00000E+00  0.00000E+00
%503      83      81       0      99       0       0  0.00000E+00  1.00000E+00  0.00000E+00
%504      83      82       0      99       0       0  0.00000E+00  1.00000E+00  0.00000E+00 
  ];
  
% Create bar element property matrix
propbar=[ ...
%  pid      mid        A           I22           I33          J        nsm/elem     nsm/length        K1           K2           C1           C2           D1           D2           E1           E2           F1           F2
1       1  A_sptop      0.47840E-05  0.16260E-05  0.37970E-05  0.00000E+00  0.00000E+00  0.63662E+00  0.24174E+00 -0.80000E-01  0.40000E-01 -0.80000E-01 -0.40000E-01  0.80000E-01 -0.40000E-01  0.80000E-01  0.40000E-01
2       1  A_spbot      I_sp         I_sp         J_sp         0.00000E+00  0.00000E+00  0.53304E+00  0.53304E+00 -0.30000E-01  0.00000E+00  0.00000E+00 -0.30000E-01  0.30000E-01  0.00000E+00  0.00000E+00  0.30000E-01
11       1  0.63150E-03  0.35500E-06  0.35500E-06  0.71000E-06  0.00000E+00  0.00000E+00  0.53238E+00  0.53238E+00 -0.35000E-01  0.00000E+00  0.00000E+00 -0.35000E-01  0.35000E-01  0.00000E+00  0.00000E+00  0.35000E-01
12       1  0.44300E-03  0.12280E-06  0.12280E-06  0.24560E-06  0.00000E+00  0.00000E+00  0.53414E+00  0.53414E+00 -0.25000E-01  0.00000E+00  0.00000E+00 -0.25000E-01  0.25000E-01  0.00000E+00  0.00000E+00  0.25000E-01
13       1  0.63150E-03  0.35500E-06  0.35500E-06  0.71000E-06  0.00000E+00  0.00000E+00  0.53238E+00  0.53238E+00 -0.35000E-01  0.00000E+00  0.00000E+00 -0.35000E-01  0.35000E-01  0.00000E+00  0.00000E+00  0.35000E-01
14       1  0.63150E-03  0.35500E-06  0.35500E-06  0.71000E-06  0.00000E+00  0.00000E+00  0.53238E+00  0.53238E+00 -0.35000E-01  0.00000E+00  0.00000E+00 -0.35000E-01  0.35000E-01  0.00000E+00  0.00000E+00  0.35000E-01
21       2  0.91420E-03  I_ss         10*I_ss      2*I_ss       0.00000E+00  0.00000E+00  0.53149E+00  0.53149E+00 -0.50000E-01  0.00000E+00  0.00000E+00 -0.50000E-01  0.50000E-01  0.00000E+00  0.00000E+00  0.50000E-01
22       3  0.31250E-01  0.66630E-07  0.42660E-02  0.12400E-03  0.00000E+00  0.79400E+01  0.90170E+00  0.84969E+00 -0.12500E-01  0.62500E+00 -0.12500E-01 -0.62500E+00  0.12500E-01 -0.62500E+00  0.12500E-01  0.62500E+00
25       2  0.91420E-03  I_ss         10*I_ss      2*I_ss       0.00000E+00  0.00000E+00  0.53149E+00  0.53149E+00 -0.50000E-01  0.00000E+00  0.00000E+00 -0.50000E-01  0.50000E-01  0.00000E+00  0.00000E+00  0.50000E-01
65       2  0.91420E-03  0.64500E-09  0.64530E-08  0.15600E-07  0.00000E+00  0.00000E+00  0.53149E+00  0.53149E+00 -0.50000E-01  0.00000E+00  0.00000E+00 -0.50000E-01  0.50000E-01  0.00000E+00  0.00000E+00  0.50000E-01
66       2  0.91420E-03  0.64500E-09  0.64530E-08  0.15600E-07  0.00000E+00  0.00000E+00  0.53149E+00  0.53149E+00 -0.50000E-01  0.00000E+00  0.00000E+00 -0.50000E-01  0.50000E-01  0.00000E+00  0.00000E+00  0.50000E-01
%99       2  1.96350E-03  I_iso        I_iso        1.00000E-03  0.00000E+00  0.00000E+00  0.53149E+00  0.53149E+00 -0.50000E-01  0.00000E+00  0.00000E+00 -0.50000E-01  0.50000E-01  0.00000E+00  0.00000E+00  0.50000E-01
];

%stiffness matrix?? something to do with springs?
%ni(:,1) is the element number
%ni(:,2) is the node number of the first end of the element
%ni(:,3) is the component number for the first node (1-6)
%ni(:,4) is the node number of the second end of the element
%ni(:,5) is the component number for the second node (1-6)
%ni(:,6) is the value for the spring stiffness 
%ni(:,7) is the value for the structural damping factor (gfac, 0 is ok) (element structural damping is Ge= gfac*Ke)
nicelas2= [ ...
%e     node1  comp1 node2  comp2   K        damping(!)
21		15		1		38		1		K_act1		c_SM_act %15-38 secondary mirror actuator
22		15		2		38		2		K_act2		c_SM_act %15-38 secondary mirror actuator
23		15		3		38		3		K_act2		c_SM_act %15-38 secondary mirror actuator
24		15		4		38		4		K_rad1		c_SM_act %15-38 secondary mirror actuator
25		15		5		38		5		K_rad2		c_SM_act %15-38 secondary mirror actuator
27		15		6		38		6		K_rad2		c_SM_act %15-38 secondary mirror actuator
28		14		1		37		1		K_act1		c_SM_act %14-37 inner secondary mirror actuator
29		16		1		42		1		K_act1		c_SM_act %16-42 secondary mirror actuator
30		17		1		41		1		K_act1		c_SM_act %17-41 secondary mirror actuator
31		19		1		40		1		K_act1		c_SM_act %19-40 secondary mirror actuator
32		18		1		39		1		K_act1		c_SM_act %18-39 secondary mirror actuator
33		20		1		31		1		K_act1		c_SM_act %20-31 inner secondary mirror actuator
34		21		1		32		1		K_act1		c_SM_act %21-32 inner secondary mirror actuator
35		25		1		33		1		K_act1		c_SM_act %25-33 inner secondary mirror actuator
36		26		1		34		1		K_act1		c_SM_act %26-34 inner secondary mirror actuator
37		24		1		35		1		K_act1		c_SM_act %24-35 inner secondary mirror actuator
38		23		1		36		1		K_act1		c_SM_act %23-36 inner secondary mirror actuator
39		14		2		37		2		K_act2		c_SM_act %14-37 inner secondary mirror actuator
40		16		2		42		2		K_act2		c_SM_act %16-42 secondary mirror actuator
41		17		2		41		2		K_act2		c_SM_act %17-41 secondary mirror actuator
42		19		2		40		2		K_act2		c_SM_act %19-40 secondary mirror actuator
43		18		2		39		2		K_act2		c_SM_act %18-39 secondary mirror actuator
44		20		2		31		2		K_act2		c_SM_act %20-31 inner secondary mirror actuator
45		21		2		32		2		K_act2		c_SM_act %21-32 inner secondary mirror actuator
46		25		2		33		2		K_act2		c_SM_act %25-33 inner secondary mirror actuator
47		26		2		34		2		K_act2		c_SM_act %26-34 inner secondary mirror actuator
48		24		2		35		2		K_act2		c_SM_act %24-35 inner secondary mirror actuator
49		23		2		36		2		K_act2		c_SM_act %23-36 inner secondary mirror actuator
50		14		3		37		3		K_act2		c_SM_act %14-37 inner secondary mirror actuator
51		16		3		42		3		K_act2		c_SM_act %16-42 secondary mirror actuator
52		17		3		41		3		K_act2		c_SM_act %17-41 secondary mirror actuator
53		19		3		40		3		K_act2		c_SM_act %19-40 secondary mirror actuator
54		18		3		39		3		K_act2		c_SM_act %18-39 secondary mirror actuator
55		20		3		31		3		K_act2		c_SM_act %20-31 inner secondary mirror actuator
56		21		3		32		3		K_act2		c_SM_act %21-32 inner secondary mirror actuator
57		25		3		33		3		K_act2		c_SM_act %25-33 inner secondary mirror actuator
58		26		3		34		3		K_act2		c_SM_act %26-34 inner secondary mirror actuator
59		24		3		35		3		K_act2		c_SM_act %24-35 inner secondary mirror actuator
60		23		3		36		3		K_act2		c_SM_act %23-36 inner secondary mirror actuator
61		14		4		37		4		K_rad1		c_SM_act %14-37 inner secondary mirror actuator
62		16		4		42		4		K_rad1		c_SM_act %16-42 secondary mirror actuator
63		17		4		41		4		K_rad1		c_SM_act %17-41 secondary mirror actuator
64		19		4		40		4		K_rad1		c_SM_act %19-40 secondary mirror actuator
65		18		4		39		4		K_rad1		c_SM_act %18-39 secondary mirror actuator
66		20		4		31		4		K_rad1		c_SM_act %20-31 inner secondary mirror actuator
67		21		4		32		4		K_rad1		c_SM_act %21-32 inner secondary mirror actuator
68		25		4		33		4		K_rad1		c_SM_act %25-33 inner secondary mirror actuator
69		26		4		34		4		K_rad1		c_SM_act %26-34 inner secondary mirror actuator
70		24		4		35		4		K_rad1		c_SM_act %24-35 inner secondary mirror actuator
71		23		4		36		4		K_rad1		c_SM_act %23-36 inner secondary mirror actuator
72		14		5		37		5		K_rad2		c_SM_act %14-37 inner secondary mirror actuator
73		16		5		42		5		K_rad2		c_SM_act %16-42 secondary mirror actuator
74		17		5		41		5		K_rad2		c_SM_act %17-41 secondary mirror actuator
75		19		5		40		5		K_rad2		c_SM_act %19-40 secondary mirror actuator
76		18		5		39		5		K_rad2		c_SM_act %18-39 secondary mirror actuator
77		20		5		31		5		K_rad2		c_SM_act %20-31 inner secondary mirror actuator
78		21		5		32		5		K_rad2		c_SM_act %21-32 inner secondary mirror actuator
79		25		5		33		5		K_rad2		c_SM_act %25-33 inner secondary mirror actuator
80		26		5		34		5		K_rad2		c_SM_act %26-34 inner secondary mirror actuator
81		24		5		35		5		K_rad2		c_SM_act %24-35 inner secondary mirror actuator
82		23		5		36		5		K_rad2		c_SM_act %23-36 inner secondary mirror actuator
83		14		6		37		6		K_rad2		c_SM_act %14-37 inner secondary mirror actuator
84		16		6		42		6		K_rad2		c_SM_act %16-42 secondary mirror actuator
85		17		6		41		6		K_rad2		c_SM_act %17-41 secondary mirror actuator
86		19		6		40		6		K_rad2		c_SM_act %19-40 secondary mirror actuator
87		18		6		39		6		K_rad2		c_SM_act %18-39 secondary mirror actuator
88		20		6		31		6		K_rad2		c_SM_act %20-31 inner secondary mirror actuator
89		21		6		32		6		K_rad2		c_SM_act %21-32 inner secondary mirror actuator
90		25		6		33		6		K_rad2		c_SM_act %25-33 inner secondary mirror actuator
91		26		6		34		6		K_rad2		c_SM_act %26-34 inner secondary mirror actuator
92		24		6		35		6		K_rad2		c_SM_act %24-35 inner secondary mirror actuator
93		23		6		36		6		K_rad2		c_SM_act %23-36 inner secondary mirror actuator
141		78		2		90		2		K_yISO		c_RWA %78-90 RWA attach point
142		74		2		86		2		K_yISO		c_RWA %74-86 RWA attach point
143		76		2		88		2		K_yISO		c_RWA %76-88 RWA attach point
144		111		2		114		2		K_yISO		c_RWA %111-114 (78-90) RWA attach point
145		109		2		112		2		K_yISO		c_RWA %109-112 (74-86) RWA attach point
146		110		2		113		2		K_yISO		c_RWA %110-113 (76-88) RWA attach point
147		73		2		85		2		K_yISO		c_RWAI %73-85 bus point?
148		77		2		89		2		K_yISO		c_RWAI %77-89 bus point?
149		75		2		87		2		K_yISO		c_RWAI %75-87 bus point?
150		103		2		104		2		K_yISO		c_RWAI %103-104 (73-85) bus point?
151		107		2		108		2		K_yISO		c_RWAI %107-108 (77-89) bus point?
152		105		2		106		2		K_yISO		c_RWAI %105-106 (75-87) bus point?
153		78		1		90		1		K_xISO		c_RWA %78-90 RWA attach point
154		74		1		86		1		K_xISO		c_RWA %74-86 RWA attach point
155		76		1		88		1		K_xISO		c_RWA %76-88 RWA attach point
156		111		1		114		1		K_xISO		c_RWA %111-114 (78-90) RWA attach point
157		109		1		112		1		K_xISO		c_RWA %109-112 (74-86) RWA attach point
164		110		1		113		1		K_xISO		c_RWA %110-113 (76-88) RWA attach point
165		73		1		85		1		K_xISO		c_RWAI %73-85 bus point?
166		77		1		89		1		K_xISO		c_RWAI %77-89 bus point?
167		75		1		87		1		K_xISO		c_RWAI %75-87 bus point?
168		103		1		104		1		K_xISO		c_RWAI %103-104 (73-85) bus point?
169		107		1		108		1		K_xISO		c_RWAI %107-108 (77-89) bus point?
170		105		1		106		1		K_xISO		c_RWAI %105-106 (75-87) bus point?
171		78		3		90		3		K_xISO		c_RWA %78-90 RWA attach point
172		74		3		86		3		K_xISO		c_RWA %74-86 RWA attach point
173		76		3		88		3		K_xISO		c_RWA %76-88 RWA attach point
174		111		3		114		3		K_xISO		c_RWA %111-114 (78-90) RWA attach point
175		109		3		112		3		K_xISO		c_RWA %109-112 (74-86) RWA attach point
176		110		3		113		3		K_xISO		c_RWA %110-113 (76-88) RWA attach point
177		73		3		85		3		K_xISO		c_RWAI %73-85 bus point?
178		77		3		89		3		K_xISO		c_RWAI %77-89 bus point?
179		75		3		87		3		K_xISO		c_RWAI %75-87 bus point?
180		103		3		104		3		K_xISO		c_RWAI %103-104 (73-85) bus point?
181		107		3		108		3		K_xISO		c_RWAI %107-108 (77-89) bus point?
182		105		3		106		3		K_xISO		c_RWAI %105-106 (75-87) bus point?
183		78		4		90		4		K_rISO		c_RWA %78-90 RWA attach point
184		74		4		86		4		K_rISO		c_RWA %74-86 RWA attach point
185		76		4		88		4		K_rISO		c_RWA %76-88 RWA attach point
186		111		4		114		4		K_rISO		c_RWA %111-114 (78-90) RWA attach point
187		109		4		112		4		K_rISO		c_RWA %109-112 (74-86) RWA attach point
188		110		4		113		4		K_rISO		c_RWA %110-113 (76-88) RWA attach point
189		73		4		85		4		K_rISO		c_RWAI %73-85 bus point?
190		77		4		89		4		K_rISO		c_RWAI %77-89 bus point?
191		75		4		87		4		K_rISO		c_RWAI %75-87 bus point?
192		103		4		104		4		K_rISO		c_RWAI %103-104 (73-85) bus point?
193		107		4		108		4		K_rISO		c_RWAI %107-108 (77-89) bus point?
194		105		4		106		4		K_rISO		c_RWAI %105-106 (75-87) bus point?
195		78		5		90		5		K_aISO		c_RWA %78-90 RWA attach point
196		74		5		86		5		K_aISO		c_RWA %74-86 RWA attach point
197		76		5		88		5		K_aISO		c_RWA %76-88 RWA attach point
198		111		5		114		5		K_aISO		c_RWA %111-114 (78-90) RWA attach point
199		109		5		112		5		K_aISO		c_RWA %109-112 (74-86) RWA attach point
200		110		5		113		5		K_aISO		c_RWA %110-113 (76-88) RWA attach point
201		73		5		85		5		K_aISO		c_RWAI %73-85 bus point?
202		77		5		89		5		K_aISO		c_RWAI %77-89 bus point?
203		75		5		87		5		K_aISO		c_RWAI %75-87 bus point?
204		103		5		104		5		K_aISO		c_RWAI %103-104 (73-85) bus point?
205		107		5		108		5		K_aISO		c_RWAI %107-108 (77-89) bus point?
206		105		5		106		5		K_aISO		c_RWAI %105-106 (75-87) bus point?
207		78		6		90		6		K_rISO		c_RWA %78-90 RWA attach point
208		74		6		86		6		K_rISO		c_RWA %74-86 RWA attach point
209		76		6		88		6		K_rISO		c_RWA %76-88 RWA attach point
210		111		6		114		6		K_rISO		c_RWA %111-114 (78-90) RWA attach point
211		109		6		112		6		K_rISO		c_RWA %109-112 (74-86) RWA attach point
212		110		6		113		6		K_rISO		c_RWA %110-113 (76-88) RWA attach point
213		73		6		85		6		K_rISO		c_RWAI %73-85 bus point?
214		77		6		89		6		K_rISO		c_RWAI %77-89 bus point?
215		75		6		87		6		K_rISO		c_RWAI %75-87 bus point?
216		103		6		104		6		K_rISO		c_RWAI %103-104 (73-85) bus point?
217		107		6		108		6		K_rISO		c_RWAI %107-108 (77-89) bus point?
218		105		6		106		6		K_rISO		c_RWAI %105-106 (75-87) bus point?
228		131		1		124		1		K_pm1		c_PM %124-131 primary mirror
229		132		1		123		1		K_pm1		c_PM %123-132 primary mirror
230		133		1		122		1		K_pm1		c_PM %122-133 primary mirror
231		134		1		115		1		K_pm1		c_PM %115-134 primary mirror
232		135		1		29		1		K_pm1		c_PM %29-135 primary mirror
233		136		1		125		1		K_pm1		c_PM %125-136 primary mirror
234		131		2		124		2		K_yPM		c_PM %124-131 primary mirror
235		132		2		123		2		K_yPM		c_PM %123-132 primary mirror
236		133		2		122		2		K_yPM		c_PM %122-133 primary mirror
237		134		2		115		2		K_yPM		c_PM %115-134 primary mirror	
238		135		2		29		2		K_yPM		c_PM %29-135 primary mirror
239		136		2		125		2		K_yPM		c_PM %125-136 primary mirror
240		131		3		124		3		K_pm3		c_PM %124-131 primary mirror
241		132		3		123		3		K_pm3		c_PM %123-132 primary mirror
242		133		3		122		3		K_pm3		c_PM %122-133 primary mirror
243		134		3		115		3		K_pm3		c_PM %115-134 primary mirror
244		135		3		29		3		K_pm3		c_PM %29-135 primary mirror
245		136		3		125		3		K_pm3		c_PM %125-136 primary mirror
246		131		4		124		4		K_pm4		c_PM %124-131 primary mirror
247		132		4		123		4		K_pm4		c_PM %123-132 primary mirror
248		133		4		122		4		K_pm4		c_PM %122-133 primary mirror
249		134		4		115		4		K_pm4		c_PM %115-134 primary mirror
250		135		4		29		4		K_pm4		c_PM %29-135 primary mirror
251		136		4		125		4		K_pm4		c_PM %125-136 primary mirror
252		131		5		124		5		K_pm5		c_PM %124-131 primary mirror
253		132		5		123		5		K_pm5		c_PM %123-132 primary mirror
254		133		5		122		5		K_pm5		c_PM %122-133 primary mirror
255		134		5		115		5		K_pm5		c_PM %115-134 primary mirror
256		135		5		29		5		K_pm5		c_PM %29-135 primary mirror
257		136		5		125		5		K_pm5		c_PM %125-136 primary mirror
258		131		6		124		6		K_pm6		c_PM %124-131 primary mirror
259		132		6		123		6		K_pm6		c_PM %123-132 primary mirror
260		133		6		122		6		K_pm6		c_PM %122-133 primary mirror
261		134		6		115		6		K_pm6		c_PM %115-134 primary mirror
262		135		6		29		6		K_pm6		c_PM %29-135 primary mirror
263		136		6		125		6		K_pm6		c_PM %125-136 primary mirror
264		137		1		126		1		K_act2		c_PM %126-137 primary mirror actuator
265		138		1		128		1		K_act2		c_PM %128-138 primary mirror actuator
266		139		1		127		1		K_act2		c_PM %127-139 primary mirror actuator
267		137		2		126		2		K_act_pm2		c_SM_act %126-137 primary mirror actuator
268		138		2		128		2		K_act_pm2		c_SM_act %128-138 primary mirror actuator
269		139		2		127		2		K_act_pm2		c_SM_act %127-139 primary mirror actuator
270		137		3		126		3		K_act_pm3		c_SM_act %126-137 primary mirror actuator
271		138		3		128		3		K_act_pm3		c_SM_act %128-138 primary mirror actuator
272		139		3		127		3		K_act_pm3		c_SM_act %127-139 primary mirror actuator
273		137		4		126		4		K_act_pm4		c_SM_act %126-137 primary mirror actuator
274		138		4		128		4		K_act_pm4		c_SM_act %128-138 primary mirror actuator
275		139		4		127		4		K_act_pm4		c_SM_act %127-139 primary mirror actuator
276		137		5		126		5		K_act_pm5		c_SM_act %126-137 primary mirror actuator
277		138		5		128		5		K_act_pm5		c_SM_act %128-138 primary mirror actuator
278		139		5		127		5		K_act_pm5		c_SM_act %127-139 primary mirror actuator
279		137		6		126		6		K_act_pm6		c_SM_act %126-137 primary mirror actuator
280		138		6		128		6		K_act_pm6		c_SM_act %128-138 primary mirror actuator
281		139		6		127		6		K_act_pm6		c_SM_act %127-139 primary mirror actuator
284		151		1		144		1		K_pm1		c_PM %144-151 primary mirror
285		152		1		143		1		K_pm1		c_PM %143-152 primary mirror
286		153		1		142		1		K_pm1		c_PM %142-153 primary mirror
287		154		1		141		1		K_pm1		c_PM %141-154 primary mirror
288		155		1		140		1		K_pm1		c_PM %140-155 primary mirror
289		156		1		145		1		K_pm1		c_PM %145-156 primary mirror
290		151		2		144		2		K_yPM		c_PM %144-151 primary mirror
291		152		2		143		2		K_yPM		c_PM %143-152 primary mirror
292		153		2		142		2		K_yPM		c_PM %142-153 primary mirror
293		154		2		141		2		K_yPM		c_PM %141-154 primary mirror
294		155		2		140		2		K_yPM		c_PM %140-155 primary mirror
295		156		2		145		2		K_yPM		c_PM %145-156 primary mirror
296		151		3		144		3		K_pm3		c_PM %144-151 primary mirror
297		152		3		143		3		K_pm3		c_PM %143-152 primary mirror
298		153		3		142		3		K_pm3		c_PM %142-153 primary mirror
299		154		3		141		3		K_pm3		c_PM %141-154 primary mirror
300		155		3		140		3		K_pm3		c_PM %140-155 primary mirror
301		156		3		145		3		K_pm3		c_PM %145-156 primary mirror
302		151		4		144		4		K_pm4		c_PM %144-151 primary mirror
303		152		4		143		4		K_pm4		c_PM %143-152 primary mirror
304		153		4		142		4		K_pm4		c_PM %142-153 primary mirror
305		154		4		141		4		K_pm4		c_PM %141-154 primary mirror
306		155		4		140		4		K_pm4		c_PM %140-155 primary mirror
307		156		4		145		4		K_pm4		c_PM %145-156 primary mirror
308		151		5		144		5		K_pm5		c_PM %144-151 primary mirror
309		152		5		143		5		K_pm5		c_PM %143-152 primary mirror
310		153		5		142		5		K_pm5		c_PM %142-153 primary mirror
311		154		5		141		5		K_pm5		c_PM %141-154 primary mirror
312		155		5		140		5		K_pm5		c_PM %140-155 primary mirror
313		156		5		145		5		K_pm5		c_PM %145-156 primary mirror
314		151		6		144		6		K_pm6		c_PM %144-151 primary mirror
315		152		6		143		6		K_pm6		c_PM %143-152 primary mirror
316		153		6		142		6		K_pm6		c_PM %142-153 primary mirror
317		154		6		141		6		K_pm6		c_PM %141-154 primary mirror
318		155		6		140		6		K_pm6		c_PM %140-155 primary mirror
319		156		6		145		6		K_pm6		c_PM %145-156 primary mirror
320		157		1		146		1		K_act2		c_SM_act %146-157 primary mirror actuator
321		158		1		148		1		K_act2		c_SM_act %148-158 primary mirror actuator
322		159		1		147		1		K_act2		c_SM_act %147-159 primary mirror actuator
323		157		2		146		2		K_act_pm2		c_SM_act %146-157 primary mirror actuator
324		158		2		148		2		K_act_pm2		c_SM_act %148-158 primary mirror actuator
325		159		2		147		2		K_act_pm2		c_SM_act %147-159 primary mirror actuator
326		157		3		146		3		K_act_pm3		c_SM_act %146-157 primary mirror actuator
327		158		3		148		3		K_act_pm3		c_SM_act %148-158 primary mirror actuator
328		159		3		147		3		K_act_pm3		c_SM_act %147-159 primary mirror actuator
329		157		4		146		4		K_act_pm4		c_SM_act %146-157 primary mirror actuator
330		158		4		148		4		K_act_pm4		c_SM_act %148-158 primary mirror actuator
331		159		4		147		4		K_act_pm4		c_SM_act %147-159 primary mirror actuator
332		157		5		146		5		K_act_pm5		c_SM_act %146-157 primary mirror actuator
333		158		5		148		5		K_act_pm5		c_SM_act %148-158 primary mirror actuator
334		159		5		147		5		K_act_pm5		c_SM_act %147-159 primary mirror actuator
335		157		6		146		6		K_act_pm6		c_SM_act %146-157 primary mirror actuator
336		158		6		148		6		K_act_pm6		c_SM_act %148-158 primary mirror actuator
337		159		6		147		6		K_act_pm6		c_SM_act %147-159 primary mirror actuator
340		171		1		164		1		K_pm1		c_PM %171-164 primary mirror
341		172		1		163		1		K_pm1		c_PM %163-172 primary mirror
342		173		1		162		1		K_pm1		c_PM %162-173 primary mirror
343		174		1		161		1		K_pm1		c_PM %161-174 primary mirror
344		175		1		160		1		K_pm1		c_PM %160-175 primary mirror
345		176		1		165		1		K_pm1		c_PM %165-176 primary mirror
346		171		2		164		2		K_yPM		c_PM %164-171 primary mirror
347		172		2		163		2		K_yPM		c_PM %163-172 primary mirror
348		173		2		162		2		K_yPM		c_PM %162-173 primary mirror
349		174		2		161		2		K_yPM		c_PM %161-174 primary mirror
350		175		2		160		2		K_yPM		c_PM %160-175 primary mirror
351		176		2		165		2		K_yPM		c_PM %165-176 primary mirror
352		171		3		164		3		K_pm3		c_PM %164-171 primary mirror
353		172		3		163		3		K_pm3		c_PM %163-172 primary mirror
354		173		3		162		3		K_pm3		c_PM %162-173 primary mirror
355		174		3		161		3		K_pm3		c_PM %161-174 primary mirror
356		175		3		160		3		K_pm3		c_PM %160-175 primary mirror
357		176		3		165		3		K_pm3		c_PM %165-176 primary mirror
358		171		4		164		4		K_pm4		c_PM %164-171 primary mirror
359		172		4		163		4		K_pm4		c_PM %163-172 primary mirror
360		173		4		162		4		K_pm4		c_PM %162-173 primary mirror
361		174		4		161		4		K_pm4		c_PM %161-174 primary mirror
362		175		4		160		4		K_pm4		c_PM %160-175 primary mirror
363		176		4		165		4		K_pm4		c_PM %165-176 primary mirror
364		171		5		164		5		K_pm5		c_PM %164-171 primary mirror
365		172		5		163		5		K_pm5		c_PM %163-172 primary mirror
366		173		5		162		5		K_pm5		c_PM %162-173 primary mirror
367		174		5		161		5		K_pm5		c_PM %161-174 primary mirror
368		175		5		160		5		K_pm5		c_PM %160-175 primary mirror
369		176		5		165		5		K_pm5		c_PM %165-176 primary mirror
370		171		6		164		6		K_pm6		c_PM %164-171 primary mirror
371		172		6		163		6		K_pm6		c_PM %163-172 primary mirror
372		173		6		162		6		K_pm6		c_PM %162-173 primary mirror
373		174		6		161		6		K_pm6		c_PM %161-174 primary mirror
374		175		6		160		6		K_pm6		c_PM %160-175 primary mirror
375		176		6		165		6		K_pm6		c_PM %165-176 primary mirror
376		177		1		166		1		K_act2		c_SM_act %166-177 primary mirror actuator
377		178		1		168		1		K_act2		c_SM_act %168-170 primary mirror actuator
378		179		1		167		1		K_act2		c_SM_act %167-179 primary mirror actuator
379		177		2		166		2		K_act_pm2		c_SM_act %166-177 primary mirror actuator
380		178		2		168		2		K_act_pm2		c_SM_act %168-170 primary mirror actuator
381		179		2		167		2		K_act_pm2		c_SM_act %167-179 primary mirror actuator
382		177		3		166		3		K_act_pm3		c_SM_act %166-177 primary mirror actuator
383		178		3		168		3		K_act_pm3		c_SM_act %168-170 primary mirror actuator
384		179		3		167		3		K_act_pm3		c_SM_act %167-179 primary mirror actuator
385		177		4		166		4		K_act_pm4		c_SM_act %166-177 primary mirror actuator
386		178		4		168		4		K_act_pm4		c_SM_act %168-170 primary mirror actuator
387		179		4		167		4		K_act_pm4		c_SM_act %167-179 primary mirror actuator
388		177		5		166		5		K_act_pm5		c_SM_act %166-177 primary mirror actuator
389		178		5		168		5		K_act_pm5		c_SM_act %168-170 primary mirror actuator
390		179		5		167		5		K_act_pm5		c_SM_act %167-179 primary mirror actuator
391		177		6		166		6		K_act_pm6		c_SM_act %166-177 primary mirror actuator
392		178		6		168		6		K_act_pm6		c_SM_act %168-170 primary mirror actuator
393		179		6		167		6		K_act_pm6		c_SM_act %167-179 primary mirror actuator
395		197		3		192		3		K_zpet		c_petal %192-197 petal hinge?
396		196		3		190		3		K_zpet		c_petal %190-196 petal hinge?
397		198		3		193		3		K_zpet		c_petal %193-198 petal hinge?
398		195		3		191		3		K_zpet		c_petal %191-195 petal hinge?
399		198		1		193		1		K_xpet		c_petal %193-198 petal hinge?
400		198		2		193		2		K_xpet		c_petal %193-198 petal hinge?
401		195		2		191		2		K_xpet		c_petal %191-195 petal hinge?
402		301		1		207		1		K_cryo 		c_cryo %207-301 new pt - instrument nicelas
403		301		2		207		2		K_cryo 		c_cryo %207-301 new pt - instrument nicelas
404		301		3		207		3		K_cryo 		c_cryo %207-301 new pt - instrument nicelas
405		301		4		207		4		K_cryo 		c_cryo %207-301 new pt - instrument nicelas
406		301		5		207		5		K_cryo 		c_cryo %207-301 new pt - instrument nicelas
407		301		6		207		6		K_cryo 		c_cryo %207-301 new pt - instrument nicelas
408		302		1		64		1		K_IA		c_IA %64-302 S/C top corner (4pts)
409		302		2		64		2		K_IA		c_IA %64-302 S/C top corner (4pts)
410		302		3		64		3		K_IA		c_IA %64-302 S/C top corner (4pts)
411		302		4		64		4		K_IA		c_IA %64-302 S/C top corner (4pts)
412		302		5		64		5		K_IA		c_IA %64-302 S/C top corner (4pts)
413		302		6		64		6		K_IA		c_IA %64-302 S/C top corner (4pts)
414		303		1		27		1		K_IA		c_IA %27-303 S/C top corner (4pts)
415		303		2		27		2		K_IA		c_IA %27-303 S/C top corner (4pts)
416		303		3		27		3		K_IA		c_IA %27-303 S/C top corner (4pts)
417		303		4		27		4		K_IA		c_IA %27-303 S/C top corner (4pts)
418		303		5		27		5		K_IA		c_IA %27-303 S/C top corner (4pts)
419		303		6		27		6		K_IA		c_IA %27-303 S/C top corner (4pts)
420		304		1		28		1		K_IA		c_IA %28-304 S/C top corner (4pts)
421		304		2		28		2		K_IA		c_IA %28-304 S/C top corner (4pts)
422		304		3		28		3		K_IA		c_IA %28-304 S/C top corner (4pts)
423		304		4		28		4		K_IA		c_IA %28-304 S/C top corner (4pts)
424		304		5		28		5		K_IA		c_IA %28-304 S/C top corner (4pts)
425		304		6		28		6		K_IA		c_IA %28-304 S/C top corner (4pts)
426		305		1		30		1		K_IA		c_IA %30-305 S/C top corner (4pts)
427		305		2		30		2		K_IA		c_IA %30-305 S/C top corner (4pts)
428		305		3		30		3		K_IA		c_IA %30-305 S/C top corner (4pts)
429		305		4		30		4		K_IA		c_IA %30-305 S/C top corner (4pts)
430		305		5		30		5		K_IA		c_IA %30-305 S/C top corner (4pts)
431		305		6		30		6		K_IA		c_IA %30-305 S/C top corner (4pts)
];
  
%==============================================================
% 4. Define Concentrated Masses
%==============================================================

%           ni(:,1)  is the element id number
%           ni(:,2)  is the grid point id number 
%           ni(:,3)  is the mass value
%           -------------------------------------
%           ni(:,4)  is the 1st offset component from the grid point to the cg (optional)
%           ni(:,5)  is the 2nd offset component from the grid point to the cg (optional)
%           ni(:,6)  is the 3rd offset component from the grid point to the cg (optional)
%           -------------------------------------
%           ni(:,7)  is the I11 mass moment of inertia about the cg (optional, see note)
%           ni(:,8)  is the I22 mass moment of inertia about the cg (optional, see note)
%           ni(:,9)  is the I33 mass moment of inertia about the cg (optional, see note)
%           -------------------------------------
%           ni(:,10) is the I21 mass product of inertia about the cg (optional, see note)
%           ni(:,11) is the I31 mass product of inertia about the cg (optional, see note)
%           ni(:,12) is the I32 mass product of inertia about the cg (optional, see note)
niconm2= [ ...
12      22  m_SM  0.00000E+00  0.00000E+00  0.00000E+00  0.22300E-01  0.42990E-01  0.22300E-01  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror RBE2
13      13  m_SMhub  0.00000E+00  0.00000E+00  0.00000E+00  I_SMhubt  I_SMhuba  I_SMhubt  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror hub
100      45  0.57000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
101      44  0.57000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
102      43  0.57000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
103      46  0.57000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
104      51  0.10000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
105      53  0.10000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
106      52  0.10000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
107      55  0.10000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
135      84  m_bus  0.00000E+00  0.00000E+00  0.00000E+00  I_bus  I_bus  I_bus  0.00000E+00  0.00000E+00  0.00000E+00 %spacecraft bus CM
136      84  m_prop  0.00000E+00  0.00000E+00  0.00000E+00 I_propt  I_propa I_propt  0.00000E+00  0.00000E+00  0.00000E+00 %spacecraft bus CM
137      79  m_RW  0.00000E+00  0.00000E+00  0.00000E+00  I_RWt  0.14339E+00  I_RWt  0.00000E+00  0.00000E+00  0.00000E+00 %reaction wheel
138      80  m_RW  0.00000E+00  0.00000E+00  0.00000E+00  I_RWt  0.14339E+00  I_RWt  0.00000E+00  0.00000E+00  0.00000E+00 %reaction wheel
139      81  m_RW  0.00000E+00  0.00000E+00  0.00000E+00  I_RWt  0.14339E+00  I_RWt  0.00000E+00  0.00000E+00  0.00000E+00 %reaction wheel
140      82  m_RW  0.00000E+00  0.00000E+00  0.00000E+00  I_RWt  0.14339E+00  I_RWt  0.00000E+00  0.00000E+00  0.00000E+00 %reaction wheel
158     116  0.70000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.56450E-04  0.20250E-02  0.20250E-02  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror
159     117  0.70000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.56450E-04  0.20250E-02  0.20250E-02  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror
160     118  0.70000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.56450E-04  0.20250E-02  0.20250E-02  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror
161     119  0.70000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.56450E-04  0.20250E-02  0.20250E-02  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror
162     120  0.70000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.56450E-04  0.20250E-02  0.20250E-02  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror
163     121  0.70000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.56450E-04  0.20250E-02  0.20250E-02  0.00000E+00  0.00000E+00  0.00000E+00 %secondary mirror
219     102  m_ISO  0.00000E+00  0.00000E+00  0.00000E+00  I_ISOt  I_ISOa  I_ISOt  0.00000E+00  0.00000E+00  0.00000E+00 %something in bus?
220     101  m_ISO  0.00000E+00  0.00000E+00  0.00000E+00  I_ISOt  I_ISOa  I_ISOt  0.00000E+00  0.00000E+00  0.00000E+00 %something in bus?
221     100  m_ISO  0.00000E+00  0.00000E+00  0.00000E+00  I_ISOt  I_ISOa  I_ISOt  0.00000E+00  0.00000E+00  0.00000E+00 %something in bus?
222      98  m_ISO  0.00000E+00  0.00000E+00  0.00000E+00  I_ISOt  I_ISOa  I_ISOt  0.00000E+00  0.00000E+00  0.00000E+00 %something in bus?
223      99  m_ISO  0.00000E+00  0.00000E+00  0.00000E+00  I_ISOt  I_ISOa  I_ISOt  0.00000E+00  0.00000E+00  0.00000E+00 %something in bus?
224      97  m_ISO  0.00000E+00  0.00000E+00  0.00000E+00  I_ISOt  I_ISOa  I_ISOt  0.00000E+00  0.00000E+00  0.00000E+00 %something in bus?
225      83  m_RWAchx  0.00000E+00  0.00000E+00  0.00000E+00  I_xRWA  I_yRWA  I_zRWA  0.00000E+00  0.00000E+00  0.00000E+00 %center of RWA
227     130  0.10780E+02  0.00000E+00  0.00000E+00  0.00000E+00  0.10271E+01  0.20347E+01  0.10156E+01  0.00000E+00  0.00000E+00  0.00000E+00 %primary mirror
283     150  0.10780E+02  0.00000E+00  0.00000E+00  0.00000E+00  0.10271E+01  0.20347E+01  0.10156E+01  0.00000E+00  0.00000E+00  0.00000E+00 %primary mirror
339     170  0.10780E+02  0.00000E+00  0.00000E+00  0.00000E+00  0.10271E+01  0.20347E+01  0.10156E+01  0.00000E+00  0.00000E+00  0.00000E+00 %primary mirror
402      91  0.59640E+02  0.00000E+00  0.00000E+00  0.16800E+00  0.13340E+02  0.40470E+02  0.27530E+02  0.00000E+00  0.00000E+00  0.00000E+00 %origin point??
403     194  0.18860E+02  0.00000E+00  0.00000E+00  0.00000E+00  0.16190E+01  0.29690E+01  0.14760E+01  0.00000E+00  0.00000E+00  0.00000E+00 %primary mirror
404      52  0.49600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
405      53  0.49600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
406      51  0.49600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
407      55  0.92600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
408      45  0.49600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
409      46  0.49600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
410      43  0.49600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
411      44  0.92600E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 %sunshield
414     207  m_instr      0.00000E+00  0.00000E+00  0.00000E+00  I_i1         I_i2         I_i3  0.00000E+00  0.00000E+00  0.00000E+00 %instrument
  ];

%TODO what is this
cmcid2= [ ...
0
0
3
3
3
3
3
3
3
3
3
3
4
5
6
7
80
79
78
77
76
75
12
11
13
10
8
9
3
17
29
34
0
17
3
3
3
3
3
3
3
3
3
];

%Make sure xyz is sorted
[s,si]=sort(xyz(:,1));
xyz=xyz(si,:);
clear s si
  
%  compute cs tranformation matrices and xyz in basic
xyzi= xyz; % save xyz before coord_in processing as xyzi
%TODO what does coord_in do
[ti,tf,xyz]=coord_in(ci,cf,xyz);
  
%==============================================================
% 4. Create Material matrix
%==============================================================
%TODO this is ugly. fix this with more intelligent variable naming

mat=0;
mid=       1;
E=   0.89800E+11;
nu=  0.31000E+00;
rho=  0.16643E+04;
alpha=  0.56300E-07;
Tref=  0.00000E+00;
mat=mat1(mid,E,nu,rho,alpha,Tref,mat);
%TODO what does mat1 do
 
% WARNING: G-IMOS IS DIFFERENT THAN G-NASTRAN
 
mid=       2;
E=  0.89800E+11;
nu=  0.31000E+00;
rho=  0.00000E+00;
alpha=  0.56300E-07;
Tref=  0.00000E+00;
mat=mat1(mid,E,nu,rho,alpha,Tref,mat);
 
% WARNING: G-IMOS IS DIFFERENT THAN G-NASTRAN
 
mid=       3;
E=  0.89800E+11;
nu=  0.31000E+00;
rho=  0.00000E+00;
alpha=  0.56300E-07;
Tref=  0.00000E+00;
mat=mat1(mid,E,nu,rho,alpha,Tref,mat);
 
%==============================================================
% 5. Define Boundary Conditions and Initialize m and k
%==============================================================

% Create initial boundary condition matrix
%  (0= constrain, 1= free)
bci=ones(273,6);
  
 %  Compute the number of dofs and dof numbers, initialize k and m
 %TODO how does bcond work
[bc,ndof]=bcond(bci);
nset=[1:ndof];
fset=nset;
k=sparse(ndof,ndof);
m=k;
g=k;

%==============================================================
% 6. Assembly of m and k
%==============================================================
%TODO what does celas2 do
if diagnostics
disp('Assembling m for conm2s')
end
m=conm2(niconm2,bc,m,xyz,cmcid2,cf);

%TODO what does beam_lump do
if diagnostics
disp('Assembling k and m for beams')
end
[k,m]=beam_lump(nibar,xyz,propbar,mat,bc,k,m,ti,tf);

if diagnostics
disp('Transforming k and m to local cs')
end
[k,m]=km2loc(k,m,ti,tf);

%TODO what does conm2 do
if diagnostics
disp( 'Assembling k for spring elements')
end
%calculates spring element matrices (stiffness matrix, structural damping matrix) 
%from connection list, node list, boundary conditions for each dof
[k,g]=celas2(nicelas2,xyz,bc,k,g);
  
%==============================================================
% 7. Determine Mass Properties and RBM's
%==============================================================
% compute rigid body modes geometrically,
% where the rotation point is the c.g. and compute mass properties
%[xyzcg]=cg_calc(m,xyz,bc,ti,tf);
%TODO double-check how cg_calc works
[xyzcg]=cg_calc(m,xyz,bc);
%mass=wtcg(bc,xyz,m,xyzcg);
Rpt=xyzcg;
%TODO what does rbmodes do
phirb=rbmodes(bc,xyz,m,Rpt);
%TODO what does normphi do
phirb=normphi(phirb,m);

%==============================================================
% 8. Process Rigid Body elements
%==============================================================
%TODO this is all so ugly! needs to be rewritten, and better commenting, and maybe looped?
  rg=[0];
  mset=[0];
  nirbe=[];
  mapdof=[];
  for ind=1:size(xyz,1)
     ids(xyz(ind,1))=ind;
  end
  bci=bc;
  % secondary mirror hub
  if diagnostics
     disp('Processing RBE2 rigid body elements')
     end
gn=[      13];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
       3
      15
      14
      18
       1
      19
      16
      17
       5
    ];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
    nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
	%TODO what does rbe2 do, its all over this section!
    [nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);

% secondary mirror RBE2
gn=[      22];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      20
      21
      25
      26
      24
      23
     202
  ];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
    nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      38];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      31
     118
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      39];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      32
     119
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      40];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      33
     120
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      41];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      34
     121
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      42];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      35
     116
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      37];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      36
     117
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
  [nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
  % Spacecraft Rigid Element
gn=[      84]; % independent node S/C CM
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      27
      28
      30
      64
      71
      72
      69
      70
      86
      88
      90
      48
      47
      50
      49
      44
      45
      46
      43
      58
     112
     113
     114
  ];
  % 27,28,30,64 S/C top corner
  % 71,72,69,70 bottom 4 attach points to Launch Vehicle
  % 86(112) , 88(113), 90(114) RWA attach points
  % 48, 50 solar array attach
  % 47, 49 mid point, S/C bottom edge (pitch)
  % 43-46   sunshield
     % 58 origin (not CM)
     % omitted for clarity nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[     111];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     103
      97
   ];
   nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      74];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      73
      98
   ];
   nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      78];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      77
      99
   ];
   nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[     110];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     107
     100
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      76];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      75
     102
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[     109];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     105
     101
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
  [nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
  % RWA Chassis RBE
gn=[      83];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      81 % replace with axle stiffness
      82
      79
      80
      85
      87
      89
     104
     106
     108
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      91];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
      11
      12
       8
       7
       9
      10
      67
      68
     173
     174
     175
     176
     171
     172
     179
     177
     178
     155
     156
     151
     152
     153
     154
     157
     158
     159
     185
     183
     184
     196
     195
     197
     198
     200
     199
     201
     207
  ];
  % omitted for clarity nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[     130];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     125
     124
     123
     122
     115
      29
     126
     128
     127
     129
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[     150];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     145
     144
     143
     142
     141
     140
     146
     148
     147
     149
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[     170];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     165
     164
     163
     162
     161
     160
     166
     168
     167
     169
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[     194];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     190
     191
     192
     193
     135
     182
     136
     131
     180
     132
     133
     181
     134
     137
     138
     139
  ];
  % omit for clarity nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      56];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     205
     206
  ];
  nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
gn=[      57];
cm=[ ...
       1
       2
       3
       4
       5
       6
  ];
gm=[ ...
     204
     203
  ];
nirbe=[nirbe; [gn*ones(length(gm),1) gm]];
    mapdof=[mapdof; [reshape(bci(ids(gm),:)',length(gm)*6,1) repmat(bci(ids(gn),:)',length(gm),1)]];
[nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf);
rg=sparse(rg);

%TODO what does mce1 do
if diagnostics
 disp('Reducing k and m to independent dofs')
end
[gm,k,m]=mce1(nset,mset,rg,k,m);
bc=bcnset(bc,nset,mset);
fset=nset;
temp=reshape(bc',prod(size(bc)),1);
mapdof(:,2)=temp(mapdof(:,2));
  
% End of FEM assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if diagnostics
   disp(['Assembled NEXUS FEM in ' num2str(toc) ' [sec]'])
end
if diagnostics
   tic
end

%eigfem is for solving for the eigenvalues and eigenvectors of this man-body spring mass problem
[phi,omeg]=eigfem(k,m);
   %disp(['Solved Eigenproblem in ' num2str(toc) ' [sec]'])
%phi are the mass normalized modeshapes.
%omeg are the modal frequencies (rad/s).

omeg=abs(omeg);
phi=real(phi);

phi=normphi(phi,m); % mass normalize
   %disp(['onormphi in ' num2str(toc) ' [sec]'])

if diagnostics
   disp(['Solved Eigenproblem in ' num2str(toc) ' [sec]'])
   end

 phi=mce_dis(nset,mset,rg,phi);
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
%  define input dofs - this is important because it tells us where disturbances enter in
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

      
 %  define displacement output dofs - this is important because it tells us where the disturbances have an impact
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
% assume baseline zeta=0.005 for elastic structural modes
damping=ones(1,length(omeg)-3);
%damping([])=ones() %IEC damping
%damping([])=ones() %bus damping
%damping([])=ones() %MTMD damping
%damping([])=ones() %cryocooler CJAA damping
%damping([])=ones() %Cryocooler CCA damping
%damping([])=ones() %Isolator assembly damping
%damping([])=ones() %RWA damping
damping([4:8 16:18])=5*ones(1,8); % sunshield damping 5%
damping([9:11 14:15 19:25])=20*ones(1,12); % isolator damping 10% %TODO which isolator?
damping([12:13 26:27])=20*ones(1,4); % solar panel damping 5%
Damping=diag(damping);
za=zeta*damping; %I don't see how the above damping percents are happening, given that they're getting multiplied by 0.005?
% remove translational rigid body modes since unobservable/uncontrollable
phi=phi(:,4:nm+3);
omeg=omeg(4:nm+3);
nrbm=3;
%  compute state space model
[Ap,Bp,Cp,Dp,lb,lc] = mode2ss2(xyz,bc,nm,ig,dg,vg,nrbm,za,phi,omeg); %TODO this is probably important for me to understand notionally
sysp=ss(Ap,Bp,Cp,Dp);
O=diag(omeg);

if diagnostics
   disp(['Created NEXUS file in XVIEW format'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_optics: Optics?
mm2m=0.001;   % conversion from mm to m
m2micron=1e6; % conversion from m to micrometers (microns)
m2nm=1e9;     % conversion from m to nanometer
R=[0 1 0; 0 0 1; 1 0 0];  % rotation matrix from MACOS frame to FEM frame
% Note in FEM frame the +Y axis is optical boresight, in MACOS it is +Z

%TODO maybe someday in the future turn these into model inputs...
load nexopt   % previously computed with MACOS
%dcdu= 
%dwdu=
%k2ij=
%wnom=
   
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_acsloop: ACS system
Kcc=Kc*ones(1,3);         % FSM-to-ACS cross coupling gain
wc=2*pi*fca;
zetac=sqrt(2)/2;   % Butterworth pattern - critical damping
% Extract Inertia's from  6x6 rigid body mass matrix
% Find principal axes of inertia 
%load mass

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
	%TODO figure out if these are input parameters...
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nexus_fsm: fast-steering mirror

% variable parameters
FgsDes = 1/Tgs;	     % desired FSM sample rate [Hz] % guider

% other fixed parameters
%TODO figure out if any of these are input parameters
sfac = FgsDes/FgsNom;		% frequency scaling
Tg = 1/FgsDes;			      % GS sample period (sec)
wlz = sfac*2*pi*0.4;		   % lead zero
wlp = sfac*2*pi*5;		   % lead pole
wlpp1 = sfac*2*pi*0.1;		% low pass poles
wlpp2 = sfac*2*pi*0.01;
FsmIntGain=0*0.7071*0.1*2*pi;  % FSM integral gain set to zero

numc = Kcf*[1 wlz];
denc = conv([1 wlp],conv([1 wlpp1],[1 wlpp2]));

Acf=[-wlp 0 0; 0 0 1; Kcf -wlpp1*wlpp2 -(wlpp1+wlpp2)];
Bcf=[(wlz-wlp) 0 Kcf]';
Ccf=[0 1 0];
Dcf=[0];
syscf=ss(Acf,Bcf,Ccf,Dcf);
   % 2 parallel uncoupled control channels
   Acf=[Acf zeros(3,3); zeros(3,3) Acf];
   Bcf=[Bcf zeros(3,1); zeros(3,1) Bcf];
   Ccf=[Ccf zeros(1,3); zeros(1,3) Ccf]; 
   Dcf=[Dcf 0; 0 Dcf];

if plotflag
f=logspace(-1,3,1000);
[magcf1,phscf1]=bode(numc,denc,2*pi*f);
[magcf,phscf]=bode(syscf,2*pi*f); magcf=squeeze(magcf); phscf=squeeze(phscf);
subplot(211)
loglog(f,magcf,'b-',f,magcf1,'r-')
title('FSM controller')
subplot(212)
semilogx(f,phscf,'b-',f,phscf1,'r-')
xlabel('Frequency [Hz]')
ylabel('phase')
end


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
%See deWeck 2004 Equation 2
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------
% Solve Lyapunov Equation
%-------------------------------------------------------------

warning off
if diagnostics
   warning on
   tstartnewlyap=cputime;
end
syszd=ss(Azd,Bzd,Czd,Dzd); scf=0;
%TODO what does ss2mod7 do
[syszdm,T]=ss2mod7(syszd,scf);
blocksize=20;
[Azdm,Bzdm,Czdm,Dzdm]=ssdata(syszdm);
%TODO what does newlyap do
Sqm = newlyap(Azdm,Bzdm,blocksize);
Szm=Czdm*Sqm*Czdm';
z1=sqrt((1/nray)*sum(diag(Szm(1:nray,1:nray)))); %RMMS WFE
z2=sqrt(sum(diag((Szm(nray+1:nray+2,nray+1:nray+2))))); %RSS LOS
if diagnostics
   disp(['Performance Evaluation (newlyap) took: ' num2str(cputime-tstartnewlyap) ' seconds'])
   disp(['Performance Results: z1= ' num2str(z1) ' [nm], z2= ' num2str(z2) ' [\mum]'])
end
warning on

disp(['Total Runtime: ' num2str(cputime-tstart) ' [sec]'])

% Lyapunov simulation results
disp('Lyapunov simulation results')
RMMS_WFE_lyap=z1; % RMMS WFE (root-mean-mean-square wavefront error)
RSS_LOS_lyap=z2; %RSS LOS (root-sum-square line-of-sight)
return
end