clear all
close all

%%%These are the inputs that are likely part of theta
%i.e. they are uncertain inputs to both system model and experiment model
%{ 
RWA micro-vibe:
input.Us=1.8;       % Us:  static wheel imbalance         [gcm]    (%0.7160 test)
input.Ud=60;        % Ud:  dynamic wheel imbalance        [gcm^2]  (%29.536 test)

Cryocooler micro-vibe:
input.C -- this represents the force in [N] associated with cryocooler disturbance somehow?
input.Qc=0.005;     % Qc:  cryocooler attenuation factor  [-] 

Transmissibility testing:
input.zeta=0.005;         % zeta: Global modal damping ratio    [-]
	IEC harness zeta TODO
	spacecraft bus zeta TODO
	sunshield zeta TODO
	OTE zeta TODO
	MTMD zeta TODO
	RWIA zeta TODO
	CCA zeta TODO
	CJAA zeta TODO
	IA zeta TODO

Stiffness testing:
IEC stuffness?? TODO
sunshield stiffness?? TODO

Flow noise verification:
Turbulent flow effect?? TODO I need to make that
%} 

%Secondary paramters implicated by testing that might be part of theta:
%{ 
RWA micro-vibe:
input.m_RW; % RW mass (ITHACO E-wheel)
input.I_RWt; % Transverse Inertia of ITHACO E reaction wheel
input.I_RWa; % Axial Inertia of ITHACO E reaction wheel

Cryocooler micro-vibe:

IEC harness damping test:

spacecraft bus damping test:

sunshield damping test:
	
OTE damping test:

MTMD damping test:

RWIA damping test:
I_xRWA=x.I_xRWA;  % Ix moment of inertia of RWA chassis
I_yRWA=x.I_yRWA;  % Iy moment of inertia of RWA chassis
I_zRWA=x.I_zRWA; % Iz moment of inertia of RWA chassis
m_ISO=x.m_ISO;  % mass of RWA isolator strut actuator
I_ISOa=x.I_ISOa;  % RWA isolator strut axial inertia
I_ISOt=x.I_ISOt;  % RWA isolator strut transverse inertia

CCA damping test:

CJAA damping test:

IA damping test:

IEC stiffness testing:

sunshield stiffness testing:

Flow noise verification:
Turbulent flow effect?? TODO I need to make that
%} 

%%%These are the inputs that are likely internally-optimized control parameters
%{ 
input.fc=30;        % fc:  cryocooler drive frequency     [Hz]
The other drive frequency?? %TODO are the two compressors optional? as in, is one one in operation at a time?
FGC sample rate?? TODO is that in this model
input.Kcf=2000;    % Kcf: FSM controller gain             [-]
%}
%TODO what about input.Kc FSM-to-ACS cross coupling gain ? I haven't seen that in literature


%These are other model inputs that should have modeled uncertainty
%but which are not inputs to any modeled experiments
%{
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THE INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------
% ---------------------   disturbance parameters  -----------------
% --------------------------------------------------------------
input.Ru=3000;      % Ru:  upper operational wheel speed  [RPM]
input.Us=1.8;       % Us:  static wheel imbalance         [gcm]    (%0.7160 test)
input.Ud=60;        % Ud:  dynamic wheel imbalance        [gcm^2]  (%29.536 test)
input.fc=30;        % fc:  cryocooler drive frequency     [Hz]
input.Qc=0.005;     % Qc:  cryocooler attenuation factor  [-] 
input.Tst=20;	     % Tst: Star Tracker update rate       [sec]
input.Srg=3e-14;    % Srg: rate gyro noise intensity      [(rad/s)^2/Hz] 
input.Sst=2;        % Sst: Star Tracker one sigma         [arcsec]
input.Tgs=0.040;     % Tgs: Guider integration time         [sec]
% --------------------------------------------------------------
% ---------------------   plant parameters  --------------------
% --------------------------------------------------------------
input.zeta=0.005;         % zeta: Global modal damping ratio    [-]
% --------------------------------------------------------------
% ---------------------   optics parameters  -------------------
% --------------------------------------------------------------
input.lambda=1e-6;	     % lambda: center optical wavelength   [m]
input.Ro=0.98;		        % Ro:  optical surface transmissivity [-]
input.QE =0.80;		     % QE:  CCD quantum efficiency    [-]
input.Mgs=15;             % Mgs: magnitude of faint guide star  [mag]
% --------------------------------------------------------------
% ---------------------   controls parameters  -----------------
% --------------------------------------------------------------
input.fca=0.01;           % fca: Target ACS control bandwidth    [Hz]
input.Kc=0.0;             % Kc:  FSM-to-ACS cross coupling gain  [0-1]
input.Kcf=2000;    % Kcf: FSM controller gain             [-]
% -------------------------------------------------------------

%Related to number of arrays for optical simulation?
input.nray= 134;

% --------------------------------------------------------------
% -------------------------   wheel stuff  ---------------------
% --------------------------------------------------------------
input.caxi = [2.95000000000000e-08 1.41600000000000e-08 7.14000000000000e-09 .36000000000000e-09 2.49100000000000e-08 1.08700000000000e-08 3.01000000000000e-08 1.36200000000000e-08 7.05000000000000e-09];
input.crad = [7.851999999999999e-08 2.984000000000000e-08 7.770000000000000e-09 1.036000000000000e-08 1.280000000000000e-08 8.490000000000000e-09 9.850000000000000e-09 1.849000000000000e-08 6.180000000000000e-09 4.020000000000000e-09];
input.ctor = [3.239000000000000e-08 7.439999999999999e-09 3.090000000000000e-09 4.719999999999999e-09 1.210000000000000e-09 4.120000000000000e-09 9.239999999999999e-09 1.027000000000000e-08 6.300000000000000e-09 8.510000000000000e-09];
input.haxi = [1 2 2.900000000000000 3.880000000000000 4.290000000000000 4.850000000000000 5.400000000000000 5.810000000000000 6.160000000000000];
input.hrad = [1 2 3 4 4.430000000000000 5 5.380000000000000 5.600000000000000 6 6.390000000000000];
input.htor = [1 2 2.900000000000000 3.880000000000000 4 4.430000000000000 5.200000000000000 5.400000000000000 5.600000000000000 5.810000000000000];
input.zeta1=0.2;     % empirical number
input.a=0.2;         % HPF corner fraction of top wheel speed
% rotation matrices from individual wheel frames to S/C frame
% coordinate locations of wheels
  input.wheel_locs=[...
  79.00000000000000  -0.47413000000000  -0.81375000000000   2.04467000000000
  80.00000000000000  -0.47413000000000  -0.96331000000000   1.72393000000000
  81.00000000000000  -0.47413000000000  -1.28405000000000   1.87349000000000
  82.00000000000000  -0.47413000000000  -1.13448000000000   2.19423000000000
  83.00000000000000  -0.56728000000000  -1.04890000000000   1.95908000000000
];

% --------------------------------------------------------------
% -------------------------   cryo stuff  ----------------------
% --------------------------------------------------------------
input.n=3;       % number of harmonics to include (including fundamental)
input.h=[1.0 2.0 3.0 4.0 5.0 6.0]; %TODO what is this?
input.C=[42 0.95 4.1 2.75 0.9 1.2;
   0.2 0.09 0.25 1.0 5.0 0.4];  %[N]
 % values from Cassini
%input.ACS_GYRO_nw=7e-19;	% gyro rate random walk (rad^2/s^3) or (rad^2/s^4)/Hz
%input.ACS_GYRO_nf=0;	   % gyro rate flicker noise (rad^2/s^2) 
 
% --------------------------------------------------------------
% -----------------------   gs optics stuff  -------------------
% --------------------------------------------------------------
input.Nsurf = 10;	         			% number of optical surfaces before detector
input.D = 2.8;		            		% aperture diameter (m) NEXUS
% lambda = 1e-6;	      	   % center of wavelength band (m)
input.BP = 0.4;		         		% bandpass (microns)
input.PH = 1.0e10;		      		% photons/m^2/micron/sec
input.R0 = 4*60;						   % detector readout noise, 60 electrons per pixel, 4 pixels

% --------------------------------------------------------------
% --------------------------   fem stuff  ----------------------
% --------------------------------------------------------------
input.m_SM=2.49;          % m_SM: mass of secondary mirror (SM) [kg]   
input.m_SMhub=4.4;  %SM hub mass
input.I_SMhubt=0.25200E+00; %SM hub inertia transverse
input.I_SMhuba=0.45900E+00; %SM hub inertia axial
input.m_RW=0.10600E+02; % RW mass (ITHACO E-wheel)
input.K_yPM=0.77400E+06; %PM outer bipod actuator transverse stiffness
input.I_xRWA=0.40187E+00;  % Ix moment of inertia of RWA chassis
input.I_yRWA=0.22445E+00;  % Iy moment of inertia of RWA chassis
input.I_zRWA=input.I_yRWA; % Iz moment of inertia of RWA chassis
input.I_RWt=0.83595E-01; % Transverse Inertia of ITHACO E reaction wheel
input.I_RWa=0.14339E+00; % Axial Inertia of ITHACO E reaction wheel
input.m_ISO=0.15E+01;  % mass of RWA isolator strut actuator
input.I_ISOa=0.11720E-03;  % RWA isolator strut axial inertia
input.I_ISOt=0.39300E-01;  % RWA isolator strut transverse inertia
input.K_yISO=0.14600E+04;  % Hexapod (i.e. secondary mirror support structure) isolator strut axial stiffness
input.K_xISO=0.14000E+12;  % Hexapod (i.e. secondary mirror support structure) isolator strut transverse stiffness
input.m_RWAchx=0.25000E+01; % mass RWA pyramid chassis
input.I_bus=0.85080E+02;    % NEXUS spacecraft bus principal inertia
input.m_bus=0.33200E+03;    % NEXUS spacecraft bus mass
input.m_prop=0.37000E+02;   % propulsion system mass (incl. propellant)
input.I_propt=0.51100E+01;  % propulsion system transverse inertia 
input.I_propa=0.74000E+00;  % propulsion system axial inertia
input.m_instr=0.50000E+02;  % instrument mass (incl. CCD detector)
input.I_i1=0.49200E+01;  % Instrument moment of inertia 1-axis
input.I_i2=0.75420E+01;  % Instrument moment of inertia 2-axis
input.I_i3=0.41280E+01;  % Instrument moment of inertia 3-axis
% SM supporting structure (spider): Tangential Bipod design
input.A_sptop=0.14040E-02;  % Cross sectional area of SM spider (top)
input.D_sp=0.060;           % SM spider (bottom) diameter
input.t_sp=0.003;         % SM spider (bottom) wall thickness  [m]
input.A_spbot=(pi/4)*(input.D_sp^2-(input.D_sp-2*input.t_sp)^2); % Cross sectional area of SM spider (bottom)
input.J_sp=(pi/32)*(input.D_sp^4-(input.D_sp-2*input.t_sp)^4);   % SM Spider torsional moment of inertia
input.I_sp=(pi/64)*(input.D_sp^4-(input.D_sp-2*input.t_sp)^4);   % SM Spider bending moment of inertia
input.I_ss=0.78350E-08;                        % Sunshield out-of-plane bending moment of inertia,
input.K_rad1=0.50000E+6;   % rotational stiffness of SM and PM actuators in 1 direction
input.K_rad2=0.30000E+6;   % rotational stiffness of SM and PM actuators in 2 direction 
input.K_rISO=3000;  % K_rISO: rot stiff RWA iso struts [Nm/rad] 
input.K_aISO=(5/3)*input.K_rISO; % axial rot RWA iso struts 
input.K_act1=0.20000E+11;  % actuator axial stiffness 
input.K_act2=0.14000E+12;  % actuator radial stiffness +12
input.I_iso=1.00000E-5;    % RW wheel isolator bending m.o.I. [m^4]
input.K_zpet=0.9000E+08;  % Deployable petal hinge stiffness [N/m] (nom:0.90000E+08)
%Note: Inertias in units of kgm^2
%Note: Stiffnesses in units of N/m  or Nm/rad
%Note: Masses in units of kg

% --------------------------------------------------------------
% ------------------------   optics stuff  ---------------------
% --------------------------------------------------------------
%Not pulling this out yet, too big

% --------------------------------------------------------------
% --------------------------   acs stuff  ----------------------
% --------------------------------------------------------------
input.mass=[752.094612601836,1.03216046820620e-15,2.23432383705813e-15,3.23726195090639e-17,1.47160061914065e-13,4.28546087505310e-14; ...
5.79397640976254e-16,752.094612601836,-4.27435864480685e-14,-1.04638520070921e-13,2.43945523111095e-15,-6.75674793892966e-15; ...
2.50494069931051e-15,-4.03219124756049e-14,752.094612601836,-1.60062087509864e-13,4.98039109952941e-15,-1.69109226052856e-15; ...
-7.47073394395928e-16,-1.47049039611602e-13,-1.31672450720544e-13,1720.96055133990,-10.5270557086614,13.3850421894362; ...
1.47160061914065e-13,2.62130691004106e-15,4.98039109952941e-15,-10.5270557086614,1352.56559956577,228.801871490305; ...
4.28546087505310e-14,-6.75674793892966e-15,-2.65496910573658e-15,13.3850421894362,228.801871490305,1508.50841532254];

% --------------------------------------------------------------
% --------------------------   fsm stuff  ----------------------
% --------------------------------------------------------------
input.FgsNom = 30;			      % nominal FSM sample rate (Hz)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL THE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Q1, Q2] = nexus_lyapunov(input)