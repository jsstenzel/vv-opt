% nexus_parfinals.m
% conversion factors
r2a = 3600*180/pi;	    % radians to arc-sec
a2r = 1/r2a;			% arc-sec to radians
nray= 134;

p_symb=str2mat('Ru','Us','Ud','fc','Qc','Tst','Srg','Sst','Tgs',...
   'm_SM','K_yPM','K_rISO','m_bus','K_zpet','t_sp','I_ss','I_propt','zeta',...
   'lambda','Ro','QE','Mgs','fca','Kc','Kcf');

p_strings=str2mat('upper operational wheel speed','static wheel imbalance',...
   'dynamic wheel imbalance','cryocooler drive frequency','cryocooler attenuation factor',...
   'star tracker update rate','rate gyro noise intensity','star tracker one sigma',...
   'guider integration time','secondary mirror mass','PM bipod act trans stiffness',...
   'RWA isolator joint stiffness','spacecraft bus mass','deploy petal hinge stiffness',...
   'SM spider wall thickness','sunshield bending m.o. inertia','propulsion sys trv inertia',...
   'proportional damping ratio','center optical wavelength','optical surface reflectivity',...
   'CCD quantum efficiency','guide star magnitude','target ACS bandwidth','FSM-ACS cross coupling gain',...
   'FSM controller gain');
p_units=str2mat('[RPM]','[gcm]','[gcm^2]','[Hz]','[-]','[sec]','[(rad/sec)^2/Hz]',...
   '[asec]','[sec]','[kg]','[N/m]','[Nm/rad]','[kg]','[N/m]','[m]','[m^4]','[kgm^2]',...
   '[-]','[m]','[-]','[-]','[mag]','[Hz]','[-]','[-]');

% --------------------------------------------------------------
% ---------------------   disturbance parameters  -----------------
% --------------------------------------------------------------
Ru=3000;      % Ru:  upper operational wheel speed  [RPM]
Us=1.8;       % Us:  static wheel imbalance         [gcm]    (%0.7160 test)
Ud=60;        % Ud:  dynamic wheel imbalance        [gcm^2]  (%29.536 test)
fc=30;        % fc:  cryocooler drive frequency     [Hz]
Qc=0.005;     % Qc:  cryocooler attenuation factor  [-] 
Tst=20;	     % Tst: Star Tracker update rate       [sec]
Srg=3e-14;    % Srg: rate gyro noise intensity      [(rad/s)^2/Hz] 
Sst=2;        % Sst: Star Tracker one sigma         [arcsec]
Tgs=0.040;     % Tgs: Guider integration time         [sec]
% --------------------------------------------------------------
% ---------------------   plant parameters  --------------------
% --------------------------------------------------------------
m_SM=2.49;          % m_SM: mass of secondary mirror (SM) [kg]
K_yPM=0.774E+06;    % K_yPM: PM bipod act trans stiffness [N/m]
K_rISO=3000;  % K_rISO: rot stiff RWA iso struts [Nm/rad]
m_bus=0.332E+03;    % m_bus:  NEXUS spacecraft bus mass   [kg]
K_zpet=0.9E+008;  % K_zpet: Deployable petal hinge stiffness [N/m]
t_sp=0.003;         % SM spider (bottom) wall thickness  [m]
I_ss=0.7835E-08;    % I_ss: Sunshield bending mof inertia [m^4]     
I_propt=0.511E+01;  % I_propt: propulsion sys trv inertia [kgm^2]
zeta=0.005;         % zeta: Global modal damping ratio    [-]
% --------------------------------------------------------------
% ---------------------   optics parameters  -------------------
% --------------------------------------------------------------
lambda=1e-6;	     % lambda: center optical wavelength   [m]
Ro=0.98;		        % Ro:  optical surface transmissivity [-]
QE =0.80;		     % QE:  CCD quantum efficiency    [-]
Mgs=15;             % Mgs: magnitude of faint guide star  [mag]
% --------------------------------------------------------------
% ---------------------   controls parameters  -----------------
% --------------------------------------------------------------
fca=0.01;           % fca: Target ACS control bandwidth    [Hz]
Kc=0.0;             % Kc:  FSM-to-ACS cross coupling gain  [0-1]
Kcf=2000;    % Kcf: FSM controller gain             [-]
% -------------------------------------------------------------
for ind=1:25
   eval(['p_final(' num2str(ind) ')=eval([p_symb(' num2str(ind) ',:)])'';' ])
end
% ------------------------- bounds -----------------------------
p_LB=[ 2400 0.1 1.0 20 0.001 2 1E-14 1 0.004 1.0 0.2E+06 30 280 .1E+08 ...
      .001 .5E-08 2 .001 .5E-6 0.8 0.5 8 0.001 0 0]';
p_UB=[ 3850 2.7 90 60 0.0250 40 1E-13 4 0.4  7.7 1.5E+06 5E3 340 1.8E+09 ...
      .005 2.0E-08 10 0.01 2E-6 .99 .99 20 .2 1 1000000]';
p_bounds=[p_LB p_UB];
p_final=p_final';
