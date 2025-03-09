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
%m_SM=x.m_SM;          % m_SM: mass of secondary mirror (SM) [kg]
%K_yPM=x.K_yPM;    % K_yPM: PM bipod act trans stiffness [N/m]
%K_rISO=x.K_rISO;  % K_rISO: rot stiff RWA iso struts [Nm/rad]
%m_bus=x.m_bus;    % m_bus:  NEXUS spacecraft bus mass   [kg]
%K_zpet=x.K_zpet;  % K_zpet: Deployable petal hinge stiffness [N/m]
%t_sp=x.t_sp;         % SM spider (bottom) wall thickness  [m]
%I_ss=x.I_ss;    % I_ss: Sunshield bending mof inertia [m^4]     
%I_propt=x.I_propt;  % I_propt: propulsion sys trv inertia [kgm^2]
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
xyz=x.xyz;
ci=x.ci;
cf=x.cf;
nibar=x.nibar;
propbar=x.propbar;
nicelas2=x.nicelas2;
niconm2=x.niconm2;		%concentrated masses
cmcid2=x.cmcid2;
mass=x.mass;
FgsNom=x.FgsNom;			      % nominal FSM sample rate (Hz)


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
SYSds=ss(Ads,Bds,Cds,Dds);

[magds,phsds]=bode(SYSds,FreqBase*2*pi); 
magds=abs(squeeze(magds));
Snn=magds.^2;
%TODO what does intrms do
rmsssrad=intrms(2*Snn,FreqBase);
rmsssarcsec=rad2arcsec*rmsssrad;
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

%TODO what does km2loc do
if diagnostics
disp('Transforming k and m to local cs')
end
[k,m]=km2loc(k,m,ti,tf);

%TODO what does conm2 do
if diagnostics
disp( 'Assembling k for spring elements')
end
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

%TODO what exactly does eigfem do
%we are solving for the eigenvalues and eigenvectors of this spring mass problem?
[phi,omeg]=eigfem(k,m);
   %disp(['Solved Eigenproblem in ' num2str(toc) ' [sec]'])

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
%TODO turn this into parameters, not magic numbers
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