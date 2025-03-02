% nexus_rwa.m
% computes the Nexus RW disturbances as a function
% of static and dynamic imbalance and upper wheel
% speed limit;
% 7 June 2001, dWo
% This code only valid with E-wheel model
% Model number: TW-50E300

% variable parameters
% Ru=3000;            % upper operational wheel speed [RPM]
% Us=0.7160;          %[gcm]
% Ud=29.536;          %[gcm^2];

%----------------------------------------------------------
dR=Ru/2;              % +/- wheel speed deviation
scale(1)=Us/0.7160;  % scale factor for Us
scale(2)=Ud/29.536;  % scale factor for Ud
f=logspace(-2,3,500);
type=0;             % uniform wheel speed pdf
model=[];
%----------------------------------------------------------
load RW_data_Ewheel
Rl=0;
Crad1=((2*pi)^2/(1000*60^2*100))*Us;
Caxi1=caxi(1);  %note axial force not affected by imbalances
Ctor1=((2*pi)^2/(1000*60^2*100^2))*Ud;
Sradm=pi*Crad1^2/(2*(Ru-Rl)*(pi/30)^5)*(2*pi*(Ru/60))^4; %PSD value
Saxim=pi*Caxi1^2/(2*(Ru-Rl)*(pi/30)^5)*(2*pi*(Ru/60))^4; %PSD value
Storm=pi*Ctor1^2/(2*(Ru-Rl)*(pi/30)^5)*(2*pi*(Ru/60))^4; %PSD value
zeta1=0.2;     % empirical number
f1max=Ru/60;
f2r=max(hrad)*f1max;
f2a=max(haxi)*f1max;
f2t=max(htor)*f1max;
krad1=4*pi*zeta1*f1max*sqrt(Sradm);
kaxi1=4*pi*zeta1*f1max*sqrt(Saxim);
ktor1=4*pi*zeta1*f1max*sqrt(Storm);
a=0.2;         % HPF corner fraction of top wheel speed
%----------------------------------------------------------
if verification
% frequency domain approximation
numhpf=[1 0];
numfr=[krad1 0]; numfa=[kaxi1 0]; numft=[ktor1 0];
numlpfr=[2*pi*f2r]; numlpfa=[2*pi*f2a]; numlpft=[2*pi*f2t];
denhpf=[1 a*2*pi*f1max];
denf=[1 2*zeta1*2*pi*f1max 4*pi^2*f1max^2];
denlpfr=[1 2*pi*f2r]; denlpfa=[1 2*pi*f2a]; denlpft=[1 2*pi*f2t];
numar=conv(conv(numhpf,numfr),numlpfr);
numaa=conv(conv(numhpf,numfa),numlpfa);
numat=conv(conv(numhpf,numft),numlpft);
denar=conv(conv(denhpf,denf),denlpfr);
denaa=conv(conv(denhpf,denf),denlpfa);
denat=conv(conv(denhpf,denf),denlpft);
[magr,phsr]=bode(numar,denar,f*2*pi);
[maga,phsa]=bode(numaa,denaa,f*2*pi);
[magt,phst]=bode(numat,denat,f*2*pi);
% compute component PSD's for fundamental only
S_rad1=grwabroad_rad(f,Ru/2,dR,type,model,scale,hrad,crad);
S_axi1=grwabroad_axi(f,Ru/2,dR,type,model,scale,haxi,caxi);
S_tor1=grwabroad_tor(f,Ru/2,dR,type,model,scale,htor,ctor);
    if plotflag
       figure
       loglog(f,magr.^2,'b:',f,maga.^2,'g:',f,magt.^2,'r:')
       hold on
       loglog(f,S_rad1,f,S_axi1,f,S_tor1)
       xlabel('Frequency [Hz]')
       axis([1 1000  1e-10 1e-1])
       ylabel('RW Disturbance PSD')
       legend('Radial Force [N^2/Hz]','Axial Force [N^2/Hz]','Radial Torque [Nm^2/Hz]',2)
       title('Reaction Wheel Disturbance PSDs in wheel frame')
    end
    
    %considered using notch filter to get "dip" after first harmonic
    %dn=0.01; % notch depth
    %numnot=[1 2*zeta1*dn*2*pi*f2max 4*pi^2*f2max^2];
    %dennot=[1 2*2*pi*f2max 4*pi^2*f2max^2];
 end
 %----------------------------------------------------------
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

if verification&plotflag
sysdw=ss(Adw,Bdw,Cdw,Ddw);
[magdw,phsdw]=bode(sysdw,f*2*pi); magdw=squeeze(magdw);
loglog(f,magdw.^2,'--')
end
%----------------------------------------------------------
% rotation matrices from individual wheel frames to S/C frame
% coordinate locations of wheels
  wheel_locs=[...
  79.00000000000000  -0.47413000000000  -0.81375000000000   2.04467000000000
  80.00000000000000  -0.47413000000000  -0.96331000000000   1.72393000000000
  81.00000000000000  -0.47413000000000  -1.28405000000000   1.87349000000000
  82.00000000000000  -0.47413000000000  -1.13448000000000   2.19423000000000
  83.00000000000000  -0.56728000000000  -1.04890000000000   1.95908000000000
];
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
%----------------------------------------------------------


