% nexus_fsm.m
% Synthesis a Lead-Lag FSm controller based on GSFC design
% variable parameters
FgsDes = 1/Tgs;	     % desired FSM sample rate [Hz] % guider
% Kcf= -2000;            % FSM controller gain [-]

% other fixed parameters
FgsNom = 30;			      % nominal FSM sample rate (Hz)
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

 %  varargout{1}=Tg;
 %  varargout{2}=Kcf;
 %  varargout{3}=wlz;
 %  varargout{4}=wlp;
 %  varargout{5}=wlpp1;
 %  varargout{6}=wlpp2;
 %  varargout{7}=FsmIntGain;
   
   %%% FSM mirror plant
   % eliminated due to numerical ill-conditioning
	% wp = 2*pi*500;			% natural frequency; Model 3A Ball FSM resonance
	% zp = 0.02;				% damping
   % Apf=[0 1; -wp^2 -2*zp*wp]; Bpf=[0 1]';
   % Cpf=[wp^2 0]; Dpf=[0];
   % causes ill-conditioning replace by unity gain (feedthrough)
   % Apf=[0]; Bpf=[0]; Cpf=[0]; Dpf=[1];
   % 2 channels in parallel, assume no coupling in FSM plant
   % Apf=[Apf zeros(1,1); zeros(1,1) Apf];
   % Bpf=[Bpf zeros(1,1); zeros(1,1) Bpf];
   % Cpf=[Cpf zeros(1,1); zeros(1,1) Cpf];
   % Dpf=[Dpf 0; 0 Dpf];
   
