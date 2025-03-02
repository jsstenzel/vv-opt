% nexus_cryo.m
% Approximation of a cryocooler disturbance with
% a state space system, where each harmonic of force
% and torque is represented by a lightly damped
% 2nd-order system. Based on BAe80K experimental
% data from GSFC.
% fc-cryocooler drive frequency [Hz]
% Qc-cryocooler attenuation factor (relative to raw test data)
% n-number of harmonics to be included in model (<= 6)
% options(1)=1 - compute time domain signals and rms
% options(2)=1 - conduct PSD analysis of cryocooler
% Outputs:
% sysdc  -  state space system with 1 white noise input
%           and three outputs: Fx, Fy, Fz, 
% dsysdc -  structure containing the matrix derivatives
%           with repect to Qc
% rmsdc  -  RMS values of force for all 3 channels
% dWo, 26 June 2000
% Modified for NEXUS: 8 June 2001

%-----------------------------------------------------------------
% fc=40;     % cryocooler drive frequency [Hz]
% Qc=0.01;  % cryocooler attenuation factor [-] 
n=3;       % number of harmonics to include (including fundamental)
%
%-----------------------------------------------------------------
wc=2*pi*fc; % drive frequency
h=[1.0 2.0 3.0 4.0 5.0 6.0];
C=[42 0.95 4.1 2.75 0.9 1.2;
   0.2 0.09 0.25 1.0 5.0 0.4];  %[N]
Q=Qc*ones(1,n);  
%-----------------------------------------------------------------
if verification&0
   dt=1/(20*wc/2/pi);
   t=[0:dt:4];
   phi=2*pi*rand(1); % random initial phase
   Ftmp=[]; F=[];
     for ind1=1:2
      for ind2=1:n,          
  	       Ftmp(ind2,:)=Q(ind2)*C(ind1,ind2)*sin(wc*h(ind2)*t+phi);
      end
          F(ind1,:)=sum(Ftmp);   % sum the contribution of all harmonics       
          rmst(ind1)=std(F(ind1,:));
       end
       if diagnostics
          disp(['RMS values for time domain: ' num2str(rmst)])
          end
   if plotflag
      figure
      comp_string=str2mat('Axial Force [N]','Lateral Force [N]');
      psd_units=['[N^2/Hz]';'[N^2/Hz]'];
      for ind1=1:2
          eval(['subplot(21' num2str(ind1) ')' ]);
          eval(['plot(t,F(' num2str(ind1) ',:));' ]);
          xlabel('time [sec]'), 
          eval(['title(''' char(comp_string(ind1,:)) ''')' ]);
          grid on
      end
   end
   f=logspace(1,3,10000); fadd=[];
      for ind2=1:n
          fadd=[ fadd linspace(0.9*wc*h(ind2)/2/pi,1.1*wc*h(ind2)/2/pi,5000) ];
      end
   f=unique(f);
   S=zeros(2,length(f));
   for ind1=1:2
      for ind2=1:n
          zij(ind1,ind2)=inv(2*wc*h(ind2));
          kij(ind1,ind2)=2*zij(ind1,ind2)*wc*h(ind2)*Q(ind2)*C(ind1,ind2);
          eval(['num' num2str(ind2) '=[kij(' num2str(ind1) ',' num2str(ind2) ') 0];']);
          eval(['den' num2str(ind2) '=[1 2*zij(' num2str(ind1) ',' num2str(ind2) ')*wc*h(' num2str(ind2) ...
               ') wc^2*h(' num2str(ind2) ')^2];']);
          eval(['[mag,phs]=bode(num' num2str(ind2) ',den' num2str(ind2) ',f*2*pi);'])
          mag=squeeze(mag); Stmp=mag.^2;
          S(ind1,:)=S(ind1,:)+Stmp'; 
      end
      rmsp(ind1)=sqrt(2*trapz(f,S(ind1,:)));
   end
   if diagnostics
      disp(['RMS values for PSD representation: ' num2str(rmsp)])
      end
   if plotflag
   figure
      for ind1=1:2
          eval(['subplot(21' num2str(ind1) ')' ]);
          eval(['loglog(f,2*pi*S(' num2str(ind1) ',:));' ]);
          xlabel('Frequency [Hz]'), 
          eval([ 'ylabel(''' psd_units(ind1,:) '''), grid' ]);
          eval([ 'title(''' comp_string(ind1,:) ''')' ]);
      end
   end
end
%-----------------------------------------------------------------
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
   %if verification
   %   sysdc=ss(Adc,Bdc,Cdc,Ddc);
   %   sq=lyap(Adc,Bdc*Bdc');
   %   sz=Cdc*sq*Cdc';
   %   rmsdc=sqrt(diag(sz))';
   %   if diagnostics
   %      disp(['RMS values for Cryo SS representation: ' num2str(rmsdc)])
   %      end
   % end
%-----------------------------------------------------------------
% compute sensitivity of rms w.r.t. Qc
% expect F/sqrt(2) as result;
%dA_dp=zeros(size(a));
%dB_dp=zeros(size(b));
%dD_dp=zeros(size(d));
%for ind1=1:4
%   drms2_dp(ind1)=...
%   trace(sq*dC_dp(ind1,:)'*c(ind1,:)+sq*c(ind1,:)'*dC_dp(ind1,:));
%   drms_dp(ind1)=(1/(2*rms2(ind1)))*drms2_dp(ind1);
% end
%
%sys_d_2=sys; dsys_d_2=ss(dA_dp,dB_dp,dC_dp,dD_dp);
% finite difference verification
