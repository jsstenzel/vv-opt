% nexus_gs.m
% Approximation of (FGS) guide star sensor noise. The NEA 
% depends on the integration time T_int and is approximated
% as a second order LPF. 
% Author: Olivier de Weck
% (c) Massachusetts Institute of Technology
% Date: 26 June 2000
% Revision: 22 August 2000
% Modified for NEXUS: 8 June 2001

%----------------------------------------------------------------
% Ro=0.90;		        % Ro:  optical surface transmissivity [-]
% QE = 0.80;		     % QE:  InSb CCD quantum efficiency    [-]
% Mgs=16.0;            % Mgs: magnitude of faint guide star  [mag]
% Tgs=0.033;           % Tgs: Guider integration time        [sec]
%----------------------------------------------------------------

Nsurf = 10;	         			% number of optical surfaces before detector
D = 2.8;		            		% aperture diameter (m) NEXUS
% lambda = 1e-6;	      	   % center of wavelength band (m)
BP = 0.4;		         		% bandpass (microns)
PH = 1.0e10;		      		% photons/m^2/micron/sec
R0 = 4*60;						   % detector readout noise, 60 electrons per pixel, 4 pixels
r2a = 3600 * 180 / pi;			% radians-to-arcseconds conversion

%----------------------------------------------------------------
% overall transmissibility 
Reff = Ro^Nsurf;
% area of primary mirror
Area = pi * D^2 / 4;
% slope of centroiding transfer function
k = (16/3/pi) * (D/lambda) / r2a;
alpha=Reff * QE * Area * BP * 10.^(-0.4.*Mgs) * PH;
% total number of detected photons
N = alpha * Tgs;
% raw NEA
NEA = sqrt(1+R0./N)/k./sqrt(N); % compare this with rms output from SS
%----------------------------------------------------------------
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
%----------------------------------------------------------------
if verification&0
   Sqc=lyap(Adg,Bdg*Bdg');
   Szc=Cdg*Sqc*Cdg';
   rmsdg=sqrt(diag(Szc));
   sysdg=ss(Adg,Bdg,Cdg,Ddg);
end
%----------------------------------------------------------------
% sensitivities
% dA_dp=[ 0 0; -2*w_fgs -2*zeta_fgs]*(-2*pi/T_int^2);
% dA_dp=[dA_dp zeros(2,2); zeros(2,2) dA_dp];
% dB_dp=zeros(size(B));
% delta=(sqrt(zeta_fgs)/k)*((1/2/pi/alpha)+(R0/2/pi/alpha^2/T_int))^(-1/2)*(-R0/2/pi/alpha^2/T_int^2);
% dC_dp=[w_fgs^2 0]*delta+[2*k_NEA*w_fgs 0]*(-2*pi/T_int^2);
% dC_dp=[dC_dp [0 0]; [0 0] dC_dp];
% dD_dp=zeros(size(D));
% dsys_d_3=ss(dA_dp,dB_dp,dC_dp,dD_dp);

%----------------------------------------------------------------
if plotflag&0
   f=logspace(-1,3,1000);
   [Gdg,tmp]=bode(sysdg(1,1),2*pi*f); Gdg=squeeze(Gdg);
   S=squeeze(Gdg.^2);
   figure
   loglog(f,S);
   xlabel('Frequency [Hz]')
   ylabel('PSD FGS Noise [asec]')
   title ('LPF Approximation of FGS Noise')
end
%----------------------------------------------------------------
if verification&0
   f=logspace(-1,3,2000);
   [mag,phs]=bode(sysdg(1,1),2*pi*f);
   mag=squeeze(mag);
   S=zeros(2,2,length(f)); varpdg=zeros(2,2);
   S(1,1,:)=mag.^2; S(2,2,:)=mag.^2;
   for ind1=1:2
      for ind2=1:2
       varpdg(ind1,ind2)=abs(2*(1/2/pi)*trapz(f*2*pi,squeeze(S(ind1,ind2,:))));
      end
   end
   fdg=f;
   Sdg=S;
   rmsSdg=sqrt(varpdg);
end
%----------------------------------------------------------------

