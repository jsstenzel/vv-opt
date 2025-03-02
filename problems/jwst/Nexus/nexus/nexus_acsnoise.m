  %  nexus_acsnoise.m
  %  The function generates a ss (and PSD) for ACS sensor noise (star tracker and
  %  gyro noise).  This is additive noise on the attitude angle
  %  outputs.  It is an EXTREMELY simplified model.
  %  One defines a frequency vector (FreqBase) and sets the hardware parameters near
  %  the top of the file and the output of interest is:
  %  ACS_PSD_noise = ACS sensor noise PSD (rad^2/Hz)
  %  ACS_PSD_noise_rms = rms of the ACS_PSD_noise PSD (rad)
  %  Depending on how the ACS loop is designed one can play
  %  around with the hardware parameters to get about 1 arcsec RMS (the present
  %  baseline number) on the attitude angles in the closed loop.
  %  Original file obtained from Bob Grogan (JPL): 13 Oct 2000
  %  Modified for use with the NEXUS Model: dWo, 8 June 2001
  
  %--------------------------------------------------------------------------------------
  %  variable hardware parameters
     ACS_ST_period=Tst; % Star Tracker update rate (sec) nom: Cassini 100 sec
     ACS_GYRO_nr=Srg;       	% gyro rate noise (rad^2/s) or (rad/s)^2/Hz nom: 3e-14ST
     ACS_ST_sigma=Sst;             % Star Tracker one sigma resolution (arcsec)
     Sstrad=Sst/60/60*pi/180;	%convert to radians
     ACS_ST_sigma=ACS_ST_sigma/60/60*pi/180;	%convert to radians

     % values from Cassini
     ACS_GYRO_nw=7e-19;	% gyro rate random walk (rad^2/s^3) or (rad^2/s^4)/Hz
     ACS_GYRO_nf=0;	   % gyro rate flicker noise (rad^2/s^2) 
     %  frequency base
     FreqBase=logspace(-5,3,1000);
     rad2arcsec=(3600*180/pi);
     f1=find(FreqBase <= 2/ACS_ST_period);
  %--------------------------------------------------------------------------------------
     if verification&1
     % angle
       S1=(ACS_ST_sigma^2)*ones(size(FreqBase(f1)));
       f2=find(FreqBase > 2/ACS_ST_period);
       S2=ACS_GYRO_nw./(FreqBase(f2).^4)+ACS_GYRO_nf./(FreqBase(f2).^3)+ACS_GYRO_nr./(FreqBase(f2).^2);
       ACS_PSD_noise=[S1 S2]; % units: rad^2/Hz
       [ACS_PSD_noise_rms]=intrms(2*ACS_PSD_noise',FreqBase');  % units: rad RMS
       ACS_PSD_noise_rms_arcsec = rad2arcsec*ACS_PSD_noise_rms; % units: arcsec RMS 
       if diagnostics
          disp(['RMS (arcsec) = ',num2str(ACS_PSD_noise_rms_arcsec)]);
          end
  
       if plotflag
          figure; loglog(FreqBase(f1),rad2arcsec*sqrt(S1),'b',...
                 FreqBase(f2),rad2arcsec*sqrt(ACS_GYRO_nw./(FreqBase(f2).^4)),'r',...
                 FreqBase(f2),rad2arcsec*sqrt(ACS_GYRO_nr./(FreqBase(f2).^2)),'g',...
                 FreqBase,rad2arcsec*sqrt(ACS_PSD_noise),'m--'); 
              grid; ylabel('PSD (arcsec/sqrt(Hz))'); 
              xlabel('Frequency (Hz)'); 
              title('ACS Sensor Noise Model')
         legend('Star Tracker',...
         'Gyro Rate Random Walk',...
         'Gyro Rate Noise',...
         'Composite');
         h=gcf; line_handle=findobj(h,'Type','line'); 
         set(line_handle,'LineWidth',2);
      end
     end
%--------------------------------------------------------------------------------------
Kn=sqrt((1/2)*(Sstrad^2+Srg./(2/Tst).^2));  % ACS noise gain

% state space filter approximation to composite
% ACS sensor noise model.
wn=2*pi*(2/Tst);
Ads=[-wn]; 
Bds=[1]';
Cds=Kn*[wn];
Dds=[0]; 
SYSds=ss(Ads,Bds,Cds,Dds);
%--------------------------------------------------------------------------------------
if verification&1
   [magds,phsds]=bode(SYSds,FreqBase*2*pi); 
   magds=abs(squeeze(magds));
   Snn=magds.^2;
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
end
%---------------------------------------
% State space for 3 axes
%---------------------------------------
Systmp=parallel(SYSds,SYSds,[],[],[],[]);
SYSds=parallel(Systmp,SYSds,[],[],[],[]);
[Ads,Bds,Cds,Dds]=ssdata(SYSds);
% Note state space outputs of SYSds are in rad !





