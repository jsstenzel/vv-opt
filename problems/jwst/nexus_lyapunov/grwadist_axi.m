function [Fa,f] = grwadist_axi(RPM,model,scale,haxi,caxi,T)
%       grwadist_axi	calculates the axial force disturbance of a RW
%
%	[Fa,f] = grwadist_axi(RPM,model,scale,haxi,caxi,T
%
%	Fa	axial forces
%	f	axial force disturbance component frequencies
%			(see below)
%	RPM  wheel speed (RPM)
%       model   Reaction wheel model that was chosen
%       scale   scale factor to scale RW disturbance level. This factor
%               can be used to account for variability between wheels
%               and/or to conduct trade studies. Default scale=1.0
%       haxi    wheel axial harmonics
%       caxi    wheel axial amplitude coeffcients
%       T     	 time base row vector for the disturbance waveforms
%
%	The wheel disturbance is calculated based on the empirical 
%	disturbance models of the above wheel models.   
%	T is the time base used to calculate the disturbances waveforms.
%	It must be a vector.  The angular velocity of the wheel, RPM, can either be
%	a single number, corresponding to a constant velocity, or it can be a time 
%	history of the velocity, having the same dimension as T.  The axial
%	force time history is returned in Fa.
%	If T is not specified, then the amplitudes (N) and
%	frequencies (Hz) of the disturbance components are returned in Fa
%	and f, respectively.
%

%  History
%  created by Jim Melody
%  22Feb93 jmelody: made public
%  10Oct93 jmelody: fixed up zeros(size()) warning
%  09Sep96 jmelody: changed print out for nargout = 0
%  12Sep96 jmelody: fixed bug for time history output
%  08Aug98 odeweck: added model variable and scale   
%  22Sep98 odeweck: added hax and Cax for E- and B-wheels
%  10Oct99 odeweck: modified version for use with RWA modeler GUI

Cax=caxi;
hax=haxi;

if length(scale)==1;
Cax=scale*Cax; % scale disturbance amplitudes for model
end



[mRPM,nRPM]=size(RPM);
w=RPM*pi/30;

if (nargin == 5)
  if ((mRPM==1)&(nRPM==1))
  %wheel speed must be constant
    amplitude=zeros(size(Cax));
    frequency=zeros(size(Cax));
    for n=1:max(size(amplitude))
      amplitude(n,1)=(RPM)^2*Cax(n);
      frequency(n,1)=hax(n)*w/2/pi;
    end		%end loop over harmonics
    if (nargout == 0)
      disp(' ');
      disp('   Ampl (N)  Freq (Hz)');
      disp([amplitude frequency]);
      disp(' ');
    elseif nargout == 2
      Fa=amplitude;
      f=frequency;
    end
    return
  else
    disp('RPM must be a single number!')
    return
  end
end

[mT,nT]=size(T);

if (mT~=1)
  if (nT==1)
    T=T';
    [mT,nT]=size(T);
  else
    'The time base is not a vector'
    return
  end
end 
  
if ((mRPM==1)&(nRPM==1))
  Fa=zeros(size(T));
  f=zeros(size(T));
  for n=1:max(size(Cax))
    Fa=(RPM)^2*(Cax(n))*sin(hax(n)*w*T)+Fa;
  end		%end loop over harmonics
  return
end

if (mRPM~=1)
  if (nRPM==1)
    RPM=RPM';
    [mRPM,nRPM]=size(RPM);
    w=w';
  else
    error('The wheel speed history is not a vector');
  end
end

if (nRPM~=nT)
  error('The vectors do not have the same length');
end

Fa=zeros(size(T));
f=zeros(size(T));
for n=1:nT
  for m=1:max(size(Cax))	%loop over harmonics
    Fa(1,n)=(RPM(1,n))^2*(Cax(m))*sin(hax(m)*w(1,n)*T(1,n))+Fa(1,n);
  end		%end loop over harmonics
end
