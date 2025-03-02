function [Tx,Ty] = grwadist_tor(RPM,model,scale,htor,ctor,T)
%       rwadist_tor	calculates the radial torque disturbances of an RWA.
%
%	[Tx,Ty] = grwadist_tor(RPM,model,T,scale,htor,ctor,T)
%
%	Tx	radial torques about the x-axis
%	Ty	radial torques about the y-axis
%	RPM	wheel speed (RPM)
%       model   Reaction wheel model that was chosen
%       scale   scale factor to scale RW disturbance level. This factor
%               can be used to account for variability between wheels
%               and/or to conduct trade studies. Default scale=1.0
%       htor    torque wheel harmonics
%       ctor    torque amplitude coefficients
%
%	     T	time base row vector for the disturbance waveforms
%
%
%	The wheel disturbances are calculated based on the empirical 
%	disturbance model of the above wheel models. 
%	T is the time base used to calculate the disturbances waveforms.
%	It must be a vector.  The angular velocity of the wheel, RPM, can either be
%	a single number, corresponding to a constant velocity, or it can be a time 
%	history of the velocity, having the same dimension as T.
%	If T is not specified, then the amplitudes (N) and
%	frequencies (Hz) of the disturbances are returned in Tx and Ty,
%	respectively.
%

%  History
%  created by Jim Melody 
%  22Feb93 jmelody: made public
%   6Aug94 jmelody: fixed up zeros(size()) warning
%  09Sep96 jmelody: changed print out for nargout = 0
%  08Aug98 odeweck: added model variable and scale 
%  22Sep98 odeweck: added htor and Ctor for E- and B-wheels
%  10Oct99 odeweck: modified version for use with RWA modeler GUI

if length(scale)==1;
   Ctor=scale*ctor; % scale disturbance amplitudes for model
else
   Ctor=ctor;
   Ctor(1)=scale(2)*Ctor(1); % scale disturbance amplitudes for Ud fundamental only
end

[mRPM,nRPM]=size(RPM);
w=RPM*pi/30;

if nargin == 5
  if ((mRPM==1)&(nRPM==1))
  %wheel speed must be constant
    amplitude=zeros(size(Ctor));
    amplitude=zeros(size(Ctor));
    for n=1:max(size(amplitude))
      amplitude(n,1)=(RPM)^2*Ctor(n);
      frequency(n,1)=htor(n)*w/2/pi;
    end		%end loop over harmonics
    Tx=amplitude;
    Ty=frequency;
    if (nargout == 0)
      disp(' ');
      disp('   Ampl (N)  Freq (Hz)');
      disp([amplitude frequency]);
      disp(' ');
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
  Tx=zeros(size(T));
  Ty=zeros(size(T));
  for n=1:max(size(Ctor))
    Tx=(RPM)^2*(Ctor(n))*sin(htor(n)*w*T)+Tx;
    Ty=(RPM)^2*(Ctor(n))*cos(htor(n)*w*T)+Ty;
  end		%end loop over harmonics
return
end

if (mRPM~=1)
  if (nRPM==1)
    RPM=RPM';
    [mRPM,nRPM]=size(RPM);
    w=w';
  else
    'The wheel speed history is not a vector'
    return
  end
end

if (nRPM~=nT)
	'The vectors do not have the same length'
	return
end

Tx=zeros(size(T));
Ty=zeros(size(T));
for n=1:nT
  for m=1:max(size(Ctor))	%loop over harmonics
    Tx(1,n)=(RPM(1,n))^2*(Ctor(m))*sin(htor(m)*w(1,n)*T(1,n))...
		+Tx(1,n);
    Ty(1,n)=(RPM(1,n))^2*(Ctor(m))*cos(htor(m)*w(1,n)*T(1,n))...
		+Ty(1,n);
  end		%end loop over harmonics
end
