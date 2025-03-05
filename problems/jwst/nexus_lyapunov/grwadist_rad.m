function [Fx,Fy] = grwadist_rad(RPM,model,scale,hrad,crad,T)
%  rwadist_rad	calculates the radial force disturbances of an RW.
%
%	[Fx,Fy] = grwadist_rad(RPM,model,scale,hrad,crad,T)
%
%	Fx	radial forces in the x-direction
%	Fy	radial forces in the y-direction
%	RPM	wheel speed (RPM)
%       model   Reaction wheel model that was chosen
%       scale   scale factor to scale RW disturbance level. This factor
%               can be used to account for variability between wheels
%               and/or to conduct trade studies. Default scale=1.0
%       hrad    radial wheel harmonics
%       crad    radial amplitude coefficients
%	     T	    time base row vector for the disturbance waveforms
%
%	The wheel disturbances are calculated based on the empirical 
%	disturbance model of the above wheel models.   
%	T is the time base used to calculate the disturbances waveforms.
%	It must be a vector.  The angular velocity of the wheel, RPM, can either be
%	a single number, corresponding to a constant velocity, or it can be a time 
%	history of the velocity, having the same dimension as T.
%	If T is not specified, then the amplitudes (N) and
%	frequencies (Hz) of the disturbances are returned in Fx and Fy,
%	respectively.
%

%History
%          jmelody: created
%  22Feb93 jmelody: made public
%  10Oct93 jmelody: fixed up zeros(size()) warning
%  09Sep96 jmelody: changed print out for nargout = 0
%  08Aug98 odeweck: added model variable and scale 
%  22Sep98 odeweck: added hrad and Crad for E- and B-wheels
%  10Oct99 odeweck: modified version for use with RWA modeler GUI

if length(scale)==1;
   Crad=scale*crad; % scale disturbance amplitudes for model
else
   Crad=crad;
   Crad(1)=scale(1)*Crad(1); % scale disturbance amplitudes for Us fundamental only
end


[mRPM,nRPM]=size(RPM);
w=RPM*pi/30;

if nargin == 5
  if ((mRPM==1)&(nRPM==1))
  %wheel speed must be constant
    amplitude=zeros(size(Crad));
    amplitude=zeros(size(Crad));
    for n=1:max(size(amplitude))
      amplitude(n,1)=(RPM)^2*Crad(n);
      frequency(n,1)=hrad(n)*w/2/pi;
    end		%end loop over harmonics
    Fx=amplitude;
    Fy=frequency;
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
  Fx=zeros(size(T));
  Fy=zeros(size(T));
  for n=1:max(size(Crad))
    Fx=(RPM)^2*(Crad(n))*sin(hrad(n)*w*T)+Fx;
    Fy=(RPM)^2*(Crad(n))*cos(hrad(n)*w*T)+Fy;
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

Fx=zeros(size(T));
Fy=zeros(size(T));
for n=1:nT
  for m=1:max(size(Crad))	%loop over harmonics
    Fx(1,n)=(RPM(1,n))^2*(Crad(m))*sin(hrad(m)*w(1,n)*T(1,n))...
		+Fx(1,n);
    Fy(1,n)=(RPM(1,n))^2*(Crad(m))*cos(hrad(m)*w(1,n)*T(1,n))...
		+Fy(1,n);
  end		%end loop over harmonics
end
