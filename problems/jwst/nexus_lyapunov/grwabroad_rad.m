function S=grwabroad_rad(f,Ro,dR,type,model,scale,hrad,crad);
%rwabroad_rad	broadband PSD of HST RWA radial forces
%
%	S=grwabroad_rad(f,Ro,dR,type,model,scale,hrad,crad)
%
%	where,	S	is the Power Spectral Density of the
%		 	reaction wheel radial force disturbance (N^2/Hz)
%		Ro	is the nominal wheel speed in RPM (optional)
%		dR	is the expected variation about the nominal (optional)
%		f	is the vector of frequencies (in Hz) for which 
%		 	S is calculated
%     type    type of distribution assumed for RWA speed
%     model   reaction wheel model to be chosen
%     scale   scale factor for disturbance magnitude
%     hrad    radial harmonics
%     crad    radial amplitude coefficients
%
%       If type == 0 or is not specified, the PSD is calculated on the
%	assumption that the RWA speed is a uniform random variable over the
%	interval [Ro-dR, Ro+dR] rpm.  If dR is not specified, the interval
%	[0, Ro] is assumed.  If neither Ro nor dR are specified, then the
%	interval is assumed to be [0 3000].
%       If type == 1, then the PSD is calculated on the assumption
%	that the RWA speed is a gaussian RV with mean of Ro and standard
%	deviation of dR.  In this case, both Ro and dR must be specified.

%History
%  12Aug94 jmelody:   created
%   5Oct94 jmelody:   fixed -2*pi factor error
%   5Oct94 jmelody:   got rid of fancy stuff in algorithm
%   5Oct94 jmelody:   added max rpm as an input
%  17Nov94 jmelody:   made interval [Ro-dR, Ro+dR]
%  20Nov94 jmelody:   made robust for either row or column frequency vector
%  26Nov94 jmelody:   added choice of RPM a Gaussian RV
%  13Dec94 jmelody:   fixed bug if ml and mu are both [], it used to barf
%  05May98 odeweck:   modified min points at high frequency for plotting
%  08Aug98 odeweck:   added model variable to choose reaction wheel model
%  10Oct99 odeweck:   modified version for use with RWA disturbance model GUI

if (nargin == 1)
  Ru=3000;
  Rl=0;
  type=0;
elseif (nargin == 2)
  Rl=0;
  Ru=Ro;
  type=0;
elseif (nargin == 3)
  Rl=Ro-dR;
  Ru=Ro+dR;
  type=0;
else
  Rl=Ro-dR;
  Ru=Ro+dR;
end

[Ci,hi]=grwadist_rad(1,model,scale,hrad,crad);

%need to convert hi to rad/s
hi=2*pi*hi;

[n,m]=size(f);
if (n > m)
  w=2*pi*f';
else
  w=2*pi*f;
end

S=zeros(size(w));
 
mini = [0] ;
for n=1:max(size(Ci))
  if (type == 0)        %Uniform RV
    ml=find((w > -hi(n)*Ru)&(w < -hi(n)*Rl));
    mu=find((w < hi(n)*Ru)&(w > hi(n)*Rl));
    m=[ml mu];
    S2=pi*Ci(n)^2/(2*(Ru-Rl)*hi(n)^5)*w.^4;
    if (length(m) > 0)
      S(m)=S(m)+S2(m);
      mini=min([S(m) mini]);
    end
  else          %Gaussian
    S=S+pi*Ci(n)^2/(hi(n)^4*sqrt(8*pi*dR^2))...
        *(w.^4.*exp(-(w/hi(n)-Ro).^2/(2*dR^2))...
        +w.^4.*exp(-(-w/hi(n)-Ro).^2/(2*dR^2)));
  end
end

if (type == 0)
%  do this for good plotting
        for n=1:size(S,1),
        if S(n)<mini
           S(n)=mini;
        end
        end
   S(n)=S(1);
   S=S+mini/1e25*ones(size(S));
end


