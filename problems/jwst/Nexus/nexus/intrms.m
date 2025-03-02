function [rmsout]=intrms(psd,f);
% intrms   calculates rms from a given psd
%          takes psd in volts^2/hz and frequency vector and
%          integrates under curve, returns rms in volts

%  Tupper Hyde 22Jun94

nf=length(f);  df=f(2:nf)-f(1:nf-1); df=df(:);
avgs=(psd(2:nf)+psd(1:nf-1))/2;
rmsout=sqrt(sum(avgs.*df));
