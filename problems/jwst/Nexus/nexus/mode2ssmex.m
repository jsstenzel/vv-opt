function [a,b,c,d]=mode2ssmex(omega,phi,zeta,lm,lb,lc)

% 
% Copyright 1992.  National Aeronautics and Space Administration,
% all rights reserved.
% Copyright 1993-2000.  California Institute of Technology, National
% Aeronautics and Space Administration.  All rights reserved.  US Government
% Sponsorship acknowledged.
% 
% THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTIBILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
% CALTECH BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
% ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
% WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE ACCESS, USE OF
% PERFORMANCE OF THIS SOFTWARE.
% 

%History
%   6Oct93 jmelody:  omega and zeta don't have to both be rows or
%		     columns
%  18Feb86 lneedels: modified for compilation
%  17Dec97 lneedels: changed max(size( to length for matlab5
%

[no,mo]=size(omega);
[nz,mz]=size(zeta);
 
if nz~=no       %if zeta and omega aren't both rows or columns
  zeta=zeta';   %make them both rows or columns
end

nlm=length(lm);
nlb=length(lb);
nlc=length(lc);
a=[zeros(nlm,nlm),eye(nlm,nlm);...
   diag(-omega(lm).*omega(lm)), diag(-2*zeta(lm).*omega(lm))];
b=[zeros(nlm,nlb);phi(lb,lm)'];
c=[zeros(nlc,nlm),phi(lc,lm)];
d=zeros(nlc,nlb);
