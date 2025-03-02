function phi=normphi(phi,M);
%NORMPHI	mass normalizes FEM modeshapes.
%
%		newphi=normphi(oldphi,m);
%
%	where,	newphi	normalized eigenvectors.
%		oldphi	non-normalized eigenvectors.
%		m	system mass matrix.
%
%Note:	oldphi can be a subset of the eigenvectors.
%

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

[m,n]=size(M);
if (m~=n)
  error('Mass Matrix must be square');
end

[p,q]=size(phi);
if (p~=m)
  error('Eigenvectors are not the proper length');
end

for n=1:q
  phi(:,n)=phi(:,n)/sqrt(phi(:,n)'*M*phi(:,n));
end

