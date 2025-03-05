function [gm,knn,mnn,varargout]=mce1(nset,mset,rg,k,m,varargin)
% MCE1 [gm,knn,mnn,cnn]=mce1(nset,mset,rg,k,m,c)
%    Creates the overall constraint equation after
%    all the MPC relations and rigid elements have been entered.
%
%    nset is the vector of dofs in the n-set (not mpc).
%    mset is the vector of dofs in the m-set (mpc's).
%    rg is the matrix of constraint coefficients.
%    k is the stiffness matrix.
%    m is the mass matrix.
%    c is the damping matrix (optional).
%    gm is the m-set sized constraint equation.
%    knn is the n-set sized stiffness matrix.
%    mnn is the n-set sized mass matrix.
%    cnn is the n-set sized damping matrix (optional).
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

%History
%    Written by Bob Norton, 9/11/93.
%   1Feb94 jmelody:  changed so that if mset is empty or zero, return
%		     without changing m and k
%   2Dec94 jmelody:  added damping matrix as input and output
%  06Jan98 lneedels: fixed bug in check to see if mset=0
%  15May98 lneedels: modified for varargin for matlab5
%

%if mset is zero or empty set, return
if (length(mset) == 0)
  mset=[0];
end
if (length(mset) == 1)&(mset(1) == 0)
  knn=k;	%set values and return
  mnn=m;
  if (nargin == 6)
    cnn=c;
  end
  gm=[];
  return
end

% Partition the rg matrix into the rm and rn sub-matrices
rm=rg(:,mset);
rn=rg(:,nset);

% Form the gm constraint equation
gm=-inv(rm)*rn;

% Partition the stiffness matrix
knb=k(nset,nset);
knm=k(nset,mset);
kmm=k(mset,mset);
knn=knb+knm*gm+gm'*knm'+gm'*kmm*gm;

% Partition the mass matrix
mnb=m(nset,nset);
mnm=m(nset,mset);
mmm=m(mset,mset);
mnn=mnb+mnm*gm+gm'*mnm'+gm'*mmm*gm;

if (nargin == 6)	%damping matrix included
  c=[varargin{1}];
  cnb=c(nset,nset);
  cnm=c(nset,mset);
  cmm=c(mset,mset);
  cnn=cnb+cnm*gm+gm'*cnm'+gm'*cmm*gm;
  varargout{1}=cnn;
end

gm=sparse(gm);

