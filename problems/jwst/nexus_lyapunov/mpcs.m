function [nset,mset,rg]=mpcs(bc,xyz,ndof,nset,mset,rg,mpc)
% MPCS [nset,mset,rg]=mpcs(bc,xyz,ndof,nset,mset,rg,mpc)
%   Creates a single mpc relationship.
%
%   bc is the (:,6) array of dof numbers.
%   xyz is the (:,3) or (:,4) array of nodal coordinates.
%   ndof is the number of dofs in the problem.
%   nset is the vector of dofs in the n-set (not mpc).
%   mset is the vector of dofs in the m-set (mpcs).
%   rg is the matrix of constraint coefficients.
%   mpc is the (:,3) matrix of the mpc definition.
%     mpc(:,1) is the grid point.
%     mpc(:,2) is the component number (1 thru 6).
%     mpc(:,3) is the coefficient.
%     The first row of this vector defines the dependent dof.  
%     The sum over all the dofs is zero:
%             sum(Ai*Ui)=0  where Ai is the coefficient
%                           and Ui is the dof.
%  The dependent dof must be in the n-set before
%  calling this function.  The independent dofs
%  can be in the n-set, the m-set, or null.
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

%History
%  Written by Bob Norton, 8/18/93.
%  17Dec93 jmelody:  got rid of quit statement for error handling
%  21Dec93 jmelody:  worked on bug for independent grid with not all
%                    six dofs active
%   4Feb97 lneedels: modified for compilation
%

[mtmp,ntmp]=size(nset);
if mtmp>1 & ntmp>1,
  error('nset must be a vector!')
end
if ntmp==1,
  nset=nset';
end

[mtmp,ntmp]=size(mset);
if mtmp>1 & ntmp>1,
  error('mset must be a vector!')
end
if ntmp==1,
  mset=mset';
end

% Check the size of [rg].  If there is only 1 row, then check
% to see if it is null.  Increment the number of constraint equations.
% If there is one row, check to see if it is just a zero.

[nrow,dum]=size(rg);
if nrow==1
  if rg*rg'==0
    ncon=0;
  else
    ncon=1;
  end
  if dum==1,
    if rg==0,
      rg=sparse(1,ndof);
    end
  end
elseif nrow==0,
  rg=sparse(1,ndof);
  ncon=nrow;
else
  ncon=nrow;
end
ncon=ncon+1;

[nset,mset,rgmat]=mpcsmex(bc,xyz,ndof,nset,mset,mpc);

rg(ncon,:)=sparse(1,rgmat(:,1),rgmat(:,2),1,ndof);

