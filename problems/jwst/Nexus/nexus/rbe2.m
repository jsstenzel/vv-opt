function [nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,varargin)
% RBE2 [nset,mset,rg]=rbe2(bc,xyz,nset,mset,rg,gn,cm,gm,ti,tf)
%   Creates a single RBE2 rigid element.
% 
%   bc is the (:,6) array of dof numbers.
%   xyz is the (:,3) or (:,4) array of nodal coordinates.
%   nset is the vector of dofs in the n-set (not mpc).
%   mset is the vector of dofs in the m-set (mpc's).
%   rg is the matrix of constraint coefficients.
%   gn is the independent grid (all 6 d.o.f.).
%   cm is a vector of the dependent component numbers. Use integers 
%     from 1 through 6.
%   gm is the list of dependent grids.
%   ti is the indexing vector for storage of cs transf data in tf (optional). 
%   tf is the matrix storing local cs transf submatrices (optional) .
%
%   The dependent dofs must be in the n-set before
%   calling this function.  The independent dofs
%   can be in the n-set, the m-set, or null.
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
%   Written by Bob Norton, 12/4/93.
%  21Dec93 jmelody:  changed so that cm and gm can be 
%  		     row or column vectors
%  21Dec93 jmelody:  worked on bug for n-set with not all 
%		     six dofs active
%  31Jan94 jmelody:  changed so that ndof is not an input variable
%  31Jan94 jmelody:  fixed the reinitialization bug with the array mpc
%   5Feb97 lneedels: modified for compilation
%  16Mar98 lneedels: major overhaul with local coord ala akissil
%  14May98 lneedels: modified for varargin for matlab5

% Check the number of dependent component numbers

ndof=max(max(bc));

if nargin<10
 ti=[];
 tf=[];
else
  ti=[varargin{1}];
  tf=[varargin{2}];
  if ~isempty(ti)&length(xyz(:,1))~=length(ti(:,1))
    error('Local cs vector ti should have same no of rows as xyz')
  end
end

[nc,dum]=size(cm);
if nc == 1	%row vector
  cm=cm';	%switch to column vector
  [nc,dum]=size(cm);
end
% Check the number of dependent grid points
[ndep,dum]=size(gm);
if ndep == 1	%row vector
  gm=gm';	%switch to column vector
  [ndep,dum]=size(gm);
end

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
else
  ncon=nrow;
end
ncon=ncon+nc*ndep;

% Check the size of [xyz].  If there are 4 columns, then the first
% column is the node label and free nodal numbering is being used.
% If only 3 columns are present,then free nodal numbering isn't being used.
[numnod,nxyz]=size(xyz);

if nxyz==3
  indgr=gn;
else
  indgr=find(xyz(:,1)==gn);
end

xyz_ind(1)=xyz(indgr,nxyz-2);	
xyz_ind(2)=xyz(indgr,nxyz-1);
xyz_ind(3)=xyz(indgr,nxyz);

if ~isempty(ti)
%  retrieve the cs transformation matrix for the independent grid
 if ti(indgr,1)==0
	tind=eye(3);
 else
  	tind=reshape(tf(ti(indgr,1),:),3,3);
 end
 tind=[tind zeros(3,3); zeros(3,3) tind];
end

% Loop over all the dependent grids
for i=1:ndep
  if nxyz==3
    depgr=gm(i);
  else
    depgr=find(xyz(:,1)==gm(i));
  end
  xyz_dep(1)=xyz(depgr,nxyz-2);
  xyz_dep(2)=xyz(depgr,nxyz-1);
  xyz_dep(3)=xyz(depgr,nxyz);

%  compute the transformation matrix (T) from independent to dependent dofs
%  where   udep(bas)= T*uind(bas)
 D=[0  xyz_dep(3)-xyz_ind(3) xyz_ind(2)-xyz_dep(2)
    xyz_ind(3)-xyz_dep(3) 0  xyz_dep(1)-xyz_ind(1)
    xyz_dep(2)-xyz_ind(2) xyz_ind(1)-xyz_dep(1)  0];
 T=[eye(3) D; zeros(3,3) eye(3)];

 if ~isempty(ti)
%  retrieve the cs transformation matrix for the dependent grid
  if ti(depgr,1)==0
	tdep=eye(3);
  else
  	tdep=reshape(tf(ti(depgr,1),:),3,3);
  end
  tdep=[tdep zeros(3,3); zeros(3,3) tdep];
%  transform to local coordinates
%  udep(loc)= tdep*T*tind'*uind(loc)
  T= tdep*T*tind';
 end

% Build the appropriate mpc equations
% loop through dependent components: one mpc equation
% for each component.
 for mc=1:nc
    cnt=1;
    mpc=zeros(2,3);
    mpc(1,1)=gm(i);
    mpc(1,2)=cm(mc); 
    mpc(1,3)=-1;

% find all nonzero indep dof coefficients for dependent component
    for j=1:6
     if T(cm(mc),j)~=0
      cnt=cnt+1;
      mpc(cnt,1)=gn;
      mpc(cnt,2)=j;
      mpc(cnt,3)=T(cm(mc),j);
     end
    end
    [nset,mset,rg]=mpcs(bc,xyz,ndof,nset,mset,rg,mpc);
 end
% end of dependent component loop

end
% end of dependent grid loop

