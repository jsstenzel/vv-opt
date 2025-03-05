function mo=conm2(ni,bc,m,varargin)
% CONM2  mo=conm2(ni,bc,m,xyz,cid,cf) Assembles concentrated masses and moments/products of inertia
%    into the system mass matrix, m.  Offsets of the cg from the grid point can be specified.
%        
%    ni     is a (:,3),(:,6),(:,9), or (:,12), each row describing a mass element
%           ni(:,1)  is the element id number
%           ni(:,2)  is the grid point id number 
%        *** when ni is a (:,3),(:,6),(:,9),(:,12), the following definitions are used:
%           ni(:,3)  is the mass value
%           -------------------------------------
%           ni(:,4)  is the 1st offset component from the grid point to the cg (optional)
%           ni(:,5)  is the 2nd offset component from the grid point to the cg (optional)
%           ni(:,6)  is the 3rd offset component from the grid point to the cg (optional)
%           -------------------------------------
%           ni(:,7)  is the I11 mass moment of inertia about the cg (optional, see note)
%           ni(:,8)  is the I22 mass moment of inertia about the cg (optional, see note)
%           ni(:,9)  is the I33 mass moment of inertia about the cg (optional, see note)
%           -------------------------------------
%           ni(:,10) is the I21 mass product of inertia about the cg (optional, see note)
%           ni(:,11) is the I31 mass product of inertia about the cg (optional, see note)
%           ni(:,12) is the I32 mass product of inertia about the cg (optional, see note)
%    bc		is the array of degree of freedom numbers.
%    m		is the system mass matrix & is returned modified in mo.
%    xyz        nodal coordinates (needed only if using arbitrary
%		node numbering)
%    cid        is a vector of coordinate system id numbers for inertia matrix definition
%               must have same number of rows as ni if present (optional, 0=basic)
%    cf         is a (:,12) array of coordinate sys definition data (see coord_in,optional)
%    mo		is the returned mass matrix.
%
%    Note:  signs of product of inertia input terms are consistent with the following
%           form of Inertia Matrix (neg signs added by func, i.e. same convention as NASTRAN):
%
%                        |  I11  -I21  -I31  |
%                        |        I22  -I32  |
%                        |  sym         I33  |

% Copyright 1993-1998.  California Institute of Technology, National
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
%  05Jan00   akissil:  corrected bug for arbitrary node numbering case
%  19May00   akissil:  took out explicit matrix input (see conm1)

[nn,mn]= size(ni);

if nargin > 3          %xyz is included
  xyz=[varargin{1}];
  [nxyz,mxyz]=size(xyz);
  if mxyz == 4          %non-sequencial node numbering
    index=nodesort(xyz);
    xyz=xyz(:,2:4);     %this line was added to correct bug, 05Jan00
    for i=1:nn
      ni(i,2)=index(ni(i,2));
    end
  end
end

% Note: we will use non-arbitrary node numbering format for conm calls (no xyz)
%       so ni(:,2) is renumbered to proper sequence and xyz will be a (:,3)


% compute transformation matrices based on coord system and grid location
cid= zeros(nn,1);
if nargin > 4
  cid= [varargin{2}];
  [nrw,nco]=size(cid);
  if nco>nrw
    cid=cid';
  end
  cf= [varargin{3}];
  if length(cid(:,1))~=nn
    error('The cid input vector for conm2 should have the same No. of rows as ni')
  end
  [ti,tf]= coord_in([ni(:,2),zeros(nn,1),cid],cf,xyz);
end
cidmx= max(cid);

if mn==3  
% no mass offsets or moments of inertia
  for j=1:nn
    m = conm(ni(j,3),ni(j,2),bc,m);
  end

elseif mn==6
% offsets present, but no moments of inertia input
  for j=1:nn
    em= zeros(3,3);	
	em(1,2)= ni(j,6)*ni(j,3);
	em(1,3)=-ni(j,5)*ni(j,3);
	em(2,3)= ni(j,4)*ni(j,3);
	em(2,1)=-em(1,2);
	em(3,1)=-em(1,3);
	em(3,2)=-em(2,3);
	ei(1,1)=  ni(j,3)*(ni(j,5)^2+ni(j,6)^2);
	ei(2,2)=  ni(j,3)*(ni(j,4)^2+ni(j,6)^2);
	ei(3,3)=  ni(j,3)*(ni(j,4)^2+ni(j,5)^2);
	ei(1,2)=-(ni(j,3)*(ni(j,4)*ni(j,5)));
	ei(1,3)=-(ni(j,3)*(ni(j,4)*ni(j,6)));
	ei(2,3)=-(ni(j,3)*(ni(j,5)*ni(j,6)));
	ei(2,1)= ei(1,2);
	ei(3,1)= ei(1,3);
	ei(3,2)= ei(2,3);
	em= [diag(ones(3,1)*ni(j,3)) em; em' ei];

% transform mass matrix to basic system if cid is present	
        if cidmx>0
          if ti(ni(j,2),1)~=0
            T= zeros(6,6);
  	    t=reshape(tf(ti(ni(j,2),1),:),3,3);
	    dofs=[1:3];
	    T(dofs,dofs)=t;
	    dofs=[4:6];
	    T(dofs,dofs)=t;
	    em= T'*em*T;
	  end
	end
	  
    m = conm(em,ni(j,2),bc,m);
  end

elseif mn==9
% offsets and moments of inertias present
  for j=1:nn
    em= zeros(3,3);	
	em(1,2)= ni(j,6)*ni(j,3);
	em(1,3)=-ni(j,5)*ni(j,3);
	em(2,3)= ni(j,4)*ni(j,3);
	em(2,1)=-em(1,2);
	em(3,1)=-em(1,3);
	em(3,2)=-em(2,3);
%    ei= zeros(3,3);	
	ei(1,1)=  ni(j,7) + ni(j,3)*(ni(j,5)^2+ni(j,6)^2);
	ei(2,2)=  ni(j,8) + ni(j,3)*(ni(j,4)^2+ni(j,6)^2);
	ei(3,3)=  ni(j,9) + ni(j,3)*(ni(j,4)^2+ni(j,5)^2);
	ei(1,2)=-(ni(j,3)*(ni(j,4)*ni(j,5)));
	ei(1,3)=-(ni(j,3)*(ni(j,4)*ni(j,6)));
	ei(2,3)=-(ni(j,3)*(ni(j,5)*ni(j,6)));
	ei(2,1)= ei(1,2);
	ei(3,1)= ei(1,3);
	ei(3,2)= ei(2,3);
	em= [diag(ones(3,1)*ni(j,3)) em; em' ei];

% transform mass matrix to basic system if cid is present	
        if cidmx>0
          if ti(ni(j,2),1)~=0
            T= zeros(6,6);
  	    t=reshape(tf(ti(ni(j,2),1),:),3,3);
	    dofs=[1:3];
	    T(dofs,dofs)=t;
	    dofs=[4:6];
	    T(dofs,dofs)=t;
	    em= T'*em*T;
	  end
	end

    m = conm(em,ni(j,2),bc,m);
  end

elseif mn==12
% offsets and moments and products of inertias present
  for j=1:nn
    em= zeros(3,3);	
	em(1,2)= ni(j,6)*ni(j,3);
	em(1,3)=-ni(j,5)*ni(j,3);
	em(2,3)= ni(j,4)*ni(j,3);
	em(2,1)=-em(1,2);
	em(3,1)=-em(1,3);
	em(3,2)=-em(2,3);
	ei(1,1)=  ni(j,7) + ni(j,3)*(ni(j,5)^2+ni(j,6)^2);
	ei(2,2)=  ni(j,8) + ni(j,3)*(ni(j,4)^2+ni(j,6)^2);
	ei(3,3)=  ni(j,9) + ni(j,3)*(ni(j,4)^2+ni(j,5)^2);
	ei(1,2)=-(ni(j,10) + ni(j,3)*(ni(j,4)*ni(j,5)));
	ei(1,3)=-(ni(j,11) + ni(j,3)*(ni(j,4)*ni(j,6)));
	ei(2,3)=-(ni(j,12) + ni(j,3)*(ni(j,5)*ni(j,6)));
	ei(2,1)= ei(1,2);
	ei(3,1)= ei(1,3);
	ei(3,2)= ei(2,3);
	em= [diag(ones(3,1)*ni(j,3)) em; em' ei];

% transform mass matrix to basic system if cid is present	
        if cidmx>0
          if ti(ni(j,2),1)~=0
            T= zeros(6,6);
  	    t=reshape(tf(ti(ni(j,2),1),:),3,3);
	    dofs=[1:3];
	    T(dofs,dofs)=t;
	    dofs=[4:6];
	    T(dofs,dofs)=t;
	    em= T'*em*T;
	  end
	end
    m = conm(em,ni(j,2),bc,m);
  end

else

  error('ni array for conm2 must be a (:,3), (:,6), (:,9), or (:,12)')
end

mo= m;
