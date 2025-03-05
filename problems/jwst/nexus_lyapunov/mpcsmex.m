function [nseto,mseto,rgmat]=mpcsmex(bc,xyz,ndof,nset,mset,mpc)

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
%  17Dec93 jmelody:  got rid of quit statement for error handling
%  21Dec93 jmelody:  worked on bug for independent grid with not all
%                    six dofs active
%   4Feb97 lneedels: modified for compilation, including moving the 
%                    node from nset to mset
%

% Check the size of [xyz].  If there are 4 columns, then the first
% column is the node label and free nodal numbering is being used.
% If only 3 columns are present, then free nodal numbering isn't being used.
[numnod,nxyz]=size(xyz);
% The first degree of freedom will be dependent upon the rest,
% so check to make sure that it is in the current n-set,
% i.e. it cannot already be dependent or be constrained.
% First check that the specified degree of freedom is
% defined in bc:
if nxyz==3
  node=mpc(1,1);
else
  node=find(xyz(:,1)==mpc(1,1));
end

ncomp=mpc(1,2);
if bc(node,ncomp)==0
  disp('The following grid and component is specified')
  disp('as the dependent degree of freedom on an mpc:')
  disp('Grid:')
  disp(mpc(1,1))
  disp('Component:')
  disp(mpc(1,2))
  disp('This degree of freedom is fixed in the bc matrix.')
  error('***Fatal Error***')
end
% Second, check to see if this degree of freedom is present
% in the n-set.
inset=find(nset(:)==bc(node,ncomp));
if size(inset)==0
  disp('The following grid and component is specified')
  disp('as the dependent degree of freedom on an mpc:')
  disp('Grid:')
  disp(mpc(1,1))
  disp('Component:')
  disp(mpc(1,2))
  disp('This degree of freedom is not present in the')
  disp('independent set of degrees-of-freedom.')
  error('***Fatal Error***')
end
% Now transfer this degree of freedom from the n-set
% to the m-set

[dum,ns]=size(nset);
nseto=[nset(1:inset-1) nset(inset+1:ns)];

[dum,ns]=size(mset);
if ns==1
  if mset(1)==0
    mset=[];
  end
end
mseto=sort([mset bc(node,ncomp)]);

% Construct the constraint equation
% First, null out this row of the constraint equation matrix
[nterms,dum]=size(mpc);
rgmat=zeros(nterms,2);
cntr=0;
for i=1:nterms
  if nxyz==3
    node=mpc(i,1);
  else
    node=find(xyz(:,1)==mpc(i,1));
  end
  dof=bc(node,mpc(i,2));
  if dof ~= 0		%This fixes bug for when independent node has
			%less than six dofs active: jmelody
    rgmat(cntr+1,:)=[dof mpc(i,3)];
    cntr=cntr+1;
  end
end
rgmat=rgmat(1:cntr,:);
