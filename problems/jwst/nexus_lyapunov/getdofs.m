function dofs=getdofs(gid,xyz,bc)
%GETDOFS dofs=getdofs(gid,xyz,bc)
%        Retrieves vector of dof numbers (dofs), given a matrix (gid) of
%        grid point numbers with six 0 or 1 flag indicators.
%        gid= (:,7) matrix, (:,1) is grid id, (:,2:7) are 0 or 1 flags
%        xyz= (:,3 or 4) grid location array
%        bc = (:,6) array of degree of freedom numbers
%        dofs= column vector of degree of freedom numbers

% Version 2.0 1997, Copyright California Institute of Technology, National
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

%  HISTORY  APR  1997   akissil   created
%
%	11Nov97  gmosier  changed 'if X == []' to 'if isempty(X)' for v5
%
	ndof=0;
	[ng,mg]=size(gid);
        [nx,mx]=size(xyz);
        if mx==4
         gv=xyz(:,1);
        else
         gv=[1:nx]';
        end

	for i=1:ng
		j=find(gv==gid(i,1));
%		if j==[]
		if isempty(j)
			error('reference to nonexistent grid point')
		end
		for k=1:6
			if gid(i,k+1)==1
if bc(j,k)==0
 error('Selected dof is not active in bc')
end
			  ndof=ndof+1;
			  dofs(ndof,1)= bc(j,k);
			end
		end
	end

