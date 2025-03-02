function [T]=gentloc(ti,tf)
%GENTLOC   [T]=gentloc(ti,tf)
%   Generate the transformation matrix from basic to local coodrinates
%   ulocal= T*ubasic
%
%   ti= #nodes x 1 index vector referencing transformation location
%   tf= #coord syst x 9 matrix for storing 3x3 transf matrices
%   T= #nodes*6 x #nodes*6 basic to local transformation matrix
%

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

%  History:
%   APR  1997   akissil   created

[nnodes,mt]=size(ti);
T=sparse(nnodes*6,nnodes*6);

for i=1:nnodes
	if ti(i,1)==0
		t=eye(3);
	else
  		t=reshape(tf(ti(i,1),:),3,3);
	end
	dofs=[((i-1)*6+1):((i-1)*6+3)];
	T(dofs,dofs)=t;
	dofs=[((i-1)*6+4):((i-1)*6+6)];
	T(dofs,dofs)=t;
end




