function [K,G]=celas2(ni,xyz,bc,K,G)
%CELASN	Scalar elastic spring element using connectivity array
%
%        [K,G]=celas2(ni,xyz,bc,K,G)
%
%    where  ni is the element connectivity array
%             ni(:,1) is the element number
%             ni(:,2) is the node number of the first end of the element
%             ni(:,3) is the component number for the first node (1-6)
%             ni(:,4) is the node number of the second end of the element
%             ni(:,5) is the component number for the second node (1-6)
%             ni(:,6) is the value for the spring stiffness 
%             ni(:,7) is the value for the structural damping factor (gfac, 0 is ok)
%                     (element structural damping is Ge= gfac*Ke)
%           xyz is the nodal coordinate array (:,3 or 4)
%		    bc	is the array of degree of freedom numbers.
%		    K	is the system stiffness matrix.
%		    G	is the system structural damping matrix (optional).
%
%	Note: In order to connect a node to ground with a spring, use zero 
%         for either ni(i,2) or ni(i,4).

% Copyright 1992.  National Aeronautics and Space Administration,
% all rights reserved.
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
%  11Jun99   akissil:  created

if length(ni(1,:))~=7
  error('ni array in celas2 elements must have 7 columns')
end

if length(xyz(1,:))==3
  xyz= [[1:length(xyz(:,1))]' xyz];
end

for j=1:length(ni(:,1))
  kelem= ni(j,6);
  gfac=  ni(j,7);
%  add spring stiffness
  Kout= celas(kelem,ni(j,2),ni(j,3),ni(j,4),ni(j,5),bc,K,xyz);

%  add structural damping if present  
  if gfac>0
    if ~exist('G')
	  disp(['celas1 element ',num2str(ni(j,1))]);
	  error('Nonzero Structural Damping Factor Defined but No G Matrix Input')
	end
    G= G + gfac*[Kout-K];
  end
  K= Kout;
end
  