function [bcn]=bcnset(bc,nset,mset)
%BCNSET  [bcn]=bcnset(bc,nset,mset)
%
%   Assembles a new dof matrix which has only nset active dofs
%
%   bc is the input array of degrees of freedom (dofs).
%   nset is the vector of independent dofs.
%   mset is the vector of dependent dofs.
%   bcn is the new array of degrees of freedom output.
%

% Copyright 1998-2000.  California Institute of Technology, National
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
%  00Apr97 akissil:  created
%  25Feb98 lneedels: changed max(size()) to length() for matlab 5
%  25Feb98 lneedels: did some cleanup work

if length(mset)<1,
  error('Don''t use bcnset if mset is empty!')
end

[nbc,mbc]=size(bc);
ndof=length(nset) + length(mset);
if ndof~=max(max(bc))
  error('The size of nset+mset does not match ndof in bc')
end

%find location of dofs that we want to zero out
nodes=findbcmex(bc,mset);
idx=find(nodes(:,1)>0);

% calculate values of location in matrix 
jdx=(nodes(idx,2)-1)*nbc+nodes(idx,1);

% reset values
bc(jdx)=zeros(length(idx),1);

%  bcond only checks to see if value is >0 to be active
[bcn,ndof]=bcondmex(bc);

