function [ug]=mce_dis(nset,mset,rg,un)
% MCE_DIS [ug]=mce_dis(nset,mset,rg,un)
%    Solves for the overall deflections, given nset displacements
%    accounting for the effects of MPCs and rigid body elements.
%
%    nset is the vector of degrees of freedom in the n-set (not mpc)
%    mset is the vector of degrees of freedom in the m-set (mpc's)
%    rg is the matrix of constrain coefficients
%    un is the nset displacement vector (or vectors
%    ug is the resultant displacement vector (or vectors)
%
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
%  17Dec93 jmelody:  adapted from Bob Norton's mce
%  12Feb97 lneedels: removed size(un) statement

%Partition the rg matrix into the rm and rn sub-matrices
rm=rg(:,mset);
rn=rg(:,nset);

%Form the gm constraint equation
gm=-inv(rm)*rn;

%Solve for the dependent displacements
um=gm*un;

%Merge the m-set and n-set displacements
ug(nset,:)=un;
ug(mset,:)=um;


