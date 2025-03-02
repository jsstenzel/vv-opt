function [a,b,c,d,lb,lc]=mode2ss2(xyz,bc,nm,ig,dg,vg,nrbm,za,phi,om)
%KM2SS [a,b,c,d,lb,lc]=mode2ss2(xyz,bc,nm,ig,dg,vg,nrbm,za,phi,om)
%       Computes the state space model from the modal data
%
%       xyz     = nodal coordinate array (n x 3, or n x 4 with free numbering)
%       bc      = degree of freedom array (n x 6, with 0 for constrained dofs)
%       nm      = number of modes to compute and use in state space model
%       ig   = (:,7) input grid numbers(:,1) and six 0 or 1 dof flags (:,2:7)
%       dg   = (:,7) displacement output grid nunmbers(:,1) and six 0 or 1 dof
%                 flags(:,2:7) ([] is ok)
%       vg   = (:,7) velocity outputs (same format as dg, [] is ok)
%       nrbm   = number of rigid body modes already in phi and om
%       za     = damping vector or value
%                if za is a vector, then it will be used as is
%                if za is a value, then it will be used for elast modes 
%                if za is [], then use 0.1% for all elastic modes
%       gm     = constraint eq coefficients (optional)
%       phi     = g-size eigenvector matrix (nm columns, mass normalized)
%       om    = column vector of modal frequencies (rad/sec, nm x 1)
%       a,b,c,d = state space model
%  Note: c matrix will have displacements first, then velocities

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

%  HISTORY    APR  1997   akissil   created

        [ng,mg]=size(xyz);
        bci=ones(ng,6);
        [bci,ndof]=bcond(bci);

	lb=getdofs(ig,xyz,bci);
        dc=[];
        vc=[];
        if ~isempty(dg)
	 dc=getdofs(dg,xyz,bci);
        end
        if ~isempty(vg)
         vc=getdofs(vg,xyz,bci);
        end
% stack output requests, displacements first
        lc=[dc; vc];
       

% generate modal damping vector
        if length(za)>1
          if length(za)~=nm
            error('za damping vector is wrong size')
          end
          zeta=za;
        elseif isempty(za)
          zeta= ones(nm,1)*0.001;
          if nrbm>0
            zeta(1:nrbm,1)=zeros(nrbm,1);
          end
        else
          zeta= ones(nm,1)*za;
          if nrbm>0
            zeta(1:nrbm,1)=zeros(nrbm,1);
          end
        end

	[a,b,c,d]=mode2ss(om,phi,zeta,[1:nm],lb,lc);

% all velocities are default for c
        if ~isempty(dc)        
          if isempty(vc)
% if only displacements are requested
     	   c(:,[nm+1:2*nm,1:nm])=c;
           return
          end
% modify c to get displacements for specific outputs
          dset=[1:max(size(dc))];
    	  c(dset,[nm+1:2*nm,1:nm])=c(dset,:);
        end

