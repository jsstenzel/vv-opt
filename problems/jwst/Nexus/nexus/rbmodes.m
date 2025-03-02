function  phi=rbmodes(bc,xyz,m,varargin)
%RBMODES	calculates rigid body modeshapes geometrically.
%
%  phi=rbmodes(bc,xyz,m,Rpt,ti,tf);
%
%  bc is the array of dofs.
%  xyz is the nodal coordinate matrix.
%  m is the system mass matrix.
%  Rpt is the rotation point (optional).
%  ti is the #nodes x 1 index vector referencing transform location (optional).
%  tf is the #coord syst x 9 matrix for storing 3x3 transf matrices (optional).
%  phi is the rigid body mode eigenvectors.
%
%	If no rotation point is specified, the cm is used as the center
%of rotation.
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
%  21Aug92 jmelody:  created
%   3Feb94 jmelody:  removed M as input and added ndof as input.
%	             This is useful for cases when xyz and bc
%		     are subsetted for component "rigid body" modes
%  07Mar97 lneedels: modified for compilation
%  14May98 lneedels: modified for varargin for matlab5
%  09Jun98 akissil : modifications for local coordinate systems-NOT BACKWARD 
%                    COMPATIBLE
%  25Jun98 lneedels: modified for varargin for matlab5
%

[mxyz,nxyz]=size(xyz);

if nxyz == 4		%non-sequencial nodal numbering
  xyz=xyz(:,2:4);	%throw out nodal number column
end

if nargin == 3,
  Rpt=[];
  ti=[];
  tf=[];
elseif nargin==4,
  Rpt=[varargin{1}];
  ti=[];
  tf=[];
elseif nargin==5,
  error('Must include both ti & tf if you want to use local coordinates')
elseif nargin==6,
  Rpt=[varargin{1}];
  ti=[varargin{2}];
  tf=[varargin{3}];
end
  
if isempty(Rpt)
  Rpt=cg_calc(m,xyz,bc,ti,tf);
end

phi=wtcgmex(bc,xyz,Rpt,ti,tf);

