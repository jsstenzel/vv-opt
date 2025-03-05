function [K]=celas(kelem,g1,c1,g2,c2,bc,K,varargin);
%CELAS	Scalar elastic element.
%
%	[K]=celas(kelem,g1,c1,g2,c2,bc,K,xyz);
%
%	where	kelem	is the spring rate of the scalar element.
%		g1	is the first grid point.
%		c1	is a vector of the affected dofs for the first 
%                       node (contains the integers 1 thru 6).
%		g2	is the second grid point.
%		c2	is a vector of the affected dofs for the second
%                       node (contains the integers 1 thru 6).
%		bc	is the array of dofs.
%		K	is the system stiffness matrix.
%		xyz	nodal coordinates (needed only if using arbitrary
%			node numbering).
%
%	In order to connect a node to ground with a spring, use zero 
%for either g1 or g2.
%
%	If c1 and c2 are vectors, then the corresponding dofs are 
%connected pairwise by a scalar springs specified by kelem.  If kelem
%is a vector, the corresponding element of kelem is used to connect
%the dofs in c1 and c2.  If kelem is a scalar, this value is used to
%connect all of the dofs specified by c1 and c2.  If kelem is a matrix,
%(size(c1)-by-size(c2)), then this matrix is added whole into system
%stiffness matrix.

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
%  17Jul94 jmelody:   created
%  28Jan97 lneedels:  modified for compilation
%   4Apr97 lneedels:  added check to make sure c1 & c2 have values 
%		      between 1-6 and modified "help" summary.
%  17Dec97 lneedels:  changed max(size( to length for matlab5
%  14May98 lneedels:  modified for varargin for matlab5


%check sizes of c1, c2 and kelem
if (min(size(c1)) ~= 1)
  error('c1 must be either a scalar or a vector');
end
if (min(size(c2)) ~= 1)
  error('c2 must be either a scalar or a vector');
end
if (length(c1) ~= length(c2))
  error('c1 and c2 must have the same dimension');
end
if (length(kelem) ~= 1)	%kelem is not a scalar
  if (min(size(kelem)) ~= 1)		%kelem is a matrix
    [nk,mk]=size(kelem);
    nc1=length(c1);
    nc2=length(c2);
    if ((nk ~= nc1)|(mk ~=nc2))
      error('The matrix kelem is the wrong size');
    end
  elseif (length(kelem) ~= length(c2))
    error('If kelem is a vector, it must have the same dimension as c1 and c2');
  end
else		%kelem is a scalar
  kelem=kelem*ones(size(c1));		%change into a vector
end

%make sure that c1 and c2 are between 1&6
if ((max(c1) > 6) | (min(c1) < 1))
  error('c1 must contain only integers with values between 1 and 6');
end
if ((max(c2) > 6) | (min(c2) < 1))
  error('c2 must contain only integers with values between 1 and 6');
end

if (g1==0 & g2==0),
    error('CELAS: warning, you are trying to connect ground to ground (g1=g2=0)');
end

%check for nodal numbering
if nargin == 8		%xyz is input
  xyz=[varargin{1}];
  [ivec,jvec,kvec]=celasmex(kelem,g1,c1,g2,c2,bc,xyz);
else
  [ivec,jvec,kvec]=celasmex(kelem,g1,c1,g2,c2,bc);
end

%  get size of K
[mdof,ndof]=size(K);

K=K+sparse(ivec,jvec,kvec,mdof,ndof);
