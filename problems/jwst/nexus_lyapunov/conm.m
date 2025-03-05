function mo=conm(em,nodes,bc,m,varargin)
% CONM    mo = conm(em,nodes,bc,m,xyz)    assembles the concentrated mass, em, 
%    into the system mass matrix, m, for the nodes listed in nodes.
%        
%    em		is the concentrated mass.
%    nodes	is a vector of node numbers where the mass is placed.
%    bc		is the array of dofs.
%    m		is the system mass matrix.
%    xyz        nodal coordinates (needed only if using arbitrary
%		node numbering).
%    mo		is the returned mass matrix.
%
%    The method handles various sizes of concentrated mass matrices.
%
%        If em is (1,1), it is added to the active translational dofs
%        for each node.
%
%        If em is a vector, it is assumed to be the diagonal of a square
%        matrix that is the element mass matrix for each node.  If the em
%        is larger than the active dofs for a node, then the
%        inactive mass is disregarded.
%
%        If em is a matrix (it must be square), then it is added to the
%        active dofs of the nodes as above.
%
%        In order to add an inertia to a single dof (or several for that
%matter), use a vector em with all zero elements, except for those dofs
%which you want to add mass to.
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

%History
%   5Oct93 jmelody:  expanded error checking on size of nodes (make sure
%		     nodes is a matrix)
%   5Oct93 jmelody:  made more robust for vector and matrrix em input
%  10Oct93 jmelody:  arbitrary node numbering capabilities added
%  29Jan97 lneedels: modified for compiler
%  14May98 lneedels: modified for varargin for matlab5
%

% nodes must be a row vector so transpose if necessary.  If nodes is a
% matrix then error will result.

[nnodes,mnodes] = size(nodes);
if min(nnodes,mnodes) > 1
  error('The input variable nodes must be a vector, not a matrix!');
end
if nnodes > mnodes;
  nodes = nodes';
  [nnodes,mnodes] = size(nodes);
end;

[nem,mem]=size(em);
if nem == 1             %if em is a vector, make sure it's a column
  if mem ~= 1
    em=em';
    [nem,mem]=size(em);
  end
end
 
if (nem ~=1 & mem~=1)
  if nem~=mem
    error('The matrix, em, is not square!');
  end
end

if nargin == 5          %xyz is included
  xyz=[varargin{1}];
  [nxyz,mxyz]=size(xyz);
  if mxyz == 4          %non-sequencial node numbering
    index=nodesort(xyz);
    for i=1:mnodes
      nodes(1,i)=index(nodes(i));
    end
  end
end

if issparse(em),  %em was passed in as a sparse matrix
  em=full(em);
end

[ivec,jvec,mvec]=conmmex(em,nodes,bc);

% get size of m
[mdof,ndof]=size(m);

mo=m+sparse(ivec,jvec,mvec,mdof,ndof);
