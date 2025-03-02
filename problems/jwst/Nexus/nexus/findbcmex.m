function nodes=findbcmex(bc,dofs,varargin)

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
%  13Sep93 jmelody:  created
%   6Oct93 jmelody:  took out arbitrary node numbering stuff
%  10Oct93 jmelody:  put in arbitrary node numbering stuff
%  13Feb97 lneedels: modified for compilation, rewrote most of function
%  17Dec97 lneedels: changed max(size( to length for matlab5
%  31Dec97 lneedels: added check to make sure dof is there
%  14May98 lneedels: modified for varargin for matlab5
%

nodes=zeros(length(dofs),2);

for n=1:length(dofs)
  [id,jd]=find(bc == dofs(n));
  if length(id)>0 & length(jd)>0,
    nodes(n,:)=[id jd];
  else,
    nodes(n,:)=[0 0];
  end
end

if nargin == 3		%non-sequencial nodal numbering is used
  xyz=[varargin{1}];
  [mxyz,nxyz]=size(xyz);
  if nxyz==4,
    nodes(:,1)=xyz(nodes(:,1),1);
  end
end
