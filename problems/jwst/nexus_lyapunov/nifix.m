function  ni=nifix(index,oldni,varargin)
%NIFIX	replaces nodal numbers in connectivity matrix
%	with nodal array row indexes.
%
%		ni=nifix(index,oldni,cols)
%
%	ni	fixed connectivities matrix.
%	index	nodal mapping index (see nodesort.m).
%	oldni	old connectivities matrix.
%	cols    specifies the columns that nifix operates on.
%		(optional)
%		
%
%If cols is not specified, nifix substitutes for all but the first
%and last columns.  These are assumed to be the element number and
%property array, respectively.
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
%   ~Jun93 jmelody:   created 
%  17Jun93 jmelody:   removed element types
%  18Mar96 jmelody:   able to deal with empty oldni
%  14Feb97 jmelody:   added optional cols input. also added node 
%		      error checking
%  15Feb97 jmelody:   got rid of error checking, causes problems 
%		      with triangular plate elements
%  12May98  lneedels: changed call to complement to the matlab5 setdiff 
%                     command.  Order of input arguments had to be swapped.
%  14May98 lneedels:  modified for varargin for matlab5
%

[r,c]=size(oldni);

if (nargin == 2)
  cols=[2:c-1]; 	%all but first and last columns
else,
  cols=[varargin{1}];
end

if (r > 0)
  ni=zeros(size(oldni));

  for n=1:r			%loop over rows
    for m=1:length(cols)	%loop over columns
      if (oldni(n,cols(m)) > 0)
        ni(n,cols(m))=index(oldni(n,cols(m)));
      end
    end
  end

%see History note dated 15Fab97
%  [x,y]=find(ni(:,cols) == 0);
%  if (length(x) > 0)
%    error('Node numbers do not exist');
%  end

  kcols=setdiff([1:c],cols);
  ni(:,kcols)=oldni(:,kcols);
else
  ni=[];
end
