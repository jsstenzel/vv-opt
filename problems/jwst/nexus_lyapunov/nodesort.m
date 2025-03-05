function index=nodesort(inp)
%NODESORT	index=nodesort(inp)
%               provides an index matrix for a matrix
%		with id numbers.
%
%		index = nodesort(inp)
%
%		inp	input matrix where inp(:,1) are
%			positive (this means non-zero also)
%			unique id numbers.
%		index	nodal index vector.
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
%   ~Jun92 jmelody:   created
%   5Jul94 jmelody:   chnged error printing from disp() to error()
%  27Dec95 lneedels:  modified so it runs much faster and added 
%                     sparse matrix benefits
%  12Jan96 lneedels:  removed size checking of incoming matrix
%  12Mar97 lneedels:  modified for compilation
%  23Dec98 lneedels:  modified so only first column of inp is passed in
%

[z,ix]=nodesortmex(inp(:,1));

index=sparse(z,1,ix,max(z),1);
