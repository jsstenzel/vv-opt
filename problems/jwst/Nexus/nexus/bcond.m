function [bco,ndof] = bcond(bc)
% BCOND    [bco,ndof] = bcond(bc)    determines equation or degree
%       of freedom numbers given fixity inputs for nodes
%
%   bc is a (:,6) array of fixity indicaters for each node
%       bc(i,j)>0 means that degree of freedom will be active
%       in the model.  Otherwise the degree of freedom is assigned
%       the value 0
%   ndof is the total number of degrees of freedom
%   bco is the (:,6) array of degree of freedom numbers

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
%   1Oct93 - jmelody   fix zeros(bc) problem
%  28Jan97 - lneedels  modified for compiled version
%
 
[bco,ndof] = bcondmex(bc);
