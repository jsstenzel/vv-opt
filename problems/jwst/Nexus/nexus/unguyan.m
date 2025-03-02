function dout=unguyan(toa,din,aset,oset)
% UNGUYAN    dout=unguyan(toa,din,aset,oset) computes the full displacement
%	vector from the guyan reduced displacement vector.
%
%   toa  is the transformation matrix from the guyan reduction
%        where  d(oset) = toa * d(aset).
%
%   din  is the d(aset), aset displacements.
%
%   aset is the vector of retained dofs.
%
%   oset is the vector of omitted dofs.
%
%   dout is the full displacement set, d(oset+aset).
%       See guyan for more information.
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

[na,temp]=size(aset);
[no,temp]=size(oset);
[temp,nlc]=size(din);
dout=zeros(no+na,nlc);
dout(aset,:)=din;
dout(oset,:)=toa*din;
