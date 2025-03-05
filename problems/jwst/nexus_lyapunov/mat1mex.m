function mat=mat1mex(mid,E,nu,rho,alpha,Tref,mat)

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

%  27Dec95  lneedels:  created
%   3Jan96  lneedels:  removed G as an input, it is now calculated
%   3Jan96  lneedels:  modified mid so that it is guaranteed to be
%			  unique (midp=10*mid+1) for mat1
%  24Jan96  lneedels:  removed the "midp=10*mid+i
%  24Jan96  lneedels:  added the "mat type" column to the mat matrix.
%  24Jan96  lneedels:  added check to make sure mid's are unique.
%  28Jan97  lneedels:  modified for compiled version
%  17Dec97  lneedels:  changed max(size( to length for matlab5
%  03Mar98  lneedels:  added description of location of variables in mat array
%  13Mar98  lneedels:  added/modified the way checks are being done
%  24Jul98  lneedels:  added ability to handle mat=[] initialization
%  21Oct98  lneedels:  added Tref
%

%   1   2   3  4   5  6    7    8
%  mid 1.0  E  G  nu rho alpha Tref

%  calculate G
G=0.5*E./(1.0+nu);

%  size of incoming materials
mrow=length(mid);
if isempty(mat),
  mat=0;
  [mmat,nmat]=size(mat);
else
  [mmat,nmat]=size(mat);
end

% check to see if mat has empty rows
idx=find(mat(:,1)==0);
nzero=length(idx);

%  mat has zeros in "mid" column (empty rows)
if nzero > 0,
  minidx=min(idx);
%  check to make sure zeros are contiguous
  if max(abs(idx-[minidx:minidx+nzero-1]'))>0,
    error('mat1 - noncontiguous zeros in mat matrix')
  end
% check for unique mids
  [z,ix]=unique([mat(1:mmat-nzero,1);mid]);
  if mrow+mmat-nzero ~= length(z)
    error('mid numbers must be unique!')
  end
  mat(minidx:minidx+mrow-1,1:8)=[mid 1.0*ones(mrow,1) E G nu rho alpha Tref];
else
% check for unique mids
  [z,ix]=unique([mat(:,1);mid]);
  if mrow+mmat ~= length(z)
    error('mid numbers must be unique!')
  end
  mat(mmat+1:mmat+mrow,1:8)=[mid 1.0*ones(mrow,1) E G nu rho alpha Tref];
end
 
