function [ivec,jvec,kvec]=celasmex(kelem,g1,c1,g2,c2,bc,varargin)

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
%  17Dec97 lneedels:  changed max(size( to length for matlab5
%  12May98 lneedels:  modified find's with row vector input to be Matlab5.1+
%                     compliant
%  14May98 lneedels:  modified for varargin for matlab5
%  15May98 lneedels:  changed "disp" to "warning" where appropriate
%  30Dec98 lneedels:  moved warning to end, so that you only get it no
%                     more than once per call to celas
%  30Dec98 lneedels:  added ability to use kelem as a matrix and 
%                     connect it to ground (g1=0 or g2=0)
%

% preallocate maximum possible length of vectors
[mtmp,ntmp]=size(kelem);

if min([mtmp ntmp])==1, % kelem is vector/scalar
  ivec=zeros(length(kelem)*4*length(c1),1);
else, % kelem is a matrix
  ivec=zeros(mtmp*ntmp*4*length(c1),1);
end
jvec=ivec;
kvec=ivec;
mvec=ivec;
cntr=0;
iwarn=0;

%check for nodal numbering
if nargin == 7		%xyz is input
  xyz=[varargin{1}];
  [temp,n]=size(xyz);
  if n == 4	%nodal numbering is used
    if (g1 ~= 0)
      g1=find(xyz(:,1) == g1);
    end
    if (g2 ~= 0)
      g2=find(xyz(:,1) == g2);
    end
  end
end

if (min(size(kelem)) == 1)	%kelem is not a matrix
  nc1=length(c1);
  for n=1:nc1
    dofs=[0 0];
    if g1>0,
      dofs(1)=bc(g1,c1(n));
    end
    if g2>0,
      dofs(2)=bc(g2,c2(n));
    end

    if ((dofs(1) == 0)&(dofs(2) == 0)),
      iwarn=1;
      sz=0;
    elseif (dofs(1)==0),
      ivec(cntr+1,1)=dofs(2);
      jvec(cntr+1,1)=dofs(2);
      kvec(cntr+1,1)=kelem(n);
      cntr=cntr+1;
    elseif (dofs(2)==0),
      ivec(cntr+1,1)=dofs(1);
      jvec(cntr+1,1)=dofs(1);
      kvec(cntr+1,1)=kelem(n);
      cntr=cntr+1;
    else,
      ivec(cntr+1:cntr+4,1)=[dofs dofs]';
      jvec(cntr+1:cntr+4,1)=[dofs(1) dofs(1) dofs(2) dofs(2)]';
      kvec(cntr+1:cntr+4,1)=kelem(n)*[1 -1 -1 1]';
      cntr=cntr+4;
    end

  end
else			%kelem is a matrix

  ek=[kelem -kelem;-kelem kelem];

  dofs=[];
  if g1>0,
    dofs=[bc(g1,c1)];
  end
  if g2>0,
    dofs=[dofs bc(g2,c2)];
  end

  [idx,jdx,vdx]=find(dofs);   % active degrees of freedom
  sz=length(jdx);

  for idx=1:sz,
    ivec(cntr+(idx-1)*sz+1:cntr+(idx-1)*sz+sz,1)=[vdx]';
    jvec(cntr+(idx-1)*sz+1:cntr+(idx-1)*sz+sz,1)=vdx(idx)*ones(sz,1);
  end
  kvec(cntr+1:cntr+sz*sz)=reshape(ek(jdx,jdx),sz*sz,1);
  cntr=cntr+sz*sz;

end

if iwarn==1,
  warning('You are trying to connect to inactive dofs');
end

%  trim output vectors
ivec=ivec(1:cntr);
jvec=jvec(1:cntr);
kvec=kvec(1:cntr);

