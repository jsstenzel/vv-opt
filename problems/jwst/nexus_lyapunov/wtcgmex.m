function dgf=wtcgmex(bc,xyz,Rpt,ti,tf)

% This function is used by the wtcg, cg_calc and the rbmodes functions.
% dgf is a f-size set of 6 rigid body mode shapes in local output cs.

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
%  02Sep92 Bob Norton - created
%  02Dec92 jmelody - comments added
%  02Dec92 jmelody - add non-sequencial nodal numbering
%  02Dec92 jmelody - input Rpt as a point in 3-space, not a node
%  02Dec92 jmelody - changed to no input Rpt --> Rpt = cg
%  09Dec92 jmelody - sign error corrected
%  29Oct96 jmelody - incorporated S.Sirlin's bug fix of ndof=max(bc)
%  13Feb97 lneedels - modified for compilation
%  09Jun98 akissil - changed +xyz to -(xyz...) [dg(bc1(i,j),4) for j=2], made 
%                    dg g-size, added accommodation of local coord systems
%  25Jun98 lneedels- made modifications since ti & tf are always defined 
%                    by calling functions
%

[n,temp]=size(bc);

bc1=reshape(1:n*6,6,n)';
ndof=n*6;
dg=zeros(ndof,6);

% define dg which is g-size and basic
for i=1:n
  for j=1:6
    dg(bc1(i,j),j)=1;
    if j==1
      dg(bc1(i,j),5)=dg(bc1(i,j),5)+xyz(i,3)-Rpt(3);
      dg(bc1(i,j),6)=dg(bc1(i,j),6)-(xyz(i,2)-Rpt(2));
    end
    if j==2
      dg(bc1(i,j),4)=dg(bc1(i,j),4)-(xyz(i,3)-Rpt(3));%!!!!!
      dg(bc1(i,j),6)=dg(bc1(i,j),6)+xyz(i,1)-Rpt(1);
    end
    if j==3
      dg(bc1(i,j),4)=dg(bc1(i,j),4)+xyz(i,2)-Rpt(2);
      dg(bc1(i,j),5)=dg(bc1(i,j),5)-(xyz(i,1)-Rpt(1));
    end
  end
end

% transform dg (g-size,basic) to local
if ~isempty(ti)
  t=gentloc(ti,tf);
  dg=t*dg;
end

% now extract the f-set from dg (dgf is f-size and local)
dgf=[];
for i=1:n
  for j=1:6
    if bc(i,j)>0
      dgf=[dgf; dg(bc1(i,j),:)];
    end
  end
end

