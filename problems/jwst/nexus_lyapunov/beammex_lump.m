function [ivec,jvec,kvec,mvec]=beammex_lump(ni,xyz,prop,mat,bc,ti,tf)

% Copyright 1992.  National Aeronautics and Space Administration,
% all rights reserved.
% Copyright 1993-1998.  California Institute of Technology, National
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
 
%HISTORY
%  14Apr92 hcb:       for sign error in local-to-global coord xform
%  17Jul92 BNorton:   use NASTRAN style consistent mass matrix
%                     and correct an error in the sign of the 2,6
%                     and 3,5 elements.
%  13Aug92 BNorton:   add shear flexibility
%  26Aug92 jmelody:   add nodal numbering
%  17Jun94 jmelody:   removed element type from calls to nifix.m
%  24Jun96 BNorton:   added pin flag capability (allows release of selected
%                     degrees of freedom at either end of the beam element)
%  16Jan96 lneedels:  major remoding, to make it suitable for matlab
%                     compiler
%  14Feb97 jmelody:   added cols input to nifix.m in order to fix bug 
%                     with simultaneous use of pin flags and arbitrary node 
%                     numbers.  Bug found by prapacz.
%  21Feb97 lneedels:  fixed bug with multiple loopings over i
%  26Mar97 lneedels:  removed local transformations
%  03Mar98 lneedels:  modified to use NASTRAN-like properties
%  14Mar98 lneedels:  modified to include local cs & orient vector ala akissil
%  12May98 lneedels:  modified find's with row vector input to be Matlab5.1+
%                     compliant
%  16Jun98 rnorton:   added nsm/length capability
%  24Sep99 akissil:   changed to lumped mass formulation
%  12Apr99 akissil:   cleaned up for formal inclusion
%  27Jul99 akissil:   corrected overzealous commenting-out of pinflag stuff
%  27Jul99 akissil:   modify = to >= (kb) per Laura Needels e-mail (April 23, 1999)

% size of ni
[msz,nni]=size(ni);

% preallocate maximum possible length of vectors
ivec=zeros(msz*12*12,1);
jvec=ivec;
kvec=ivec;
mvec=ivec;
cntr=0;

for i=1:msz;
	
  ek=zeros(12,12);
  np=ni(i,2);
  nm=ni(i,3);

  if np==0 | nm==0
    disp('Problem with element number ')
    disp(ni(i,1))
    error('Beam element references non-existant grid number')
  end
  n3=ni(i,4);
  if nni >=6
    ka=ni(i,6);
  else
    ka=0;
  end
  if nni>=7
    kb=ni(i,7);
  else
    kb=0;
  end
  v1=[xyz(nm,:)-xyz(np,:)];
  len=sqrt(v1*v1');
  if len==0
    fprintf('The problem is with element Number %7.0f\n',ni(i,1))
    error('Zero length beam element encountered')
  end
  v1=v1'/len;
  if n3==0
    tni=0;
    if ~isempty(ti)
      tni = ti(np);
    end
    temp=ni(i,8:10);
    if tni > 0
      t1=reshape(tf(tni,:),3,3);
      temp=[t1'*temp']';
    end
  else
    temp=[xyz(n3,:)-xyz(np,:)];
  end
  if norm(temp)==0
    fprintf('The problem is with element Number %7.0f\n',ni(i,1))
    error('Zero length orientation vector for beam element encountered')
  end
  temp=temp'/sqrt(temp*temp');
  v3=[v1(2)*temp(3)-v1(3)*temp(2);...
  v1(3)*temp(1)-v1(1)*temp(3);...
  v1(1)*temp(2)-v1(2)*temp(1)];
  if norm(v3)==0
    fprintf('The problem is with element Number %7.0f\n',ni(i,1))
    error('Orientation vector is in same direction as beam element')
  end
  v3=v3/sqrt(v3'*v3);
  v2=[v3(2)*v1(3)-v3(3)*v1(2);...
  v3(3)*v1(1)-v3(1)*v1(3);...
  v3(1)*v1(2)-v3(2)*v1(1)];
  v2=v2/sqrt(v2'*v2);
  t=[v1';v2';v3'];

  c=zeros(12,12);
  c([1:3],[1:3])=t;
  c([4:6],[4:6])=t;
  c([7:9],[7:9])=t;
  c([10:12],[10:12])=t;
  
  pid=ni(i,5);
%  find which row of prop has pid
  proprow=find(prop(:,1)==pid);
%  find which row of mat has mid
  matrow=find(mat(:,1)==prop(proprow,2));
%  check to make sure this is a mat1 card
  if mat(matrow,2)~=1,
    error('beam.m - beam elements must use mat1 cards')
  end

  nsm=prop(proprow,7);  %  nsm/element
  nsm2=prop(proprow,8);  %  nsm/length
  E=mat(matrow,3);
  G=mat(matrow,4);
  A=prop(proprow,3);
  rhoA=mat(matrow,6)*A;
  k1=prop(proprow,9);
  k2=prop(proprow,10);
  I22=prop(proprow,4);
  I33=prop(proprow,5);
  J=prop(proprow,6);

% define some common terms:
  lsq=len*len;
  lcube=len*len*len;
  ei1=E*I22;
  ei2=E*I33;
  gak1=G*A*k1;
  gak2=G*A*k2;
  if k1==0
    r1=12*ei1/lcube;
  else
    r1=(12*ei1*gak1)/(gak1*lcube+12*len*ei1);
  end
  if k2==0 
    r2=12*ei2/lcube;
  else
    r2=(12*ei2*gak2)/(gak2*lcube+12*len*ei2);
  end;
  sk1=r1*lsq/4 + ei1/len;
  sk2=r2*lsq/4 + ei2/len;
  sk3=r1*lsq/4 - ei1/len;
  sk4=r2*lsq/4 - ei2/len;
  ael=A*E/len;
  lr1=len*r1/2;
  lr2=len*r2/2;
  gjl=G*J/len;

  ek(1,1)=ael;
  ek(1,7)=-ael;
  ek(2,2)=r1;
  ek(2,6)=lr1;
  ek(2,8)=-r1;
  ek(2,12)=lr1;
  ek(3,3)=r2;
  ek(3,5)=-lr2;
  ek(3,9)=-r2;
  ek(3,11)=-lr2;
  ek(4,4)=gjl;
  ek(4,10)=-gjl;
  ek(5,3)=-lr2;
  ek(5,5)=sk2;
  ek(5,9)=lr2;
  ek(5,11)=sk4;
  ek(6,2)=lr1;
  ek(6,6)=sk1;
  ek(6,8)=-lr1;
  ek(6,12)=sk3;
  ek(7,1)=-ael;
  ek(7,7)=ael;
  ek(8,2)=-r1;
  ek(8,6)=-lr1;
  ek(8,8)=r1;
  ek(8,12)=-lr1;
  ek(9,3)=-r2;
  ek(9,5)=lr2;
  ek(9,9)=r2;
  ek(9,11)=lr2;
  ek(10,4)=-gjl;
  ek(10,10)=gjl;
  ek(11,3)=-lr2;
  ek(11,5)=sk4;
  ek(11,9)=lr2;
  ek(11,11)=sk2;
  ek(12,2)=lr1;
  ek(12,6)=sk3;
  ek(12,8)=-lr1;
  ek(12,12)=sk1;
  
%  const=((rhoA+nsm2)*len+nsm)/420;
% Use the consistent mass matrix formulation from NASTRAN:
%  bl22=22*len;
%  bl13=13*len;
%  blsq4=4*len*len;
%  blsq3=3*len*len;
%  em=zeros(12,12);
%  em(1,1)=175;
%  em(1,7)=35;
%  em(7,1)=35;
%  em(2,2)=156;
%  em(2,6)=bl22;
%  em(6,2)=bl22;
%  em(2,8)=54;
%  em(8,2)=54;
%  em(2,12)=-bl13;
%  em(12,2)=-bl13;
%  em(3,3)=156;
%  em(3,5)=-bl22;
%  em(5,3)=-bl22;
%  em(3,9)=54;
%  em(9,3)=54;
%  em(3,11)=bl13;
%  em(11,3)=bl13;
%  em(5,5)=blsq4;
%  em(5,9)=-bl13;
%  em(9,5)=-bl13;
%  em(5,11)=-blsq3;
%  em(11,5)=-blsq3;
%  em(6,6)=blsq4;
%  em(6,8)=bl13;
%  em(8,6)=bl13;
%  em(6,12)=-blsq3;
%  em(12,6)=-blsq3;
%  em(7,7)=175;
%  em(8,8)=156;
%  em(8,12)=-bl22;
%  em(12,8)=-bl22;
%  em(9,9)=156;
%  em(9,11)=bl22;
%  em(11,9)=bl22;
%  em(11,11)=blsq4;
%  em(12,12)=blsq4;
% Multiply the matrix so far by the constant.
%  em=const*em;
% To get the torsional inertia about the longitudinal axis, 
% assume that all the mass is located at the radius of gyration
% of the cross-section.  Since the radius of gyration is the 
% square root of the moment of inertia/area and there are two
% moments of inertia, use the average of the two.
%  if A==0
%   rad2=0;
%  else
%   rad2=(I22+I33)/(2*A);
%  end
%  em(4,4)=210*const*rad2;
%  em(10,10)=210*const*rad2;

% Use the pin flags, if present, to modify the stiffness matrix
  if (ka+kb) > 0
    ipin=zeros(1,10);
    ekp=zeros(12,12);
    emp=zeros(12,12);
    for idx1=1:5
      ipin(idx1)=rem(ka,10);
      ipin(idx1+5)=rem(kb,10)+6;
      if ipin(idx1+5)==6
        ipin(idx1+5)=0;
      end
      ka=fix(ka/10);
      kb=fix(kb/10);
    end
    for idx1=1:10
      if ipin(idx1)>0
        ii=ipin(idx1);
        if ek(ii,ii)==0
          for j=1:12
            ek(j,ii)=0;
            ek(ii,j)=0;
          end
        else
          for j=1:12
            for ll=1:12
              ekp(ll,j)=ek(ll,j)-(ek(ll,ii)/ek(ii,ii))*ek(ii,j);
            end
            ekp(j,ii)=0;
            ekp(ii,j)=0;
          end
          ek=ekp;
        end
      end
    end
% Now modify the mass matrix
%    for j=1:10
%      if ipin(j)>0
%        jj=ipin(j);
%        if ek(jj,jj)~=0
%          for idx1=1:12
%            for l1=1:12
%              emp(l1,idx1)=em(l1,idx1)-ek(jj,l1)*em(idx1,jj)/ek(jj,jj)...
%                       -ek(idx1,jj)*em(jj,l1)/ek(jj,jj)...
%                       +ek(jj,l1)*em(idx1,jj)*em(jj,jj)/ek(jj,jj)^2;
%            end
%          end
%          em=emp;
%        end
%        for kk=1:12
%          em(jj,kk)=0;
%          em(kk,jj)=0;
%        end
%      end
%    end
  end
% akissil mod- added next 8 lines to give lumped mass form
  const=((rhoA+nsm2)*len+nsm)*0.5;
  em=zeros(12,12);
  em(1,1)= const;
  em(2,2)= const;
  em(3,3)= const;
  em(7,7)= const;
  em(8,8)= const;
  em(9,9)= const;
% Transfrom the mass and stiffness matrices from local
% to global

  ek=c'*ek*c;
  em=c'*em*c; 
    
% assembly of vectors which store m&k infomation
  id=[bc(np,1:6) bc(nm,1:6)];
  [idx,jdx,vdx]=find(id);   % active degrees of freedom
  sz=length(jdx);

  for idx=1:sz,
    ivec(cntr+(idx-1)*sz+1:cntr+(idx-1)*sz+sz,1)=[vdx]';
    jvec(cntr+(idx-1)*sz+1:cntr+(idx-1)*sz+sz,1)=vdx(idx)*ones(sz,1);
  end
  kvec(cntr+1:cntr+sz*sz)=reshape(ek(jdx,jdx),sz*sz,1);
  mvec(cntr+1:cntr+sz*sz)=reshape(em(jdx,jdx),sz*sz,1);
  cntr=cntr+sz*sz;

end;

%  trim output vectors
ivec=ivec(1:cntr);
jvec=jvec(1:cntr);
kvec=kvec(1:cntr);
mvec=mvec(1:cntr);
