function [ti,tfo,xyz]=coord_in(ci,cf,xyzi)
%COORD_IN  [ti,tfo,xyz]=coord_in(ci,cf,xyzi)
%  Computes the local coordinate system transformations used for 
%  input of grid point locations as well as for grid displacements.
%
%  ci is a (:,3) array of grid point and coordinate system ids.
%   ci(:,1) is the node number.
%   ci(:,2) is the input coordinate sys id (for grid locations). (0=basic).
%   ci(:,3) is the output coordinate sys id (for grid displacements). (0=basic).
%
%  cf is a (:,12) array of coordinate sys definition data.
%   cf(:,1) is the coordinate system id.
%   cf(:,2) is the coordinate system type.
%          1= rectangular (x,y,z)
%          2= cylindrical (r,theta,z)
%          3= spherical (r,theta,phi) (theta from z-axis and
%                                      phi from x-axis in xy plane)
%   cf(:,3) is the reference coordinate system id used for a,b,c
%          point locations (0=basic).
%          a= the origin of the new coordinate system.
%          b= a point on the z-axis of the new coordinate system.
%          c= a point in the +zx quadrant.
%   cf(:,4) is a1, 1st coordinate of point a in reference coord system.
%   cf(:,5) is a2, 2nd coordinate of point a in reference coord system.
%   cf(:,6) is a3, 3rd coordinate of point a in reference coord system.
%   cf(:,7) is b1, 1st coordinate of point b in reference coord system.
%   cf(:,8) is b2, 2nd coordinate of point b in reference coord system.
%   cf(:,9) is b3, 3rd coordinate of point b in reference coord system.
%   cf(:,10) is c1, 1st coordinate of point c in reference coord system.
%   cf(:,11) is c2, 2nd coordinate of point c in reference coord system.
%   cf(:,12) is c3, 3rd coordinate of point c in reference coord system.
%
%  xyzi is the (:,3) or (:,4) array of nodal coordinates given 
%         in the input coordinate systems.
%
%  ti is a (#nodes,1) array of indices to the rows of tf for output.
%          ti(:,1) is the row number of tf containing the transform data.
%          0= basic coordinate system.
%
%  tfo is a (:,9) array of coordinate system transformation data.
%   tfo(:,1:9) is the 3x3 transformation matrix stored row-wise such
%            that when t=reshape() is used, ulocal=t*ubasic.
%  xyz is the (:,3) or (:,4) array of nodal coordinates transformed 
%         to the basic coordinate system.

% 
% Copyright 1998-2000.  California Institute of Technology, National
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
%  04Oct96  akissil:  created
%  13May98 lneedels:  changed "if tu==[]" to be matlab5 compliant
%  30Aug99  akissil:  added mods to handle undefined of tfo
% 

dr=pi/180;
[cin,cim]=size(ci);
[cfn,cfm]=size(cf);
[xyzn,xyzm]=size(xyzi);
%   error checking
if cim~=3
  error('coord_in ****  ci should have 3 columns  ****')
end
if cfm~=12
  error('****  cf should have 12 columns  ****')
end
if max(cf(:,2))>3 | min(cf(:,2))<1
  error('coord_in ****  coordinate system types can only be 1, 2 or 3')
end

xyzin=xyzi(:,1:3);

%if arbitrary nodal numbering is used then coordinates are shifted
if xyzm==4
  xyzin=xyzi(:,2:4);
end

%count the coordinate systems as they are stored
cnt=0;

%rows of cfxtf correspond to rows of cf and contain tf row numbers
cfxtf=zeros(cfn,1);

%loop through coordinate systems iteratively to
%handle daisey chains

iter=0;
while cnt<cfn

  iter= iter+1;
  if iter>cfn
    error('check cf reference coordinate system definitions')
  end
  for i=1:cfn

%check if cs is not stored yet
    if cfxtf(i,1)==0

%check if reference coordinate system is defined yet
      cref=find(cf(:,1)==cf(i,3));

      if (cf(i,3)==0) | (cfxtf(cref,1)~=0)

%process new coordinate system and store the system type
        cnt= cnt + 1;
        cfxtf(i,1)=cnt;
        tf(cnt,10)=cf(i,2);

%if the reference coordinate system is not the basic, then
%transform a, b & c to basic
        if cf(i,3)~=0

%retrieve reference coordinate system transformation matrix
%transpose on reshape gives t such that ubasic=t*ulocal

          t=reshape(tf(cfxtf(cref,1),1:9),3,3)';

%if the reference coordinate system is rectangular, then 
%proceed with the transformation

          if cf(cref,2)==1
            d=cf(i,4:6)';
            a=t*d;
            d=cf(i,7:9)';
            b=t*d;
            d=cf(i,10:12)';
            c=t*d;

%if the ref coord system is cylindrical then we must convert
%a, b & c from (r, theta, z) to (x, y, z)

          elseif cf(cref,2) ==2
            d=[cf(i,4)*cos(cf(i,5)*dr)
            cf(i,4)*sin(cf(i,5)*dr)
            cf(i,6)];
            a=t*d;
            d=[cf(i,7)*cos(cf(i,8)*dr)
            cf(i,7)*sin(cf(i,8)*dr)
            cf(i,9)];
            b=t*d;
            d=[cf(i,10)*cos(cf(i,11)*dr)
            cf(i,10)*sin(cf(i,11)*dr)
            cf(i,12)];
            c=t*d;

%if the new coord system is spherical then we must convert
%a, b & c from (r, theta, phi) to (x, y, z) using sph function

          elseif cf(cref,2)==3
            d=[cf(i,4)*sin(cf(i,5)*dr)*cos(cf(i,6)*dr)
            cf(i,4)*sin(cf(i,5)*dr)*sin(cf(i,6)*dr)
            cf(i,4)*cos(cf(i,5)*dr)];
            a=t*d;
            d=[cf(i,7)*sin(cf(i,8)*dr)*cos(cf(i,9)*dr)
            cf(i,7)*sin(cf(i,8)*dr)*sin(cf(i,9)*dr)
            cf(i,7)*cos(cf(i,8)*dr)];
            b=t*d;
            d=[cf(i,10)*sin(cf(i,11)*dr)*cos(cf(i,12)*dr)
            cf(i,10)*sin(cf(i,11)*dr)*sin(cf(i,12)*dr)
            cf(i,10)*cos(cf(i,11)*dr)];
            c=t*d;
            
          end

%add the basic coordinates of the reference origin and
%save the converted a, b & c back in the original cf array

          d=cf(cref,4:6);
          cf(i,4:12)=[a'+d b'+d c'+d];
        end

%compute the basic to local transformation matrix
%ahat is along local x, bhat along local y, chat along local z
%dhat from the local origin to c
        a=cf(i,4:6);
        b=cf(i,7:9);
        c=cf(i,10:12);
        if norm(b-a)==0
          error('points a and b of cs definition are coincident')
        end
        chat=(b-a)/norm(b-a);
        dhat=c-a;
        if dhat==0
          error('points a and c of cs definition are coincident')
        end
%bhat is chat cross dhat

        bhat=[chat(1,2)*dhat(1,3)-chat(1,3)*dhat(1,2);
        chat(1,3)*dhat(1,1)-chat(1,1)*dhat(1,3);
        chat(1,1)*dhat(1,2)-chat(1,2)*dhat(1,1)]';
        if norm(bhat)==0
          error('zero length for local y-axis')
        end
        bhat=bhat/norm(bhat);
%ahat is bhat cross chat
        ahat=[bhat(1,2)*chat(1,3)-bhat(1,3)*chat(1,2);
        bhat(1,3)*chat(1,1)-bhat(1,1)*chat(1,3);
        bhat(1,1)*chat(1,2)-bhat(1,2)*chat(1,1)]';
        if norm(ahat)==0
          error('zero length for local x-axis')
        end
        ahat=ahat/norm(ahat);

%store transformation matrix
        t=[ahat
        bhat 
        chat]';
        tf(cnt,1:9)=[t(1,:) t(2,:) t(3,:)];
      end
    end
  end
end

%determine nodal output coordinate system indices (0= basic)

xyz=xyzin;
ti=zeros(xyzn,1);
tused=zeros(xyzn,2);
cnto=0;
for i=1:cin
  if xyzm==4
    nd=find(xyzi(:,1)==ci(i,1));
  else
    nd=ci(i,1);
  end
%		if ci(i,3)~=0
%			ti(nd,1)=cfxtf(find(cf(:,1)==ci(i,3)));
%		end

%retrieve the input coordinate system transformation matrix
%transpose on reshape gives t such that ubasic=t*ulocal

  if ci(i,2)~=0
    cfin=find(cf(:,1)==ci(i,2));
    if isempty(cfin)
      error('Reference to Undefined Input Coord System')
    end
    tfin=cfxtf(cfin,1);
    t=reshape(tf(tfin,1:9),3,3)';
    cstyp=tf(tfin,10);

%perform nodal transformations

    d=xyzin(nd,:)';
    if cstyp==2
      d=[xyzin(nd,1)*cos(xyzin(nd,2)*dr)
      xyzin(nd,1)*sin(xyzin(nd,2)*dr)
      xyzin(nd,3)];
    elseif cstyp==3
      d=[xyzin(nd,1)*sin(xyzin(nd,2)*dr)*cos(xyzin(nd,3)*dr)
      xyzin(nd,1)*sin(xyzin(nd,2)*dr)*sin(xyzin(nd,3)*dr)
      xyzin(nd,1)*cos(xyzin(nd,2)*dr)];
    end
    xyz(nd,:)=[t*d]' + cf(cfin,4:6);
  else
    xyz(nd,:)=xyzin(nd,:);
%	    cstyp=0;
  end

%compute the output transformations

%		if (cstyp==2) | (cstyp==3)
%			ti(nd,2)=xyzin(nd,1);
%		end
  if ci(i,3)~=0
    cfout=find(cf(:,1)==ci(i,3));
    if isempty(cfout)
      error('Reference to Undefined Output Coord System')
    end
    tfout=cfxtf(cfout,1);
    t=reshape(tf(tfout,1:9),3,3);
    cstyp=tf(tfout,10);

%perform nodal transformations

    if cstyp==1
      tu=find(tused(:,1)==tfout);
      if isempty(tu),
        cnto= cnto + 1;
        tfo(cnto,1:9)=tf(tfout,1:9);
        tused(nd,:)=[tfout cnto];
        ti(nd,1)= cnto;
      else
        ti(nd,1)= tused(tu,2);
      end 
    elseif cstyp==2
      cnto= cnto + 1;
      ti(nd,1)= cnto;
      d=t*[xyz(nd,:)-cf(cfout,4:6)]';
      ahat=[t'*[d(1:2,1); 0]]';
      if norm(ahat)==0
        error('grid point located on axis of cyl output coord sys (ambiguity)')
      end
      ahat=ahat/norm(ahat);
      bhat=[t'*[-d(2,1); d(1,1); 0]]';
      bhat=bhat/norm(bhat);
      chat=[ahat(1,2)*bhat(1,3)-ahat(1,3)*bhat(1,2);
      ahat(1,3)*bhat(1,1)-ahat(1,1)*bhat(1,3);
      ahat(1,1)*bhat(1,2)-ahat(1,2)*bhat(1,1)]';
      chat=chat/norm(chat);
%store transformation matrix
      t=[ahat; bhat; chat]';
      tfo(cnto,1:9)=[t(1,:) t(2,:) t(3,:)];
    elseif cstyp==3
      cnto= cnto + 1;
      ti(nd,1)= cnto;
      d=t*[xyz(nd,:)-cf(cfout,4:6)]';
      ahat=[t'*[xyz(nd,:)-cf(cfout,4:6)]']';
      if norm(chat)==0
        error('grid point located on origin of sph output coord sys (ambiguity)')
      end
      ahat=ahat/norm(ahat);
      chat=[t'*[-d(2,1); d(1,1); 0]]';
      if norm(chat)==0
        error('grid point located on axis of sph output coord sys (ambiguity)')
      end
      chat=chat/norm(chat);
      bhat=[chat(1,2)*ahat(1,3)-chat(1,3)*ahat(1,2);
      chat(1,3)*ahat(1,1)-chat(1,1)*ahat(1,3);
      chat(1,1)*ahat(1,2)-chat(1,2)*ahat(1,1)]';
      bhat=bhat/norm(bhat);
%store transformation matrix
      t=[ahat; bhat; chat]';
      tfo(cnto,1:9)=[t(1,:) t(2,:) t(3,:)];
    end

  end
end

%convert back to arbitrary node numbering
if xyzm==4
  xyz=[xyzi(:,1) xyz];
end

% check if tfo has been defined
if exist('tfo')==0,
  tfo=[];
end
