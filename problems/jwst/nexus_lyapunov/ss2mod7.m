function [sysm,T] = ss2mod7(sys,scf)
% function [sysm,T] = ss2mod7(sys,scf)
%
%  file: ss2mod7.
%  written by:  mathieu mercadal 8/13/91
%  coded by:    ken lazarus 8/13/91
%  modified by: doug macmartin 10/1/91, 3/18/92
%               Simon Grocott 20 Jul 1992 to correct transformation for balancing
%               S.G.  5 Dec 1992, put in sort
%               S.G.  fixed sort for  
%  convert state space system to modal tri-diagonal system
%               S.G. 8 Mar 1993 option for second order canonical form (scf=1)
%  NOTE: only works for diagonalizable a, not checked.
%
% ss2mod7: assumes am has an even number of states, puts all eigenvalues 
%          (including real ones) in 2 X 2 blocks. Greg Mallory, 1999/04/20
% scf=0 -> real modal form w real ev's on diagonal
% scf=1 -> 2 X 2 blocks with [0 1; a b]
%
% The returned T is such that sysm=ss2ss( sys,inv(T) )
%
epsilon = 1e-14;
%
% get evals and evecs
%
[a,b,c,d]=ssdata(sys);
[n,n] = size(a);

if (2*floor(n/2) ~= n)
   return
end

[ev,eval] = eig(a);
[lam,k]=sort(diag(eval));
ev=ev(:,k);
T=zeros(n,n);
T2=zeros(n,n);
%
% set pointer and check evals
%
ii = 1;
%
% find real e-values and stack vectors
%
%
% find complex e-values and stack repeated vectors
%
am=zeros(size(a));
indreal=[];
indsep=[];

for i = 1:n,
   if imag(lam(i))>epsilon, % only gets +i part of evalue
      if real(lam(i))>-epsilon*100
         indsep=[indsep;ii;ii+1];
      end
      T(:,ii)    = real(ev(:,i));
      T(:,ii+1)  = imag(ev(:,i));
      am(ii:ii+1,ii:ii+1)=[real(lam(i)) imag(lam(i));-imag(lam(i)) real(lam(i))];
      T2(ii:ii+1,ii:ii+1)=[1 0;real(lam(i)) imag(lam(i))];
      ii = ii +  2;
   elseif imag(lam(i))>-epsilon
      indreal=[indreal;i];
   end
end

%
nr=length(indreal);
ki=1;
while (ki<nr)
   %for i=indreal
   i=indreal(ki);
   ip=indreal(ki+1); % ip=i+1
   if abs(imag(lam(i)))>0,
      if imag(lam(i))>0,
        if real(lam(i))>-epsilon*100
           indsep=[indsep;ii;ii+1];
        end
        T(:,ii)    = real(ev(:,i));
        T(:,ii+1)  = imag(ev(:,i));
        am(ii:ii+1,ii:ii+1)=diag([real(lam(i)) real(lam(i))]);
        T2(ii:ii+1,ii:ii+1)=eye(2);
        ii = ii +  2;
     end
     ki=ki+1;
  else
     %if real(lam(i))>-epsilon*100
      %    indsep=[indsep;ii];
      %end
      %T(:,ii)=real(ev(:,i));
      %am(ii,ii)=real(lam(i));
      %T2(ii,ii)=1;
      %ii = ii + 1;
      if real(lam(i))>-epsilon*100
          indsep=[indsep;ii;ii+1];
      end
      T(:,ii)=real(ev(:,i));
      T(:,ii+1)  = real(ev(:,ip));
      am(ii:ii+1,ii:ii+1)=diag([real(lam(i)) real(lam(ip))]);
      T2(ii:ii+1,ii:ii+1)=[1 1; real(lam(i)) real(lam(ip))];
      ii = ii + 2;
      ki=ki+2;
   end
end

% transform system
%
rT=rank(T);
if rT<n
   %fix for ill-conditioned transformation T
   cm = c*T;
   [U,S,V]=svd(T);
   S=diag(S); iS=[1./S(1:rT)' zeros(1,n-rT)]; iS=diag(iS);
   Ti=V*iS*U'; %pseudoinverse transformation obtained via svd
   bm = Ti*b;
else
  bm = T\b;
  cm = c*T;
end
if scf
  rT=rank(T2);
  if rT<n
  %temporary fix for ill-conditioned transformation T2
   [U,S,V]=svd(T2);
   S=diag(S); iS=[1./S(1:rT)' zeros(1,n-rT)]; iS=diag(iS);
   T2i=V*iS*U';
   am = T2*am*T2i;
   bm = T2*bm;
   cm = cm*T2i;
  else
    am = T2*am/T2;
    bm = T2*bm;
    cm = cm/T2;
  end   
end

normalize=1;
if (normalize)
	%normalize b and c
	indred=[];
	for ii = 1:2:n,
   	 if norm(cm(:,ii:ii+1))>epsilon&norm(bm(ii:ii+1,:))>epsilon,
     	 	fac=sqrt(norm(bm(ii:ii+1,:))/norm(cm(:,ii:ii+1)));
     		bm(ii:ii+1,:)=bm(ii:ii+1,:)/fac;
     		cm(:,ii:ii+1)=cm(:,ii:ii+1)*fac;
      	T(:,ii:ii+1)=T(:,ii:ii+1)*fac;
    	else
     		indred=[indred;ii;ii+1];
    	end  
	end
end

if scf &(rT==n)
   T=T/T2;
elseif scf&(rT<n)
   T=T*T2i;
end
sysm=ss(am,bm,cm,d);
%am=real(am);bm=real(bm);cm=real(cm);







