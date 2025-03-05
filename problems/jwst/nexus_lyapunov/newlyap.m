function [S] = newlyap(A,B,varargin)
%  S = newlyap(A,B,varargin)   
%
%	Solves for the state covariance matrix S of the matrix state equation
%  q_dot = A*q + B*u using sparseness, by assuming that the A matrix 
%  is 2x2 block diagonal.  Use a third argument to specify the block
%  size of the solutions, otherwise selected for you.
%  Choose to solve the exact solution on 2x2 blocks if the third argument = 1.

% Written by:   Dave Miller   April 11, 2000
% Modified by:  Dave Miller April 26, 2000  -> exact solution added
% 					 Scott Uebelhart June 1, 2000 -> factoring added
%					 Scott Uebelhart July 6, 2000 -> keep only even factors
%               Olivier de Weck July 30, 2000 -> modified from newgram

dbstop if error

[n,n]=size(A);
BB=B*B';

exact = 0 ;

if length(varargin)~=0 & varargin{1}==1,   % if using exact solution on 2x2 blocks
   exact=varargin{1} ;
elseif length(varargin)~=0 & varargin{1}>1,
   m=varargin{1} ;
else,
   % Find possible block sizes and order them for best time to solve each
   % The ideal block sizes (starting between m=20 to m=100) were decided upon after
   % examining empirical tests of different sizes.
   fact = factors(n) ;
   fact = fact(find(mod(fact,2)==0)) ;		% Keep only even values
   fact_ordered = fact([ find(fact>=20 & fact<=100) ...
         find(fact>4 & fact<20) ...
         find(fact>100) ...
         find(fact<4) ]) ;
   m=fact_ordered(1) ;
   %disp(['Chosen block size is m=',num2str(m)])
end



% -------------------  LYAP  ----------------- %
if 0,
   tic
   S=lyap(A,BB);
   time_lyap=toc ;
end


% ----------------- NEW_LYAP ----------------- %
if ~exact,
   S=[];
   tic
   for j=1:n/m
      for jj=j:n/m
         At=(A(m*(j-1)+1:m*j,m*(j-1)+1:m*j));
         Bt=(A(m*(jj-1)+1:m*jj,m*(jj-1)+1:m*jj)');
         Ct=(BB(m*(j-1)+1:m*j,m*(jj-1)+1:m*jj));
         X=(lyap(At,Bt,Ct));
         S(m*(j-1)+1:m*j,m*(jj-1)+1:m*jj)=X;
         S(m*(jj-1)+1:m*jj,m*(j-1)+1:m*j)=X';
      end
   end
   time_newlyap=toc;
   %disp(['Newlyap solution took: ' num2str(time_newlyap) ' [sec]']);
end

% ----------------- EXACT-------------------- %
if exact,
   S=[];
   tic
   m=2;
   for j=1:n/2
      for jj=j:n/2
         At=(A(m*(j-1)+1:m*j,m*(j-1)+1:m*j));
         Bt=(A(m*(jj-1)+1:m*jj,m*(jj-1)+1:m*jj)');
         Ct=(BB(m*(j-1)+1:m*j,m*(jj-1)+1:m*jj));
         a11=At(1,1); a12=At(1,2); a21=At(2,1); a22=At(2,2);
         b11=Bt(1,1); b12=Bt(1,2); b21=Bt(2,1); b22=Bt(2,2);
         c11=Ct(1,1); c12=Ct(1,2); c21=Ct(2,1); c22=Ct(2,2);
         %x11=-(c11*b22*a22^2-a12*c21*b22^2+c11*a11*a22^2+c11*b11*b22^2-b21*c12*a22^2+c11*a11*a22*b22+c11*a11*b11*a22+c11*a11*b11*b22-c11*a11*b21*b12+c11*b22*b11*a22-c11*b22*b21*b12-c11*a21*a12*a22-c11*a21*a12*b11-b21*c12*a22*b22-b21*c12*b11*a22-b21*c12*b11*b22-b21*c12*a21*a12-a12*c21*b21*b12-a12*c21*a11*a22-a12*c21*a11*b22-a12*c21*a22*b22+a12*b21*c22*a11+a12*b21*c22*b11+a12*b21*c22*a22+a12*b21*c22*b22+c11*a22*b22^2+b21^2*c12*b12+a12^2*c21*a21)/(a11*b22*a22^2+a11^2*a22*b22+a11*b11*b22^2+b11*a22*b22^2-a21*a12*b22^2+2*a11*b22*b11*a22-a11*b22*b21*b12-2*a11*a21*a12*a22-a11*a21*a12*b11-b11*a11*b21*b12-2*b11*b22*b21*b12-b11*a21*a12*a22-b12*b21*a22*b22-b12*b21*b11*a22-2*b12*a21*a12*b21-a21*a12*a11*b22-a21*a12*a22*b22+a11^2*b11*a22+a11^2*b11*b22-a11^2*b21*b12+a11*b11^2*a22+a11*b11^2*b22+b22*b11^2*a22-a21*a12*b11^2+a11*a22*b22^2+b11*a11*a22^2+b11*b22*a22^2+a11^2*a22^2+b11^2*b22^2+b21^2*b12^2+a21^2*a12^2-b12*b21*a22^2);
         %x12=-(c12*b11*a22^2+c12*a11*a22^2-a12*c22*b11^2-c11*b12*a22^2+c11*b12^2*b21+c12*b11^2*a22+c12*b11^2*b22+a12^2*c22*a21-c11*b12*a22*b22-c11*b12*b11*a22-c11*b12*b11*b22-c11*b12*a21*a12+c12*a11*a22*b22+c12*a11*b11*a22+c12*a11*b11*b22-c12*a11*b21*b12+c12*b22*b11*a22-c12*b11*b21*b12-c12*a21*a12*a22-c12*a21*a12*b22+a12*c21*b12*a11+a12*c21*b12*b11+a12*c21*b12*a22+a12*c21*b12*b22-a12*c22*a11*a22-a12*c22*a11*b11-a12*c22*b11*a22-a12*c22*b21*b12)/(a11*b22*a22^2+a11^2*a22*b22+a11*b11*b22^2+b11*a22*b22^2-a21*a12*b22^2+2*a11*b22*b11*a22-a11*b22*b21*b12-2*a11*a21*a12*a22-a11*a21*a12*b11-b11*a11*b21*b12-2*b11*b22*b21*b12-b11*a21*a12*a22-b12*b21*a22*b22-b12*b21*b11*a22-2*b12*a21*a12*b21-a21*a12*a11*b22-a21*a12*a22*b22+a11^2*b11*a22+a11^2*b11*b22-a11^2*b21*b12+a11*b11^2*a22+a11*b11^2*b22+b22*b11^2*a22-a21*a12*b11^2+a11*a22*b22^2+b11*a11*a22^2+b11*b22*a22^2+a11^2*a22^2+b11^2*b22^2+b21^2*b12^2+a21^2*a12^2-b12*b21*a22^2);
         %x21=-(-b21*c22*a11^2+c21*b11*b22^2+c21*a11^2*a22+c21*a11^2*b22-c11*a21*b21*b12-c11*a21*a11*a22-c11*a21*a11*b22-c11*a21*a22*b22+b21*a21*c12*a11+b21*a21*c12*b11+b21*a21*c12*a22+b21*a21*c12*b22+c21*a11*a22*b22-c21*a11*a21*a12+c21*a11*b11*a22+c21*a11*b11*b22+c21*b22*b11*a22-c21*a21*a12*b11-c21*b21*b12*a22-c21*b22*b21*b12-b21*c22*a11*b22-b21*c22*a11*b11-b21*c22*b11*b22-b21*c22*a21*a12-c11*a21*b22^2+c11*a21^2*a12+c21*a11*b22^2+b21^2*c22*b12)/(a11*b22*a22^2+a11^2*a22*b22+a11*b11*b22^2+b11*a22*b22^2-a21*a12*b22^2+2*a11*b22*b11*a22-a11*b22*b21*b12-2*a11*a21*a12*a22-a11*a21*a12*b11-b11*a11*b21*b12-2*b11*b22*b21*b12-b11*a21*a12*a22-b12*b21*a22*b22-b12*b21*b11*a22-2*b12*a21*a12*b21-a21*a12*a11*b22-a21*a12*a22*b22+a11^2*b11*a22+a11^2*b11*b22-a11^2*b21*b12+a11*b11^2*a22+a11*b11^2*b22+b22*b11^2*a22-a21*a12*b11^2+a11*a22*b22^2+b11*a11*a22^2+b11*b22*a22^2+a11^2*a22^2+b11^2*b22^2+b21^2*b12^2+a21^2*a12^2-b12*b21*a22^2);
         %x22=-(c22*a11*b11^2+c22*a11^2*b11-c21*b12*a11^2-c12*a21*b11^2+c12*a21^2*a12+c21*b12^2*b21+c22*a11^2*a22+c22*b11^2*b22+c11*a21*b12*a11+c11*a21*b12*b11+c11*a21*b12*a22+c11*a21*b12*b22-c12*a21*a11*a22-c12*a21*a11*b11-c12*a21*b11*a22-c12*a21*b21*b12-c21*b12*a11*b22-c21*b12*a11*b11-c21*b12*b11*b22-c21*b12*a21*a12+c22*a11*a22*b22+c22*a11*b11*b22+c22*a11*b11*a22+c22*b22*b11*a22-c22*b21*b12*a22-c22*b11*b21*b12-c22*a11*a21*a12-c22*a21*a12*b22)/(a11*b22*a22^2+a11^2*a22*b22+a11*b11*b22^2+b11*a22*b22^2-a21*a12*b22^2+2*a11*b22*b11*a22-a11*b22*b21*b12-2*a11*a21*a12*a22-a11*a21*a12*b11-b11*a11*b21*b12-2*b11*b22*b21*b12-b11*a21*a12*a22-b12*b21*a22*b22-b12*b21*b11*a22-2*b12*a21*a12*b21-a21*a12*a11*b22-a21*a12*a22*b22+a11^2*b11*a22+a11^2*b11*b22-a11^2*b21*b12+a11*b11^2*a22+a11*b11^2*b22+b22*b11^2*a22-a21*a12*b11^2+a11*a22*b22^2+b11*a11*a22^2+b11*b22*a22^2+a11^2*a22^2+b11^2*b22^2+b21^2*b12^2+a21^2*a12^2-b12*b21*a22^2);
         %X=[x11 x12;x21 x22];
         Q=[a11+b11 b21 a12 0;b12 a11+b22 0 a12;a21 0 a22+b11 b21;0 a21 b12 a22+b22];
         XX=inv(Q)*[-c11;-c12;-c21;-c22];
         X=[XX(1,1) XX(2,1);XX(3,1) XX(4,1)];
         S(m*(j-1)+1:m*j,m*(jj-1)+1:m*jj)=X;
         S(m*(jj-1)+1:m*jj,m*(j-1)+1:m*j)=X';
      end
   end
   time_newexact=toc;
   %disp(['Newlyap solution took: ' num2str(time_newexact) ' [sec]']);
end


% Consider quality of grammians
% -----------------------------
res_new = A*S + S*A' + BB ;   %  The resultants SHOULD be zero.

max_res = max(max(abs(res_new))) ;
disp(['Maximum Lyapunov resultant: ',num2str(max_res)])