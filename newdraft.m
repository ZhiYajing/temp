function x=newdraft(L,S1,S2,b)
[n,n]=size(L);[m1,m1]=size(S1);[m2,m2]=size(S2);m=m1*m2;
% n=4;m1=2;m2=4;m=m1*m2;
% a=10*(rand(n,1));L=tril(toeplitz(a));
% S1=toeplitz([2 -1 zeros(1, m1-2)]);
% S2=toeplitz([2 -1 zeros(1, m2-2)]);
% b=10*(rand(n*m,1));
% solve: S=(kron(eye(m1),S2) + kron(S1,eye(m2)));
% (kron(L,eye(m)) + kron(eye(n),(kron(eye(m1),S2) + kron(S1,eye(m2))))) x = b
% xx=inv(kron(L,eye(m)) + kron(eye(n),(kron(eye(m1),S2) + kron(S1,eye(m2)))))*b

% S1= Q1'*A1*Q1, A1=diag(eig(S1));
% S2= Q2'*A2*Q2, A2=diag(eig(S2));
% S_hat=(kron(Q1',Q2'))*(kron(eye(m1),A2)+kron(A1,eye(m2)))*(kron(Q1,Q2))
%      = Q_hat'*A_hat*Q_hat;

%a1=eig(S1);a2=eig(S2);
for i=1:m1
    a1(m1-i+1)=2+2*cos(i*pi/(m1+1));
end
for i=1:m2
    a2(m2-i+1)=2+2*cos(i*pi/(m2+1));
end
for i = 1:m1
   A(:,i)= a2+a1(i);
end
a=A(:);
% A1=diag(a1); A2=diag(a2);


% f=P'c;c=b_tilde;c=kron(eye(n),Q_hat)b;Q_hat=kron(Q1,Q2);f=P'kron(eye(n),Q_hat)b;
% c=kron(eye(n),Q_hat)b=vec(Q_hat*BI)=vec(Q_hat*B); B is m*n;

% Q_hat*B=kron(Q1,Q2)*[b1 b2 ... bn]; bi=B(:,i);
% kron(Q1,Q2)*bi=vec(Q2*B1*Q1'); B1= reshape(bi,[m2,m1]);
% E=Q2*B1=[dst(bi1) dst(bi2) ... dst(bim1)];
% E*Q1'=(Q1*E')'=(Q1*F)'=[dst(f1) dst(f2) ... dst(fm2)]';F=E'=[f1 f2 fm2];
B=reshape(b,[m,n]);
for i=1:n
   %B(:,i)=b((i-1)*m+1 : i*m);
   B1= reshape(B(:,i),[m2,m1]); % B1= reshape(bi,[m2,m1]);
   %E=Q2*B1
   for j=1:m1
      E(:,j)= sine_transform_data(B1(:,j));
   end
   %H=E*Q1'=(Q1*E')'
   F=E';
   for k=1:m2
      G(:,k)= sine_transform_data(F(:,k));
   end
   H=G';
   BB(:,i)=H(:);
end

c=BB(:); %c=b_tilde;
%QB=[Qb1 Qb2 ... Qbn];
%Qb(i)=fst(b(i));
%c=QB(:);

%f=P'c=vec(C'); f=b_hat;
C=reshape(c,[m,n]);
% for i=1:n
%     C(:,i)=c((i-1)*m+1 : i*m);
% end
C=C';
f=C(:);
% f=b_hat;

% ( kron(eye(m),L)+ kron(A, eye(n)) ) z = f
% x_hat=z; b_hat=f;
% solve (L+ lamda*I)x'=b'; length(a)=m;
for i=1:m
    lamda=a(i);
    T=LowToeplitzInv(L+ lamda*eye(n));
    y=f((i-1)*n+1 : i*n);
    Z(:,i)= ToeplitzMatVec(T,y);
    %Z(:,i)=LowToeplitzInv(L+ lamda*eye(n)) * f((i-1)*n+1 : i*n);
end 
z=Z(:);%z=x_hat

%z=P'y; y=x_tilde;y=kron(eye(n),Q_hat)x; 
%x=kron(eye(n),Q_hat')Pz;
%y=Pz=vec(Z');
Z=Z';
y=Z(:);

% x=kron(eye(n),Q_hat')y=vec(Q_hat'YI)=vec(Q_hat'Y);Y=reshape(y,[m,n]);
% Q_hat'Y=kron(Q1',Q2')*[y1 y2 ... yn]; yi=Y(:,i);
% kron(Q1',Q2')*yi=Q2'*Y1*Q1; Y1=reshape(yi,[m2,m1]);
% E=Q2'*Y1=[idst(yi1) idst(yi2) ... idst(yim1)];
% E*Q1=(Q1'*E')'=(Q1'*F)'=[idst(f1) idst(f2) ... idst(fm2)]';F=E'=[f1 f2 fm2];
% E*Q1'=(Q1*E')'=(Q1*F)'=[dst(f1) dst(f2) ... dst(fm2)]';F=E'=[f1 f2 fm2];
Y=reshape(y,[m,n]);
for i=1:n
  Y1=reshape(Y(:,i),[m2,m1]); %Y1=reshape(yi,[m2,m1]);
  %E=Q2'*Y1
  for j=1:m1
     E(:,j)= sine_transform_data(Y1(:,j));
  end
  %H=E*Q1=(Q1'*E')'
  F=E';
  for k=1:m2
     G(:,k)= sine_transform_data(F(:,k));
  end
  H=G';
  YY(:,i)=H(:);
end
x=YY(:);
%Q'Y=[Q'y1 Q'y2 ... Q'yn];
%Q'y(i)=fst(y(i));
%x=Q'Y(:);
end

function X= LowToeplitzInv(A)
% a=10*(rand(9,1));A=toeplitz(a);A=tril(A);

[n,n]= size (A ); X= zeros (n,n);
     %if dividable (A)==1
      
     if n>=2
         m=pow2(ceil(log2(n)));p=n;
        if m~= n
            a=[A(:,1);zeros(m-n,1)];
            A=toeplitz(a);A=tril(A);
            n=m;
        end
        [B ,C]= divide (A);
        T=LowToeplitzInv (B);
%         X (1: n /2 ,1: n /2)= T;
%         X (1: n/2 ,n /2+1: n )= 0;
%         X(n /2+1: n ,1: n /2)= -1* T* C* T;
%         X(n /2+1: n ,n /2+1: n )= T;
        
        x(1: n /2)=T*eye(length(T),1);
        %x(1: n /2)=T(:,1);
%         x(n /2+1: n)= -T* (C* (T*eye(length(T),1)));
        w=ToeplitzMatVec(C,T(:,1));
        w=ToeplitzMatVec(-T,w);
        x(n /2+1: n)=w;
        
        X (1: n /2 ,1: n /2)= T;
        X (1: n/2 ,n /2+1: p )= 0;
        %X(n /2+1: n ,1: n /2)= toeplitz(x(n/2+1 : n),x(n/2+1:-1:2));
        Y=toeplitz(x(n/2+1 : n),x(n/2+1:-1:2));
        X(n /2+1: p ,1: n /2)= Y(1:(p-n/2),:);
        %X(n /2+1: n ,n /2+1: n )= T;
        X(n /2+1: p ,n /2+1: p )= T(1:(p-n/2),1:(p-n/2));
    else X=1/A;
    end
end

function [B ,C] = divide (A)
[n ,n ]= size (A );

B=A(1:n /2 ,1: n /2);
C=A(n /2+1: n ,1: n /2 );
end

function x=ToeplitzMatVec(T,y)
%n=8
%a=10*(rand(n,1));b=[a(1);10*(rand((n-1),1))];T=toeplitz(a,b);
%y=10*(rand(n,1));
[n,n]=size(T);
c=T(1,:); %1st row
r=T(:,1); %1st column
c=[0 c(n:-1:2)];
r=[0;r(n:-1:2)];
B=toeplitz(c,r);
C=[T B;B T]; % circulant matrices
z=[y;zeros(n,1)];

%Cz=F*AFz
lamda=fft(C(:,1));
%A=diag(lamda);
x=fft(z);
x=lamda.*x;
x=ifft(x);
x=x(1:n);
end






