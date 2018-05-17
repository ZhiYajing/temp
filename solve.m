function x=solve(L,S,b)
[n,n]=size(L);[m,m]=size(S)
% n=32;m=8;
% a=10*(rand(n,1));L=tril(toeplitz(a));
% S=toeplitz([2 -1 zeros(1, m-2)]);
% b=10*(rand(n*m,1));
% solve (kron(L,eye(m)) + kron(eye(n),S)) x = b

% S= Q'AQ, A=diag(eig(S));
a=eig(S);
A=diag(a);

% f=P'c;c=b_tilde;c=kron(eye(n),Q)b;f=P'kron(eye(n),Q)b;
% c=kron(eye(n),Q)b=vec(QBI)=vec(QB);
for k=1:n
   B(:,k)=b((k-1)*m+1 : k*m);
   BB(:,k)=dst(B(:,k))*sqrt(2/(m+1));
end
c=BB(:);
%QB=[Qb1 Qb2 ... Qbn];
%Qb(i)=fst(b(i));
%c=QB(:);

%f=P'c=vec(C');
for h=1:n
    C(:,h)=c((h-1)*m+1 : h*m);
end
D=C';
f=D(:);

% ( kron(eye(m),L)+ kron(A, eye(n)) ) z = f
% x_hat=z; b_hat=f;
% solve (L+ lamda*I)x'=b'; length(a)=m;
for i=1:m
    lamda=a(i);
    Z(:,i)=LowToeplitzInv(L+ lamda*eye(n)) * f((i-1)*n+1 : i*n);
end 
z=Z(:);

%z=P'y; y=x_tilde;y=kron(eye(n),Q)x; 
%x=kron(eye(n),Q')Pz;
%y=Pz=vec(Z');
W=Z';
y=W(:);

%x=kron(eye(n),Q')y=vec(Q'YI)=vec(Q'Y);
for j=1:n
   Y(:,j)=y((j-1)*m+1 : j*m); 
   YY(:,j)=dst(Y(:,j))*sqrt(2/(m+1));
end
x=YY(:);
%Q'Y=[Q'y1 Q'y2 ... Q'yn];
%Q'y(i)=fst(y(i));
%x=Q'Y(:);
end

function X= LowToeplitzInv(A)
% a=10*(rand(8,1));A=toeplitz(a);A=tril(A);

[n,n]= size (A ); X= zeros (n,n);
    if dividable (A)==1
        [B ,C]= divide (A );
        T=LowToeplitzInv (B);
         
%         X (1: n /2 ,1: n /2)= T;
%         X (1: n/2 ,n /2+1: n )= 0;
%         X(n /2+1: n ,1: n /2)= -1* T* C* T;
%         X(n /2+1: n ,n /2+1: n )= T;
        
        x(1: n /2)=T*eye(length(T),1);
%         x(n /2+1: n)= -T* (C* (T*eye(length(T),1)));
        w=ToeplitzMatVec(C,T(:,1));
        w=ToeplitzMatVec(-T,w);
        x(n /2+1: n)=w;
        
        X (1: n /2 ,1: n /2)= T;
        X (1: n/2 ,n /2+1: n )= 0;
        X(n /2+1: n ,1: n /2)= toeplitz(x(n/2+1 : n),x(n/2+1:-1:2));
        X(n /2+1: n ,n /2+1: n )= T;
    else X=1/A;
    end
end
function X= dividable (A)
[n,n ]= size (A);
X =((n >=2)&&( mod (n ,2)==0));
end
function [B ,C] = divide (A)
[n ,n ]= size (A );

B=A (1:n /2 ,1: n /2);
C=A(n /2+1: n ,1: n /2 );

end

function x=ToeplitzMatVec(T,y)
%n=8
%a=10*(rand(n,1));b=[a(1);10*(rand((n-1),1))];T=toeplitz(a,b);
%y=10*(rand(n,1));
[n,n]=size(T);
t1=T(1,:); %1st row
t2=T(:,1); %1st column
c=[0 t1(n:-1:2)];
r=[0;t2(n:-1:2)];
B=toeplitz(c,r);
C=[T B;B T]; % circulant matrices
z=[y;zeros(n,1)];

%Cz=F*AFz
lamda=fft(C(:,1));
%A=diag(lamda);
x1=fft(z);
x2=lamda.*x1;
x3=ifft(x2);
x=x3(1:n);
end







