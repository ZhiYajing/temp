function x=draft(L,S,b)
[n,n]=size(L);[m,m]=size(S)
% a=10*(rand(n,1));L=tril(toeplitz(a));
% S=toeplitz([2 -1 zeros(1, m-2)]);
% b=10*(rand(n*m,1));
% solve (kron(L,eye(m)) + kron(eye(n),S)) x = b

% S= Q'AQ, A=diag(eig(S));
a=eig(S);
A=diag(a);

%f=P'c;c=b_tilde;c=kron(eye(n),Q)b;f=P'kron(eye(n),Q)b;
%c=kron(eye(n),Q)b=vec(QBI)=vec(QB);
for k=1:n
   B(:,k)=b((k-1)*m+1 : k*m);
   BB(:,k)=sine_transform_data(B(:,k));
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
% what is f ?????
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
   YY(:,j)=sine_transform_data(Y(:,j));
end
x=YY(:);
%Q'Y=[Q'y1 Q'y2 ... Q'yn];
%Q'y(i)=fst(y(i));
%x=Q'Y(:);






