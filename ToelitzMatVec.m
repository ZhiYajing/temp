function x=ToelitzMatVec(T,y)
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

