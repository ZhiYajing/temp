function x=ToeplitzMatVec(t1,t2,y)
%n=8
%t1=10*(rand(n,1));t2=[t1(1);10*(rand((n-1),1))];y=10*(rand(n,1));
%T=toeplitz(t1,t2); x=T*y;

n=length(t1);
m=length(t2);
%c=T(1,:); %1st row
if m==0
    t2=[t1(1);zeros((n-1),1)];
end
%r=t1 1st column

z=[y;zeros(n,1)];

%Cz=F*AFz;lamda=fft(C(:,1));A=diag(lamda);
lamda=fft([t1;0;t2(n:-1:2)]);
x=fft(z);
x=lamda.*x;
x=ifft(x);
x=x(1:n);
end
