function x=testToepSys(l,s1,s2,b)
n=length(l);m1=s1(4);m2=s2(4);m=m1*m2;
% INPUT:
% n=9;m1=5;m2=4;m=m1*m2;
% l=10*(rand(n,1));L=tril(toeplitz(l));
% s1=[-1;2;-1;m1];s2=[-1;2;-1;m2];
% S1=toeplitz([2 -1 zeros(1, m1-2)]);
% S2=toeplitz([2 -1 zeros(1, m2-2)]);
% b=10*(rand(n*m,1));
% solve:
% (kron(L,eye(m)) + kron(eye(n),(kron(eye(m1),S2) + kron(S1,eye(m2))))) x = b
% xx=inv(kron(L,eye(m)) + kron(eye(n),(kron(eye(m1),S2) + kron(S1,eye(m2)))))*b

% S1= Q1'*A1*Q1, A1=diag(eig(S1));
% S2= Q2'*A2*Q2, A2=diag(eig(S2));
% S_hat=(kron(Q1',Q2'))*(kron(eye(m1),A2)+kron(A1,eye(m2)))*(kron(Q1,Q2))
%      = Q_hat'*A_hat*Q_hat;

% a1=eig(S1);a2=eig(S2);
for i=1:m1
    a1(m1-i+1)=2+2*cos(i*pi/(m1+1));
end
for i=1:m2
    a2(m2-i+1)=2+2*cos(i*pi/(m2+1));
end

% f=P'c;c=b_tilde;c=kron(eye(n),Q_hat)b;Q_hat=kron(Q1,Q2);f=P'kron(eye(n),Q_hat)b;
% c=kron(eye(n),Q_hat)b=vec(Q_hat*BI)=vec(Q_hat*B); B is m*n;

% Q_hat*B=kron(Q1,Q2)*[b1 b2 ... bn]; bi=B(:,i);
% kron(Q1,Q2)*bi=vec(Q2*B1*Q1'); B1= reshape(bi,[m2,m1]);
% E=Q2*B1=[dst(bi1) dst(bi2) ... dst(bim1)];
% E*Q1'=(Q1*E')'=(Q1*F)'=[dst(f1) dst(f2) ... dst(fm2)]';F=E'=[f1 f2 fm2];
B=reshape(b,[m,n]);
for i=1:n
   %B(:,i)=b((i-1)*m+1 : i*m);
   %E=Q2*B1
   %************
   %[b(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2)]
   %G=E*Q1'=(Q1*E')'
   %*************
   for j=1:m1
       %bb((j-1)*m2+1:j*m2)=sine_transform_data(b(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2));
       bb((j-1)*m2+1:j*m2)=sqrt(2/(m2+1))*dst(b(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2));
   end
   for k=1:m
       if  rem(k,m1)==0
           cc(k)=bb((m1-1)*m2+ceil(k/m1));
       else cc(k)=bb((rem(k,m1)-1)*m2+ceil(k/m1)); 
       end 
   end
   %*********
   %**********
   for j=1:m2
       %bb((j-1)*m1+1:j*m1)=sine_transform_data(cc(((j-1)*m1+1):+j*m1)); 
       bb((j-1)*m1+1:j*m1)=sqrt(2/(m1+1))*dst(cc(((j-1)*m1+1):+j*m1)); 
   end
  
   for k=1:m
       if  rem(k,m2)==0
           cc(k)=bb((m2-1)*m1+ceil(k/m2));
       else cc(k)=bb((rem(k,m2)-1)*m1+ceil(k/m2)); 
       end 
   end
   a((i-1)*m+1 : i*m)=cc;
end
a;
%c=B(:); %c=b_tilde;
%QB=[Qb1 Qb2 ... Qbn];
%Qb(i)=fst(b(i));
%c=QB(:);

%f=P'c=vec(C'); f=b_hat;
% for i=1:n
%     C(:,i)=c((i-1)*m+1 : i*m);
% end

for i=1:m*n
       if  rem(i,n)==0
           f(i)=a((n-1)*m+ceil(i/n));
       else f(i)=a((rem(i,n)-1)*m+ceil(i/n)); 
       end 
end
%%%%%%
f;
% f=b_hat;

% ( kron(eye(m),L)+ kron(A, eye(n)) ) z = f
% x_hat=z; b_hat=f;
% solve (L+ lamda*I)x'=b'; length(a)=m;
for i=1:m
    if rem(i,m2)==0
        lamda=a2(m2)+a1(ceil(i/m2));
    else lamda=a2(rem(i,m2))+a1(ceil(i/m2));
    end
    t=LowToeplitzInv([l(1)+lamda;l(2:length(l))]);
    y=f((i-1)*n+1 : i*n);
%******************    
    bb((i-1)*n+1:i*n)= ToeplitzMatVec(t,[],y);
    %Z(:,i)=LowToeplitzInv(L+ lamda*eye(n)) * f((i-1)*n+1 : i*n);
end 
%z=Z(:);%z=x_hat
%z=P'y; y=x_tilde;y=kron(eye(n),Q_hat)x; 
%x=kron(eye(n),Q_hat')Pz;
%y=Pz=vec(Z');
%*************

for i=1:m*n
       if  rem(i,m)==0
           y(i)=bb((m-1)*n+ceil(i/m));
       else y(i)=bb((rem(i,m)-1)*n+ceil(i/m)); 
       end 
end
y;

% x=kron(eye(n),Q_hat')y=vec(Q_hat'YI)=vec(Q_hat'Y);Y=reshape(y,[m,n]);
% Q_hat'Y=kron(Q1',Q2')*[y1 y2 ... yn]; yi=Y(:,i);
% kron(Q1',Q2')*yi=Q2'*Y1*Q1; Y1=reshape(yi,[m2,m1]);
% E=Q2'*Y1=[idst(yi1) idst(yi2) ... idst(yim1)];
% E*Q1=(Q1'*E')'=(Q1'*F)'=[idst(f1) idst(f2) ... idst(fm2)]';F=E'=[f1 f2 fm2];
% E*Q1'=(Q1*E')'=(Q1*F)'=[dst(f1) dst(f2) ... dst(fm2)]';F=E'=[f1 f2 fm2];
%Y=reshape(y,[m,n]);
for i=1:n
  %E=Q2'*Y1
  %G=E*Q1=(Q1'*E')';
  %*************
   for j=1:m1
       bb((j-1)*m2+1:j*m2)=sqrt(2/(m2+1))*dst(y(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2)); 
       %bb((j-1)*m2+1:j*m2)=sine_transform_data(y(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2)); 
   end
  for k=1:m
       if  rem(k,m1)==0
           cc(k)=bb((m1-1)*m2+ceil(k/m1));
       else cc(k)=bb((rem(k,m1)-1)*m2+ceil(k/m1)); 
       end 
  end
   %*********
  %**********
   for j=1:m2
       bb((j-1)*m1+1:j*m1)=sqrt(2/(m1+1))*dst(cc(((j-1)*m1+1):+j*m1));
       %bb((j-1)*m1+1:j*m1)=sine_transform_data(cc(((j-1)*m1+1):+j*m1)); 
   end
  
   for k=1:m
       if  rem(k,m2)==0
           cc(k)=bb((m2-1)*m1+ceil(k/m2));
       else cc(k)=bb((rem(k,m2)-1)*m1+ceil(k/m2)); 
       end 
   end
   x((i-1)*m+1 : i*m)=cc;
end
%x=a;
%Q'Y=[Q'y1 Q'y2 ... Q'yn];
%Q'y(i)=fst(y(i));
%x=Q'Y(:);
end

function x= LowToeplitzInv(t)
% t=10*(rand(9,1));A=toeplitz(t);A=tril(A);
n=length(t);x=zeros(n,1);     
     if n>=2
         m=pow2(ceil(log2(n)));p=n;
        if m~= n
            t=[t;zeros(m-n,1)];
            n=m;
        end
        tt=LowToeplitzInv(t(1:n /2));
%         X (1: n /2 ,1: n /2)= T;
%         X (1: n/2 ,n /2+1: n )= 0;
%         X(n /2+1: n ,1: n /2)= -1* T* C* T;
%         X(n /2+1: n ,n /2+1: n )= T;
        
        x(1: n /2)=tt;
%         x(n /2+1: n)= -T* (C* (T*eye(length(T),1)));

        w=ToeplitzMatVec(t(n/2+1:n),t(n/2+1:-1:2),tt);
        w=ToeplitzMatVec(-tt,[],w);
        x(n /2+1: n)=w;
        x=x(1:p);
    else x=1/t;
    end
end

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

z=[y(:);zeros(n,1)];

%Cz=F*AFz;lamda=fft(C(:,1));A=diag(lamda);
lamda=fft([t1;0;t2(n:-1:2)]);
x=fft(z);
x=lamda.*x;
x=ifft(x);
x=x(1:n);
end
