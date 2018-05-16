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
        x(n /2+1: n)= -T* (C* (T*eye(length(T),1)));
        
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