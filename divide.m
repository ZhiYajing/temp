function [B ,C] = divide (A)
[n ,n ]= size (A );
B= zeros (n/2 ,n /2);
C=B;
   B=A (1:n /2 ,1: n /2);
C=A(n /2+1: n ,1: n /2 );
end