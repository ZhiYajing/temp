function x= doit13 (a)
[n,n]= size (a ); x= zeros (n,n);
    if dividable (a )==1
        [b ,c ,d ,f ]= divide (a );
        x (1: n /2 ,1: n /2)= doit13 (b);
        x (1: n/2 ,n /2+1: n )= c;
        x(n /2+1: n ,1: n /2)= -1* doit13 (d )* f* doit13 (b );
        x(n /2+1: n ,n /2+1: n )= doit13 (d );
    else
        if (n ==1)
            x= inv (a );
        else
            disp('Hello World.');

        end
    end
end
function x= dividable (a)
[n,n ]= size (a);
x =((n >=2)&&( mod (n ,2)==0));
end
function [b ,c ,d ,f] = divide (a)
[n ,n ]= size (a );
b= zeros (n/2 ,n /2);
c=b;
d=b;
f=b;
for i =1: n /2
for j=n /2: n
b=a (1:i ,1: i);
f=a(n /2+1: j ,1: i );
d=a(n /2+1: j ,n /2+1: j );
c= zeros (n/2 ,n /2) ;
end
end
end