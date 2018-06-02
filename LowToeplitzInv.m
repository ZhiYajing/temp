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

        w=ToeplitzMatVec(t(n/2+1:n),t(1:n /2),tt);
        w=ToeplitzMatVec(-tt,[],w);
        x(n /2+1: n)=w;
        x=x(1:p);
    else x=1/t;
    end
end
