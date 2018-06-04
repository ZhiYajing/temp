%%%%%%%%%%%%%%%%%%%%%%test eigenvalue
 a1=eig(S1);a2=eig(S2);
 m=6;Q=zeros(m);
 for i=1:m
    for j=1:m
       Q(i,j)=sqrt(2/(m+1))*sin(pi*i*j/(m+1)); 
    end
 end
 for i=1:m 
     a1(i)=2+2*cos(i*pi/(m+1));
 end
 A1=diag(a1);
 for i=1:m
     a2(m-i+1)=2+2*cos(i*pi/(m+1));
% end
 A2=diag(a2);
 Q'*A1*Q;
 Q'*A2*Q;
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%vector
 B=reshape(b,[m,n]);
for i=1:n
   B1= reshape(B(:,i),[m2,m1]); % B1= reshape(bi,[m2,m1]);
   B1=10*B1;
   B1=B1';%E=E';
   B1=10*B1;
   B1=B1';%G=G';
   BB(:,i)=B1(:);%B(:,i)=G(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
   for j=1:m1
       bb((j-1)*m2+1:j*m2)=10*b(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2);
   end
   for k=1:m
       if  rem(k,m1)==0
           cc(k)=bb((m1-1)*m2+ceil(k/m1));
       else cc(k)=bb((rem(k,m1)-1)*m2+ceil(k/m1)); 
       end 
   end

   for j=1:m2
       bb((j-1)*m1+1:j*m1)=10*cc(((j-1)*m1+1):j*m1); 
   end
   for k=1:m
       if  rem(k,m2)==0
           cc(k)=bb((m2-1)*m1+ceil(k/m2));
       else cc(k)=bb((rem(k,m2)-1)*m1+ceil(k/m2)); 
       end 
   end
   a((i-1)*m+1 : i*m)=cc;
end
for i=1:m*n
       if  rem(i,n)==0
           f(i)=a((n-1)*m+ceil(i/n));
       else f(i)=a((rem(i,n)-1)*m+ceil(i/n)); 
       end 
end
f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%vector
B=reshape(b,[m,n]);
for i=1:n
   %B(:,i)=b((i-1)*m+1 : i*m);
   B1= reshape(B(:,i),[m2,m1]); % B1= reshape(bi,[m2,m1]);
   %E=Q2*B1
   %************
   B1=dst(B1)*sqrt(2/(m2+1));%E=dst(B1)*sqrt(2/(m2+1));
   %[b(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2)]
   %G=E*Q1'=(Q1*E')'
   %*************
   for j=1:m1
       bb((j-1)*m2+1:j*m2)=dst(b(((i-1)*m+(j-1)*m2+1):(i-1)*m+j*m2))*sqrt(2/(m2+1)); 
   end
   B1=B1';%E=E';
   for k=1:m
       if  rem(k,m1)==0
           cc(k)=bb((m1-1)*m2+ceil(k/m1));
       else cc(k)=bb((rem(k,m1)-1)*m2+ceil(k/m1)); 
       end 
   end
   %*********
   B1=dst(B1)*sqrt(2/(m1+1));%G=dst(E)*sqrt(2/(m1+1));
   %**********
   for j=1:m2
       bb((j-1)*m1+1:j*m1)=dst(cc(((j-1)*m1+1):+j*m1))*sqrt(2/(m1+1)); 
   end
   B1=B1';%G=G';
   
   for k=1:m
       if  rem(k,m2)==0
           cc(k)=bb((m2-1)*m1+ceil(k/m2));
       else cc(k)=bb((rem(k,m2)-1)*m1+ceil(k/m2)); 
       end 
   end
   a((i-1)*m+1 : i*m)=cc;
   B(:,i)=B1(:);%B(:,i)=G(:);
end
