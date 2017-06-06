function [x,U]=gausselim(a,b)
% x = gausselim(a,b) computes the solution of an nxn matrix using gaussian
% elimination
%
% inputs
% ------
% a: an n x n matrix
% b: an n x 1 column vector, equivalent to 'y' in y = a * x 
%
% outputs
% -------
% x: the n x 1 solution to the system a*x = b

n=length(b);
m=zeros(n,1);
x=zeros(n,1);

for k =1:n-1
      m(k+1:n) = a(k+1:n,k)/a(k,k);
      for i=k+1:n
          a(i, k+1:n) = a(i,k+1:n)-m(i)*a(k,k+1:n);
      end;
      b(k+1:n)=b(k+1:n)-b(k)*m(k+1:n);
end
U= triu(a);
x(n)=b(n)/a(n,n);
for k =n-1:-1:1
      b(1:k)=b(1:k)-x(k+1)* U(1:k,k+1);
      x(k)=b(k)/U(k,k);
end
end