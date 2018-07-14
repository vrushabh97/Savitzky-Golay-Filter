function [p,mu] = plofit(x,y,n)
%  x = input samples
%  y  = output function,n = order
m = length(x); %number of rows in the Vandermonde Matrix
V = zeros(m,n);
a = n;
for i = 1:m
    v = zeros(1,n);
    for j = a:-1:1
      v(n+1-j) =  realpow(x(i),j);
    end
    V(i,:) = v;
end
V(:,n+1)=ones(m,1);% adding 1 column to ones to the vandermonde matrix
%% QR method to compute the least squares solution for the coefficients,'p'
[Q,R] = qr(V,0);
p = transpose(R \ (transpose(Q) * y'));
f = polyval(p,x);
%% to find mean
mean = sum(x)/length(x);
sq = 0;
for i =1:length(x)
    sq = sq + (x(i)-mean)^2;
end
sd = (sq/length(x))^0.5;
mu = [mean;sd];
end
