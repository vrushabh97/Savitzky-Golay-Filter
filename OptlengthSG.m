function[SNRIN,SNROP,y,OPT_M,MSE] = OptlengthSG(x,val,lp,M_Max,Ord)
% val is the ecg signal to which noise is to be added,to add noise the
% below steps are followed
% w = randn(1,L);
% w = (w-mean(w))/std(w);
% sigmas = sum(val.^2)/L;
% lp = sqrt(sigmas/realpow(10,0.1*3));
% x = (val + lp*w);
%Ord is the order to which the polynomial is approximated
L = length(val);
lambda = 1.2*lp*lp;%regularization parameter
%% construction of the filter folloowed by optimizing
M_Min = ceil(0.5*(Ord+1));%the minimum half window length
M_Vec= M_Min:M_Max';%%the maximum half window length
OPT_M = zeros(L,1);% the vector of optimum half winodw length for every point 
lambda = 1.2*lp*lp;%Regularization prameter
for i = 1:L
   count = 1;
   REG_RISK = zeros(length(M_Vec),1);
   for M = M_Min:M_Max
       x_n = [zeros(1,M),x,zeros(1,M)];%the vector constructed such that every element of x is a center value of an interval
       T = (-M:M)';
       A = zeros(length(T),Ord+1);
       WIND = x_n(1,i:2*M+i)';
       A(:,1) = 1;
       for j=1:Ord
           A(:,j+1)=T(:,1).^j;
       end
    H = pinv(A'*A)*A';
    RISK = (sum((A*H*WIND).^2) - 2*WIND'*H'*A'*WIND + 2*(lp^2)*sum(diag(transpose(A*H))) + sum((WIND).^2))/(2*M+1) - (lp^2);
    REG_RISK(count) = RISK + (lambda/(2*M+1))*sum(diag((A*H).^2));
    count = count+1;        
    [~,INDEX] = min(REG_RISK);
    OPT_M(i,1)=M_Vec(1,INDEX); 
   end
end
%% constucting the approximate clean signal with the aid of optimum half window length
y = zeros(1,length(x));
for k=1:L
    WIND = x_n(1,(M_Max-OPT_M(k,1)+k):(M_Max+OPT_M(k,1)+k))';
    d = (-OPT_M(k,1):OPT_M(k,1))';
    A = zeros(length(d),Ord+1);
    A(:,1)=1;
    for j=1:Ord
        A(:,j+1) = d(:,1).^j;
    end
    H = pinv(A'*A)*A';
    g = A*H*WIND;
    y(k) = g(OPT_M(k,1)+1);
end
%% Performance
noise1 = y - val;
SNROP = snr(y,noise1);
MSE = 10*log10((sum((x-y).^2)/length(y)));
subplot(311)
plot(x);axis([0 7200 -150 150]); xlabel('number of samples');ylabel('amplitude');
subplot(312)
plot(y);axis([0 7200 -150 150]); xlabel('number of samples');ylabel('amplitude');
subplot(313)
stem(OPT_M);axis([0 7200 0 30]);xlabel('number of samples');ylabel('Length');
end
