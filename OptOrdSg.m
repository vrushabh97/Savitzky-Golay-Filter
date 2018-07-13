function[y,OPT_ORD,MSE,SNRIN,SNROP] = OptOrdSg(x,val,lp,Ord_Max,M)
tic
% val is the ecg signal to which noise is to be added,to add noise the
% below steps are followed
% w = randn(1,L);
% w = (w-mean(w))/std(w);
% sigmas = sum(val.^2)/L;
% lp = sqrt(sigmas/realpow(10,0.1*3));
% x = (val + lp*w);
L = length(val);
lambda = 1.2*lp*lp;%regularization parameter
%% construction of Savitzky-Golay Filter
%M = 5;
x_n = [zeros(1,M),x,zeros(1,M)];%the vector constructed such that every element of x is a center value of an interval
Ord_Min = 1;% minimum order
%Ord_Max = 7;%maximum order
T = (-M:M)'; 
A = zeros(2*M+1,Ord_Max+1);
A(:,1) = 1;
REG_RISK = zeros(1,Ord_Max);
for j=1:Ord_Max
    A(:,j+1)=T(:,1).^j;
end
Ord_Vec = Ord_Min:Ord_Max;
OPT_ORD = zeros(L,1);
for i=1:L
    WIND = x_n(1,i:2*M+i);
    count = 1;
    for Ord = Ord_Min:Ord_Max
        a = A(:,1:Ord+1);
        H = pinv(a'*a)*a';
        H2 = (a*H);
        H3 = a*H*WIND';
        div = 1/(2*M+1);
        RISK = div*(sum(H3.^2) - 2*H3'*WIND' + 2*(lp^2)*sum(diag(H2')) + sum((WIND).^2)) - (lp^2);
        REG_RISK(count) = RISK + lambda*sum(diag(H2.^2))/(2*M+1);
        count = count+1;        
    end
    [~,INDEX] = min(REG_RISK);%caluculating the order at which risk is minimum
    OPT_ORD(i)=Ord_Vec(1,INDEX);
end
y = reconstruct(x,OPT_ORD,M);%this fucntion constructs the signal element by element
MSE = 10*log10((sum((val-y).^2)/length(y)));%Mean squared error
noise1 = y' - val;
SNROP = snr(y',noise1);%SNR at the output
toc
end


