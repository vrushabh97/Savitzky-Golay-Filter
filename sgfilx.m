function[y,MSE,Bias,Var,SNRIN,SNROP] = sgflix(x,val,M,Ord)
tic
% loading of ecg signal and addition of WG noise
% x = (awgn(val,SNR,'measured'));x is the noisy signal
% M is the half window length
% Ord is the order to whcih the polynomial can be approximated
 L=length(x);
%% construction of Savitzky-Golay Filter
%M = 30;   %in samples
%Ord = 1; % order of the filter
T = (-M:1:M)';% column vector
x_n = [zeros(1,M),x,zeros(1,M)];
A = zeros(length(T),Ord+1);
A(:,1) = 1;
for j=1:Ord
    A(:,j+1)=T(:,1).^j;
end
H = pinv(A'*A)*A';
a = zeros(Ord+1,L);
for i =1:L 
    WIND = x_n(1,i:2*M+i)';
    a(:,i) = H*WIND;
end
y = A*a;
y = y(M+1,:);% considering M+1 row
%% to calculate mean square error
MSE = 10*log10(sum((val-y).^2)/L); % x is the signal with AWGN , polvalues is the recovered signal 
Bias = sum(y-val)/length(val);
Var = sum(y.^2)/length(y)-(sum(y)/length(y))^2;
noise1 = y - val;
SNROP = snr(x,noise1);
figure(1)
plot(val);%ylabel('Amplitude'),xlabel('Number of samples');%title('Clean signal');
%axis([0 7200 -100 100])
% figure(2)
% rishabh_spectro('FALK0_te.wav',0.04,0.025,1024);
% figure(3)
% rishabh_spectro('xyz.wav',0.04,0.025,1024);
figure(2)
plot(x); %ylabel('Amplitude'),xlabel('Number of samples');%title('Noisy signal');%axis([0 2000 -100 100])
%axis([0 7200 -150 150])
figure(3)
plot(y); %ylabel('Amplitude'),xlabel('Number of samples');%title('Recovered signal');%axis([0 2000 -100 100])
%xis([0 7200 -100 100])
toc
end

