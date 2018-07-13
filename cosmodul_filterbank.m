function [y,Vm,H,F]=cosmodul_filterbank(h,f,M,x)
N = length(h); 
D=N-1; 
m1 = [0:M-1]; 
n1 = [0:N-1];  
arg = repmat((-1).^m1.'*pi/4,1,N);
theta1 = 2*cos(pi/M*(m1+0.5).'*(n1-D/2)+arg); 
theta2 = 2*cos(pi/M*(m1+0.5).'*(n1-D/2)-arg); 
H = repmat(h,M,1).*theta1;
F = repmat(f,M,1).*theta2; 
X_buf=zeros(M,N); % X_buf and W_buf have the states of filters H and G.
W_buf=zeros(M,N); 
Vm=[]; 
lx=length(x); K=ceil(lx/M); L=K*M;
for n=1:L+N-1
   if n<=lx   
       xn=repmat(x(n),M,1);  
   else
       xn=zeros(M,1);  
   end
   X_buf = [xn  X_buf(:,1:N-1)];
   Hx = sum(H.*X_buf,2); % Analysis filter bank outputs
   if mod(n-1,M)==0
     Vm=[Vm Hx]; Wm=Hx; % Downsampling
    else   
     Wm=zeros(M,1); % Upsampling
   end
   W_buf = [Wm  W_buf(:,1:N-1)];
   Fw = sum(F.*W_buf,2); % Synthesis filter bank outputs
   y(n) = sum(Fw);
end
end