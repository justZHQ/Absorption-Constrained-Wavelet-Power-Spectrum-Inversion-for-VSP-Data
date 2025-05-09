function [log_ampx_Correction,ampx_Correction,Gaussian]=new_Fast(sigma,lambda,average,GausiannVariance,log_ampx)
% n: Number of geophones
% m: Number of frequency samples
% lambda:Spectrum Constraint (SC) Weight Coefficient. The smaller the lambda, the greater the weight of the Absorption Constraint (AC) term
% sigma:Variance of Gaussian function used for lambda weighting
[m,n]=size(log_ampx);
%%
L1=zeros(m-1,m);
for i=1:m-1
    L1(i,i)=-1;
    L1(i,i+1)=1;
end

L2=zeros(m-2,m);
for i=1:m-2
    L2(i,i)=1;
    L2(i,i+1)=-2;
    L2(i,i+2)=1;
end
%%
I=eye(m*n);
y=reshape(log_ampx,m*n,1);

%%
Gall=zeros(m*n,m*n);
K=1;
N=n-K;

% lamda=eye(n)*lamda;
if sigma==1
mid=round(average*n);
T=mid-n:mid-1;
Gaussian=exp(-T.^2/(2*GausiannVariance^2));
k=5;
s1=10;s2=5;
k1=4;k2=1;
% s1=10;s2=5;
% k1=1.5;k2=1;
x1=T([1,s1:s1+k]);
y1=[k1,Gaussian(s1:s1+k)];
x2=T([end-s2-k:end-s2,end]);
y2=[Gaussian(end-s2-k:end-s2),k2];
Gaussian(1:s1)=interp1(x1,y1,T(1:s1),'spline');
Gaussian(end-s2:end)=interp1(x2,y2,T(end-s2:end),'spline');
elseif sigma==0
%     lamda=lamda;
    Gaussian=ones(1,n);
end
lambda=lambda*(Gaussian);


dig=[];
for i=1:n
    d=ones(m,1)*lambda(i);
    dig=[dig;d];
end
dig=diag(dig);


for i=1:N
    Pi=zeros(m,m*n);
    for j=1:m
        Pi(j,(i-1)*m+j)=-1;
        Pi(j,(i-1+K)*m+j)=1;
    end
    Gii=L1*Pi;
%     Gi=L2*Pi;
Gall=Gall+(Gii'*Gii);
%     Gall=Gall+lamda1*(Gi'*Gi) +lamda2*(Gii'*Gii);   
end

A=[];
for i=1:n-1
    Pi=zeros(m,m*n);
    for j=1:m
        Pi(j,(i-1)*m+j)=-1;
        Pi(j,i*m+j)=1;
    end
    A=[A;L2*Pi];
end

% PTP=zeros(m*n,m*n);
% for i=1:N
%     Pi=zeros(m,m*n);
%     for j=1:m
%         Pi(j,(i-1)*m+j)=-1;
%         Pi(j,(i-1+K)*m+j)=1;
%     end
% PPP=Pi'*Pi;
% PTP=PPP+PTP;
% end
% y=PTP*y;
y=dig*y;
invG=inv(Gall+dig*I);
z=inv(A*invG*A')*A*invG*y;
log_ampx_Correction=invG*(y-A'*z);
    
log_ampx_Correction=reshape(log_ampx_Correction,m,n);
ampx_Correction=[];
for i=1:n
amp=exp(log_ampx_Correction(:,i));
ampx_Correction=[ampx_Correction,amp];
end