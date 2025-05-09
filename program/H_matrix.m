function [H]=H_matrix(log_ampx)
% n: Number of geophones
% m: Number of frequency samples
% lamda:Absorption Constraint Weight Coefficient
[m,n]=size(log_ampx);
%%
L1=zeros(m-1,m);
for i=1:m-1
    L1(i,i)=-1;
    L1(i,i+1)=1;
end

%%
y=reshape(log_ampx,m*n,1);
%%

N=n-1;
Pi=cell(1,N+1);
for i=1:N
    Pi{i}=zeros(m,m*n);
    for j=1:m
        Pi{i}(j,(i-1)*m+j)=-1;
        Pi{i}(j,i*m+j)=1;% 
    end
end
H=zeros(m-1,N);
for i=1:N
    Gii=L1*Pi{i};
    H(:,i)=Gii*y;
    
end
H=H';
% H=(H.^2)';