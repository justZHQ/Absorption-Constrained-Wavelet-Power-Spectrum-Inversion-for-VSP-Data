function [DATA,DATA_noise,dh,Time]=SyntheticData(wavelet_type,snr,dt,Nt)
%% Problem setting 
t=dt:dt:(Nt)*dt;
fsample=1/dt;   % number of samples
if wavelet_type=='Ricker'
    fp=50;          % peak frequency of ricker in hertz
    t0=1.6/fp;      % timeshift of ricker
    Srct=(1-0.5*(2*pi*fp*(t-t0)).^2).*exp(-0.25*(2*pi*fp*(t-t0)).^2); %雷克子波
elseif wavelet_type=='Ormsby'
    F_L1=20;
    F_L2=30;
    F_H1=60;
    F_H2=70;
    t0=4/45;
    t=t-t0+0.0001;
    s1=(((pi*F_H2).^2/(pi*F_H2-pi*F_H1)).*((sin(pi*F_H2.*t)./(pi*F_H2.*t)).^2)-((pi*F_H1).^2/(pi*F_H2-pi*F_H1)).*((sin(pi*F_H1.*t)./(pi*F_H1.*t)).^2))...
        -(((pi*F_L2).^2/(pi*F_L2-pi*F_L1)).*((sin(pi*F_L2.*t)./(pi*F_L2.*t)).^2)-((pi*F_L1).^2/(pi*F_L2-pi*F_L1)).*((sin(pi*F_L1.*t)./(pi*F_L1.*t)).^2));
    Srct=s1/max(s1);
end

Srcf=fft(Srct);
Q_define=[30,60,80,40,60,70];
c=[2500,2600,2800,2550,2600,2700];
p=[2.3,2.4,2.5,2.35,2.45,2.5];


dh=20;
% h1=[20,40,60,80,100,120,140,160,180,200,220,240,260];
% h2=[20,40,60,80,100,120,140,160,180,200,220,240,260];
% h3=[20,40,60,80,100,120,140,160,180,200,220,240,260];
% h4=[20,40,60,80,100,120,140,160,180,200,220,240,260];
% h5=[20,40,60,80,100,120,140,160,180,200,220,240,260];
% h6=[20,40,60,80,100,120,140,160,180,200,220,240,260];

h1=[20,40,60,80,100,120,140,160,180,200,220,240,260];
h2=[20,40,60,80,100,120,140,160,180,200,220];
h3=[20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320];
h4=[20,40,60,80,100,120,140,160,180,200,220,240,260,280,300];
h5=[20,40,60,80,100,120,140,160,180,200,220,240];
h6=[20,40,60,80,100,120,140,160,180,200,220];

Layer_number=length(Q_define);
for i=1:Layer_number
    eval(['Receiver_spacing(i)=length(h',num2str(i),');']);
end

f=0:fsample/Nt:fsample/2;
DATA=[Srct];
Recf=Srcf;
T=1;
Time=[];
for j=1:Layer_number
    Rect=zeros(Receiver_spacing(j),length(Srct));%介质1中时域接收信号矩阵，每一行代表一�?2000采样�?
    Recf_now=zeros(Receiver_spacing(j),length(Srcf));%介质1中频域接收信号矩阵，每一行代表一�?2000采样�?
    eval(['h=h',num2str(j),';']);
    TravelTime=h/c(j);
    for k=1:1:Receiver_spacing(j)
        Recf_now(k,1:Nt/2+1)=T*Recf(1:Nt/2+1).*exp(-pi*TravelTime(k).*f/Q_define(j)).*exp(-1i*2*pi.*f*TravelTime(k));%接收信号衰减
        Recf_now(k,Nt/2+2:end)=conj(Recf_now(k,Nt/2:-1:2));
        Rect(k,:)=real(ifft(Recf_now(k,:)));
    end
    if j<Layer_number
        T=(2*p(j)*c(j))/(p(j)*c(j)+p(j+1)*c(j+1)); %界面1的�?�射系数
    end
    Recf=Recf_now(Receiver_spacing(j),:);%介质1中频域接收信号矩阵，每一行代表一�?2000采样�?
    DATA=[DATA;Rect];
    Time=[Time,TravelTime./(1:length(h))];
end
Time=Time';

% L = 5;
% dmax  = max(max(abs(DATA)));
% op = hamming(L);
% Noise = conv2(randn(size(DATA)),op,'same');
% Noisemax = max(max(abs(Noise)));
% % Rect_noise = Rect.*(1+Noise*(dmax/Noisemax)/snr);
% DATA_noise = DATA+Noise*(dmax/Noisemax)/snr;

DATA_noise= awgn(DATA,snr,'measured') ;
