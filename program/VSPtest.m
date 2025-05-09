clc
clear all;
close all;
%%
dt=0.001;       % sample time
Nt=1000;        % number of samples
fsample=1/dt;   % number of samples
df=fsample/Nt;
f=0:fsample/Nt:fsample/2;
K=1;
wavelet_type='Ricker';
% wavelet_type='Ormsby';

if wavelet_type=='Ricker'
    Begf=15;Endf=65;
    sigma=1;
    average=1;
    GausiannVariance=50;
    snr_my=17;
    lambda=0.001; % Spectrum Constraint (SC) Weight Coefficient. The smaller the lambda, the greater the weight of the Absorption Constraint (AC) term
elseif wavelet_type=='Ormsby'
    Begf=25;Endf=65;
    sigma=1;
    average=1;
    GausiannVariance=35;
    snr_my=15;
    lambda=0.001; % Spectrum Constraint (SC) Weight Coefficient. The smaller the lambda, the greater the weight of the Absorption Constraint (AC) term
end
%%
[DATA,DATA_noise,dh,Time]=SyntheticData(wavelet_type,snr_my,dt,Nt);
NUMBER=length(DATA(:,1))-1;
h_all=0:dh:NUMBER*dh;h_all=h_all*0.001;
T=0:1:Nt-1;
%%
[BegNum,Bandf,ff,log_ampx_orignal,ampx_orignal]=Effective_Amplitude_Spectrum(df,f,DATA,Begf,Endf);
[~,~,~,log_ampx,ampx]=Effective_Amplitude_Spectrum(df,f,DATA_noise,Begf,Endf);
[~,~,~,log_ampx_ture,~]=Effective_Amplitude_Spectrum(df,f,DATA,Begf,Endf);

[log_ampx_Correction,ampx_Correction,Gaussian]=new_Fast(sigma,lambda,average,GausiannVariance,log_ampx);
% figure;
% plot(Gaussian)

ampx_orignal1=zeros(size(ampx_orignal));
ampx1=zeros(size(ampx));
ampx_Correction1=zeros(size(ampx_Correction));
for i=1:length(DATA(:,1))
ampx_orignal1(:,i)=ampx_orignal(:,i)./max(ampx_orignal(:,i));
ampx1(:,i)=ampx(:,i)./max(ampx(:,i));
ampx_Correction1(:,i)=ampx_Correction(:,i)./max(ampx_Correction(:,i));
end
[m,n]=size(log_ampx);
% n: Number of geophones
% m: Number of frequency samples
%%
fit_order=5;
ampx_fit=zeros(size(ampx));
p=zeros(n,fit_order+1);
for i=1:n
p(i,:) = polyfit(ff,ampx(:,i),fit_order);
ampx_fit(:,i) = polyval(p(i,:),ff);
% ampx_fit(:,i)=ampx_fit(:,i)./max(ampx_fit(:,i));
end
log_ampx_fit=log(abs(ampx_fit));

ampx_fit1=zeros(size(ampx_fit));
for i=1:length(DATA(:,1))
ampx_fit1(:,i)=ampx_fit(:,i)./max(ampx_fit(:,i));
end

%% 2. Centered frequency shift method
type=2;
[Qcfs_Orignal,SpectrumRatio_Orignal]=GetQ(type,K,BegNum,Bandf,df,ampx_orignal,Time);
[Qcfs,SpectrumRatio]=GetQ(type,K,BegNum,Bandf,df,ampx,Time);
[Qcfs_Correction,SpectrumRatio_Correction]=GetQ(type,K,BegNum,Bandf,df,ampx_Correction,Time);
[Qcfs_fit,SpectrumRatio_fit]=GetQ(type,K,BegNum,Bandf,df,ampx_fit,Time);
%% 3. Log spectral ratio method
type=3;
[Qlsr_Orignal,~]=GetQ(type,K,BegNum,Bandf,df,ampx_orignal,Time);
[Qlsr,~]=GetQ(type,K,BegNum,Bandf,df,ampx,Time);
[Qlsr_Correction,~]=GetQ(type,K,BegNum,Bandf,df,ampx_Correction,Time);
[Qlsr_fit,~]=GetQ(type,K,BegNum,Bandf,df,ampx_fit,Time);
%%
gcf1=figure;
% wigb_my(-DATA_noise',1,h,T*0.001,1);
imagesc(T*0.001,h_all,DATA_noise);
xlabel('Time/s');
ylabel('Depth/km');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(gray);
cb=colorbar;
set(gcf,'unit','normalized','position',[0.1,0.1,0.25,0.5]);
set(gca,'position',[0.17 0.14 0.65 0.84]);
caxis([-0.1 0.1]);
% set(cb,'YTickLabel',{'-0.1','0','0.1'});
cb.Ticks=linspace(-0.1,0.1,3)';
% xticks([0,0.5,1.0]);
axis([0 1,0 h_all(end)]);

DATA_noise_extract=DATA_noise([1,40,79],:);
snr1=snr(DATA(1,:),DATA_noise(1,:)-DATA(1,:));
snr2=snr(DATA(40,:),DATA_noise(40,:)-DATA(40,:));
snr3=snr(DATA(79,:),DATA_noise(79,:)-DATA(79,:));
DATA_noise_extract=DATA_noise_extract./max(DATA_noise_extract,[],2);

gcf2=figure;
wigb_my(-DATA_noise_extract',1,h_all([1,40,79]),T*0.001,1);
xlabel('Time/s');
ylabel('Depth/km');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gcf,'unit','normalized','position',[0.1,0.1,0.2,0.5]);
set(gca,'position',[0.215 0.14 0.77 0.85]);
axis([0 1,-0.7 2.2]);
yticks([0,0.5,1,1.5]);
text(0.67,0.2,['SNR=',num2str(floor(snr1))],'FontName','Times New Roman','FontSize',20,'linewidth',5,'FontWeight','bold');
text(0.67,1,['SNR=',num2str(floor(snr2))],'FontName','Times New Roman','FontSize',20,'linewidth',5,'FontWeight','bold');
text(0.67,1.8,['SNR=',num2str(floor(snr3))],'FontName','Times New Roman','FontSize',20,'linewidth',5,'FontWeight','bold');

gcf3=figure;box on;
imagesc(h_all,ff,ampx_orignal1);hold on;
xlabel('Depth/km');
ylabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
c=colorbar;
c.Ticks = 0:1:1;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.4]);
set(gca,'position',[0.155 0.16 0.72 0.82]);
caxis([0 1]);
stem(h_all(20),ff(end),'--k','linewidth',4,'Marker','none');
% stem(h(40),ff(end),'k','linewidth',2,'Marker','none');
stem(h_all(60),ff(end),'--k','linewidth',4,'Marker','none');

gcf4=figure;box on;
imagesc(h_all,ff,ampx1);hold on;
xlabel('Depth/km');
ylabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
c=colorbar;
c.Ticks = 0:1:1;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.4]);
set(gca,'position',[0.155 0.16 0.72 0.82]);
caxis([0 1]);
stem(h_all(20),ff(end),'--k','linewidth',4,'Marker','none');
% stem(h(40),ff(end),'k','linewidth',2,'Marker','none');
stem(h_all(60),ff(end),'--k','linewidth',4,'Marker','none');

gcf5=figure;box on;
imagesc(h_all,ff,ampx_fit1);hold on;
xlabel('Depth/km');
ylabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
c=colorbar;
c.Ticks = 0:1:1;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.4]);
set(gca,'position',[0.155 0.16 0.72 0.82]);
caxis([0 1]);
stem(h_all(20),ff(end),'--k','linewidth',4,'Marker','none');
% stem(h(40),ff(end),'k','linewidth',2,'Marker','none');
stem(h_all(60),ff(end),'--k','linewidth',4,'Marker','none');

gcf6=figure;box on;
imagesc(h_all,ff,ampx_Correction1);hold on
xlabel('Depth/km');
ylabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
c=colorbar;
c.Ticks = 0:1:1;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.4]);
set(gca,'position',[0.155 0.16 0.72 0.82]);
caxis([0 1]);
stem(h_all(20),ff(end),'--k','linewidth',4,'Marker','none');
% stem(h(40),ff(end),'k','linewidth',2,'Marker','none');
stem(h_all(60),ff(end),'--k','linewidth',4,'Marker','none');
%%
cr1=diag(corr(ampx1,ampx_orignal1));
cr2=diag(corr(ampx_fit1,ampx_orignal1));
cr3=diag(corr(ampx_Correction1,ampx_orignal1));
crone=ones(size(cr1));
gcf7=figure;
hold on; box on;
plot(h_all,cr1,'b','linewidth',3);
plot(h_all,cr2,'g','linewidth',3);
plot(h_all,cr3,'-.k','linewidth',3);
plot(h_all,crone,'r','linewidth',3);
legend('Measured value','Polynomial fitting','AWPSI','Location','SouthWest');
ylabel('Correlation Coefficient');
xlabel('Depth/km');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);    
set(gcf,'unit','normalized','position',[0.1,0.1,0.34,0.35]);
set(gca,'position',[0.21 0.16 0.75 0.82]);
% axis([0 h(end),0.95 1.01]);
axis([0 h_all(end),0.975 1.001]);yticks([0.98,0.99,1]);
%%
i1=20;
i2=60;
gcf8=figure;hold on;box on;
plot(ff,ampx_orignal(:,i1),'r','linewidth',2.5);
plot(ff,ampx(:,i1),'b','linewidth',2.5);
plot(ff,ampx_fit(:,i1),'g','linewidth',2.5);
plot(ff,ampx_Correction(:,i1),'-.k','linewidth',2.5);
% legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','SouthWest');
legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','south');
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);    
set(gcf,'unit','normalized','position',[0.3,0.3,0.22,0.4]);
set(gca,'position',[0.19 0.13 0.78 0.85]);
xlim([ff(1) ff(end)]);

gcf9=figure;hold on;box on;
plot(ff,ampx_orignal(:,i1+1),'r','linewidth',2.5);
plot(ff,ampx(:,i1+1),'b','linewidth',2.5);
plot(ff,ampx_fit(:,i1+1),'g','linewidth',2.5);
plot(ff,ampx_Correction(:,i1+1),'-.k','linewidth',2.5);
% legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','SouthWest');
legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','south');
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);    
set(gcf,'unit','normalized','position',[0.3,0.3,0.22,0.4]);
set(gca,'position',[0.19 0.13 0.78 0.85]);
xlim([ff(1) ff(end)]);

gcf10=figure;hold on;box on;
plot(ff,ampx_orignal(:,i2),'r','linewidth',2.5);
plot(ff,ampx(:,i2),'b','linewidth',2.5);
plot(ff,ampx_fit(:,i2),'g','linewidth',2.5);
plot(ff,ampx_Correction(:,i2),'-.k','linewidth',2.5);
% legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','northeast');
legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','south');
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);    
set(gcf,'unit','normalized','position',[0.3,0.3,0.22,0.4]);
set(gca,'position',[0.19 0.13 0.78 0.85]);
xlim([ff(1) ff(end)]);

gcf11=figure;hold on;box on;
plot(ff,ampx_orignal(:,i2+1),'r','linewidth',2.5);
plot(ff,ampx(:,i2+1),'b','linewidth',2.5);
plot(ff,ampx_fit(:,i2+1),'g','linewidth',2.5);
plot(ff,ampx_Correction(:,i2+1),'-.k','linewidth',2.5);
% legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','northeast');
legend('Theoretical value','Measured value','Polynomial fitting','AWPSI','Location','south');
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);    
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.4]);
set(gca,'position',[0.19 0.13 0.78 0.85]);
xlim([ff(1) ff(end)]);

gcf12=figure;hold on;box on;
FitPar=polyfit(ff,log_ampx(:,21)-log_ampx(:,20),1);
plot(ff,log_ampx_orignal(:,21)-log_ampx_orignal(:,20),'r','linewidth',2.5);
plot(ff,log_ampx(:,21)-log_ampx(:,20),'b','linewidth',2.5);
plot(ff,FitPar(2)+FitPar(1)*ff,'m','linewidth',2.5);
plot(ff,log_ampx_fit(:,21)-log_ampx_fit(:,20),'g','linewidth',2.5);
plot(ff,log_ampx_Correction(:,21)-log_ampx_Correction(:,20),'-.k','linewidth',2.5);
legend('Theoretical value','Measured value','Traditional','Polynomial fitting','AWPSI','Location','SouthWest');
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);    
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.4]);
set(gca,'position',[0.19 0.13 0.78 0.85]);
xlim([ff(1) ff(end)])
yticks([-0.3,-0.2,-0.1,0,0.1,0.2,0.3]);

gcf13=figure;hold on;box on;
FitPar=polyfit(ff,log_ampx(:,61)-log_ampx(:,60),1);
plot(ff,log_ampx_orignal(:,61)-log_ampx_orignal(:,60),'r','linewidth',2.5);
plot(ff,log_ampx(:,61)-log_ampx(:,60),'b','linewidth',2.5);
plot(ff,FitPar(2)+FitPar(1)*ff,'m','linewidth',2.5);
plot(ff,log_ampx_fit(:,61)-log_ampx_fit(:,60),'g','linewidth',2.5);
plot(ff,log_ampx_Correction(:,61)-log_ampx_Correction(:,60),'-.k','linewidth',2.5);
legend('Theoretical value','Measured value','Traditional','Polynomial fitting','AWPSI','Location','SouthWest');
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);    
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.4]);
set(gca,'position',[0.19 0.13 0.78 0.85]);
xlim([ff(1) ff(end)]);

%%
h=h_all(K+1:end)';
% h=h(1:K:end)
gcf14=figure;
hold on;box on;
plot(h,Qcfs_Orignal,'r','linewidth',2);
plot(h,Qcfs,'b','linewidth',2);
plot(h,Qcfs_fit,'-.g','linewidth',2);
plot(h,Qcfs_Correction,'k','linewidth',2);
xlabel('Depth/km');
ylabel('1/Q value');
axis([0 h(end),-0.015 0.05]);
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
axis([0 h(end),-0.015 0.05]);
lgd=legend({'Real value','Traditional CFS','Polynomial fitting','AWPSI'},...
    'FontSize',15,'TextColor','black','Location','north');
% lgd.NumColumns = 2

gcf15=figure;
hold on;box on;
plot(h,Qlsr_Orignal,'r','linewidth',2);
plot(h,Qlsr,'b','linewidth',2);
plot(h,Qlsr_fit,'-.g','linewidth',2);
plot(h,Qlsr_Correction,'k','linewidth',2);
xlabel('Depth/km');
ylabel('1/Q value');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
axis([0 h(end),-0.015 0.05]);
lgd=legend({'Real value','Traditional LSR','Polynomial fitting','AWPSI'},...
    'FontSize',15,'TextColor','black','Location','north');
% lgd.NumColumns = 2;
%%
[H_orignal]=H_matrix(log_ampx_orignal);
[H]=H_matrix(log_ampx);
[H_fit]=H_matrix(log_ampx_fit);
[H_Correction]=H_matrix(log_ampx_Correction);
gcf16=figure;box on;
imagesc(ff(1:end-1),h(1:end),H_orignal);
ylabel('Depth/km');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
colorbar;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.35]);
set(gca,'position',[0.18 0.175 0.6 0.8]);
caxis([-0.003 0.001]);
hold on;
stem(ff(20),h(end),'--k','linewidth',4,'Marker','none');
plot(ff(1)-1:ff(end),ones(size(ff(1)-1:ff(end)))*h(47),'--k','linewidth',4);

gcf17=figure;box on;
imagesc(ff(1:end-1),h(1:end),H);
ylabel('Depth/km');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
colorbar;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.35]);
set(gca,'position',[0.18 0.175 0.6 0.8]);
% caxis([-0.002 0.0001]);
caxis([-1 1]);
hold on;
stem(ff(20),h(end),'--k','linewidth',4,'Marker','none');
plot(ff(1)-1:ff(end),ones(size(ff(1)-1:ff(end)))*h(47),'--k','linewidth',4);

gcf18=figure;box on;
imagesc(ff(1:end-1),h(1:end),H_fit);
ylabel('Depth/km');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
colorbar;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.35]);
set(gca,'position',[0.18 0.175 0.6 0.8]);
% caxis([-0.1 0.1]);
% caxis([-1 1]);
caxis([-0.003 0.001]);
hold on;
stem(ff(20),h(end),'--k','linewidth',4,'Marker','none');
plot(ff(1)-1:ff(end),ones(size(ff(1)-1:ff(end)))*h(47),'--k','linewidth',4);

gcf19=figure;box on;
imagesc(ff(1:end-1),h(1:end),H_Correction);
ylabel('Depth/km');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
colormap(jet);
colorbar;
set(gcf,'unit','normalized','position',[0.1,0.1,0.22,0.35]);
set(gca,'position',[0.18 0.175 0.6 0.8]);
caxis([-0.003 0.001]);
% caxis([-1 1]);
hold on;
stem(ff(20),h(end),'--k','linewidth',4,'Marker','none');
plot(ff(1)-1:ff(end),ones(size(ff(1)-1:ff(end)))*h(47),'--k','linewidth',4);

k=47;
gcf20=figure;
positionVector1 = [0.12, 0.67, 0.87, 0.3];
subplot('Position',positionVector1);
hold on;box on;
plot(ff(1:end-1),H_orignal(k,:),'r','linewidth',2);
plot(ff(1:end-1),H(k,:),'b','linewidth',2);
plot(ff(1:end-1),H_fit(k,:),'-.g','linewidth',2);
plot(ff(1:end-1),H_Correction(k,:),'--k','linewidth',2);
lgd=legend({'Real value','Measured value','Polynomial fitting','AWPSI'},...
    'FontSize',18,'TextColor','black','Location','northoutside');
legend('boxoff');lgd.NumColumns = 4;
ylabel('Magnitude');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'xticklabel',[]);
positionVector1 = [0.12, 0.18, 0.87, 0.4];
subplot('Position',positionVector1);
hold on;box on;
plot(ff(1:end-1),H_orignal(k,:),'r','linewidth',2);
plot(ff(1:end-1),H(k,:),'b','linewidth',2);
plot(ff(1:end-1),H_fit(k,:),'-.g','linewidth',2);
plot(ff(1:end-1),H_Correction(k,:),'--k','linewidth',2);
ylabel('Magnitude');
xlabel('Frequency/Hz');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
axis([0 h(end),-0.02 0.02]);
set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.35]);
xlim([ff(1) ff(end-1)]);

g=20;
gcf21=figure;
positionVector1 = [0.12, 0.67, 0.87, 0.33];
subplot('Position',positionVector1);
hold on;box on;
plot(h(1:end),H_orignal(:,g),'r','linewidth',2);
plot(h(1:end),H(:,g),'b','linewidth',2);
plot(h(1:end),H_fit(:,g),'-.g','linewidth',2);
plot(h(1:end),H_Correction(:,g),'--k','linewidth',2);
lgd=legend({'Real value','Measured value','Polynomial fitting','AWPSI'},...
    'FontSize',18,'TextColor','black','Location','northoutside');
legend('boxoff');lgd.NumColumns = 4;
ylabel('Magnitude');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
set(gca,'xticklabel',[])
positionVector1 = [0.12, 0.18, 0.87, 0.4];
subplot('Position',positionVector1);
hold on;box on;
plot(h(1:end),H_orignal(:,g),'r','linewidth',2);
plot(h(1:end),H(:,g),'b','linewidth',2);
plot(h(1:end),H_fit(:,g),'-.g','linewidth',2);
plot(h(1:end),H_Correction(:,g),'--k','linewidth',2);
ylabel('Magnitude');
xlabel('Depth/km');
set(gca,'FontName','Times New Roman','FontSize',20,'linewidth',3,'FontWeight','bold');
set(gca,'TickLength',[0 0.001]);
axis([0 h(end),-0.0025 0.0015]);
set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.35]);

