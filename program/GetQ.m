function [Q,SpectrumRatio]=GetQ(type,K,BegNum,Bandf,df,ampx,Time)
% type：这个是选择计算Q的方法%
% K：这个是隔几道算一个Q，常用的就是每相邻道算一个，K=1
% BegNum：起始频率的检索
% Bandf：计算Q时所用的频带
% df：频率采样间隔
% ampx：算Q的振幅谱，每一列表示一道的振幅谱
% Time：时间序列，表示各道的初值时间差
TraT=Time;
dm=length(ampx(1,:));
Q=[];
% TravelTime=K*TravelTime;

SpectrumRatio=[];
for i=1:dm-K
    TravelTime=TraT(i);
    SpectrumRatio_i=log(ampx(:,i+K))-log(ampx(:,i));
    SpectrumRatio=[SpectrumRatio,SpectrumRatio_i];
    BandAbsSrcf=ampx(:,i);
    BandAbsRecf=ampx(:,i+K);
    if type==1
        %% 1. Peak frequency shift method
        % [2] Zhang, C. and T.J., Ulrych, 2002, Estimation of quality factors from
        % CMP records:  Gephysics, 67, 1542-1547.        
        [maxValSrc, LocSrcPf]=max(BandAbsSrcf);
        [maxValRec, LocRecPf]=max(BandAbsRecf);
        SrcPf=df*(LocSrcPf+BegNum-1-1);
        RecPf=df*(LocRecPf+BegNum-1-1);
        Qpfs=pi*TravelTime*RecPf*SrcPf^2/2.0/(SrcPf^2-RecPf^2);
        Q_mid=Qpfs;
        
    elseif type==2
        %% 2. Centered frequency shift method
        % [2] Quan, Y. and J.M., Harris, 1997, Seismic attenuation tomography using
        % the frequency shift method: Gephysics, 62, 895-905.        
        SrcCf=sum(Bandf.*BandAbsSrcf)/sum(BandAbsSrcf);
        RecCf=sum(Bandf.*BandAbsRecf)/sum(BandAbsRecf);
        Var2=sum((Bandf-SrcCf).^2.*BandAbsSrcf)/sum(BandAbsSrcf);
        Qcfs=pi*TravelTime*Var2/(SrcCf-RecCf);
        Q_mid=Qcfs;
        
    elseif type==3
        %% 3. Log spectral ratio method
        % [3] White, R.E., 1992, The accuracy of estimating Q from seismic data:
        % Gephysics, 57, 1508-1511.        
        Specratio=log(BandAbsRecf./BandAbsSrcf);
        FitPar=polyfit(Bandf,Specratio,1);
        Qlsr=-pi*TravelTime/FitPar(1);     
%         figure(2)
%         plot(Bandf,Specratio,'o')
%         hold on
%         fitcur=FitPar(1)*Bandf+FitPar(2);
%         plot(Bandf,fitcur,'r')
%         title(['Qlsr=',num2str(Qlsr)])
        Q_mid=Qlsr;
        
    elseif type==4
        %% 4. Improved frequency shift method (IFS)
        % [4] Hu. C., N. Tu, and W. Lu, 2013, Seismic attenuation estimation using
        % an improved frequency shift method: IEEE Geoscience and remote sensing
        % letters, 10, 1026-1030.      
        SrcCf=sum(Bandf.*BandAbsSrcf)/sum(BandAbsSrcf);
        RecCf=sum(Bandf.*BandAbsRecf)/sum(BandAbsRecf);
        SrcIpf=pi^0.5*SrcCf/2.0;
        RecIpf=pi^0.5*RecCf/2.0;
        Qifs=pi*TravelTime*RecIpf*SrcIpf^2/2.0/(SrcIpf^2-RecIpf^2);
        Q_mid=Qifs;
        
    elseif type==5
        %%  5. Domaint-central frequency shift method (DCFS)
        % [5]Li, F., H. Zhou, N. J, J. Bi and K.J. Marfurt,2015, Q estiamtion from
        % reflection seismic data for hydrocarbon detection using a modified
        % frequency shift method: Journal of Geophysics and engeering, 557-586.
        %        
        [maxValSrc, LocSrcDf]=max(BandAbsSrcf);
        SrcDf=df*(LocSrcDf+BegNum-1-1);
        PowerRecCf=sum(Bandf.*BandAbsRecf.^2)/sum(BandAbsRecf.^2);
        Qdcfs=-SrcDf^2*pi*TravelTime/(4*SrcDf*(2/pi)^0.5-3*PowerRecCf*pi^0.5);
        Q_mid=Qdcfs;   
        
    elseif type==6
        %% 6. frequency shift method by FWE
        %[6] Li, C. and X. Liu, 2015,A new method for interval Q-factor inversion from
        %seismic refelction data: Geophyscis, 361-373.      
        SrcCf=sum(Bandf.*BandAbsSrcf)/sum(BandAbsSrcf);
        RecCf=sum(Bandf.*BandAbsRecf)/sum(BandAbsRecf);
        VarSrc=sum((Bandf-SrcCf).^2.*BandAbsSrcf)/sum(BandAbsSrcf);
        VarRec=sum((Bandf-RecCf).^2.*BandAbsRecf)/sum(BandAbsRecf);
        nsrc=SrcCf^2/VarSrc-1;
        nrec=RecCf^2/VarRec-1;
        avern=(nsrc+nrec)/2;
        corrSrcFb=SrcCf/(avern+1);
        corrRecFb=RecCf/(avern+1);
        Qfwe=pi*TravelTime*corrSrcFb*corrRecFb/(corrSrcFb-corrRecFb);
        Q_mid=Qfwe;  
        
    elseif type==7
        %% 7. logarithmic spectral area difference
        %[7] Wang, S., D. Yang, J. Li, and H. Song, 2015,Q factor estimation based on
        %the method of logarithmic spectral are difference: Geophysics, 157-171.        
        SrcArea=sum(log(BandAbsSrcf));
        RecArea=sum(log(BandAbsRecf));
        AreaDiff=SrcArea-RecArea;
        BandArea=sum(Bandf);
        Qlsad=pi*TravelTime*BandArea/AreaDiff;
        Q_mid=Qlsad;  
        
    elseif type==8
        %% 8. Weighted centroid frequency shift method
        %[8] Li, J., D. Yang, S. Wang, S. Wang, Z. Chen, and P. Li, 2015, An improved Q
        %estiamtion approach: weighted centroid frequency shift: SEG New Orleans
        %Annual Meeting, 5605-5609.
        SrcCf=sum(Bandf.*BandAbsSrcf)/sum(BandAbsSrcf);
        SrcVar=sum((Bandf-SrcCf).^2.*BandAbsSrcf)/sum(BandAbsSrcf);
        Weight=exp(-(Bandf-SrcCf).^2/2/SrcVar);
        WeightSrcCf=sum(Weight.*Bandf.*BandAbsSrcf)/sum(Weight.*BandAbsSrcf);
        WeightRecCf=sum(Weight.*Bandf.*BandAbsRecf)/sum(Weight.*BandAbsRecf);
        WeightSrcVar=sum(Weight.*(Bandf-SrcCf).^2.*BandAbsSrcf)/sum(Weight.*BandAbsSrcf);
        Qwcfs=pi*TravelTime*WeightSrcVar/(WeightSrcCf-WeightRecCf);       
        Q_mid=Qwcfs;  
    end
    Q=[Q,Q_mid];
end
% Q=1./Q;  
Q=1./Q/2;    