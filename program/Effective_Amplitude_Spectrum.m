function [BegNum,Bandf,ff,log_ampx,ampx]=Effective_Amplitude_Spectrum(df,f,DATA,Begf,Endf)
[dm,~]=size(DATA);
% frequency band parameter slection
BegNum=floor(Begf/df)+1; 
EndNum=floor(Endf/df)+1;
Bandf=((BegNum-1):(EndNum-1))*df;
Bandf=Bandf';
log_ampx=[];
ampx=[];
for i=1:dm
AbsSrcf=abs(fft(DATA(i,:)));
AbsSrcf=AbsSrcf.^2;%功率�?
BandAbsSrcf=AbsSrcf(BegNum:EndNum);
log_ampx=[log_ampx,log(BandAbsSrcf')];
ampx=[ampx,BandAbsSrcf'];
end
ff=f(BegNum:EndNum)';
