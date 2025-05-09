function [t,ricker,f,amplitude_spectrum]=Ricker_my(dt,L,rick_mainf,fa,fmax)
% dt;%时间采样/s
% L;%子波长度；%wavelength=2*L+1
% rick_mainf;%雷克子波的主频/HZ
% fa;%相位
% fmax;%最大显示频率/Hz
% t;%输出时间序列/s
% ricker;%输出时间域雷克子波
% f;%输出频率序列
% amplitude_spectrum;%输出雷克子波的振幅谱

t=-L*dt:dt:L*dt;
s1=(1-2*(pi*rick_mainf*t).^2).*exp(-(pi*rick_mainf*t).^2);
s2=hilbert(s1);
ricker=real(s2)*cos(fa)+imag(s2)*sin(fa);
[f,amplitude_spectrum]=Amplitude_spectrum_my(dt,ricker,fmax);