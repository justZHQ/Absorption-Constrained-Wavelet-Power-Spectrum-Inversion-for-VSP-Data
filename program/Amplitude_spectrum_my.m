function [f,amplitude_spectrum]=Amplitude_spectrum_my(dt,signal,fmax)
% dt;%时间采样/s
% signal;%输入时域信号
% fmax;%最大显示频率/Hz
% f;%输出频率序列/Hz
% amplitude_spectrum;%输入信号signal的振幅谱
F=fft(signal,2^ceil(log2(1000)));
N=length(F);
amplitude_spectrum=abs(F);%abs(Cn)=abs(X(k))=sqrt(realX(k)^2+imagX(k)^2)
f=(0:N-1)*(1/(N*dt));%频率采样间隔等于基波频率f0=1/T;T=dt*N=N/fs
fmax_number=ceil(fmax/(1/(N*dt)));
f=f(1:fmax_number);
amplitude_spectrum=amplitude_spectrum(1:fmax_number)*2/N;%频率与真实振幅之间的关系——振幅谱；由于abs(Cn)=abs(Xn)=1/2*An,且在离散时间序列傅里叶变换：X(n)到X(w)时是为了方便在X(w)前加了1/N;因此为了恢复各频率对应的真振幅，
                              %应将mag=abs(y)乘以2/N; 对于实数序列，离散时间序列的傅里叶变换的X(k)关于N/2具有对称性;因此只要写出Nyquist频率之前的N/2个频率的振幅
