clear all;
close all;
clc;
t=0:1/1023:1;
x=cos(2*pi*50*t)+cos(2*pi*100*t)+(eps*i);
x=x';
fs=1/1024;
mex test_mex.c;

tic;[a,b]=test_mex(x);toc
tic;y=fft(x,1024);toc
tic;z=ifft(y,1024);toc
f=linspace(0,1,1024)*2*pi;

t_1=[1/fs:1/fs:length(x)/fs];
subplot(3,2,[1,2]);plot(t_1,x);title('signal');
subplot(3,2,3);plot(f,abs(y));title('matlab fft');
subplot(3,2,4);plot(f,abs(a));title('mex fft');
subplot(3,2,5);plot(t_1,z);title('matlab ifft');
subplot(3,2,6);plot(t_1,b);title('mex ifft');
