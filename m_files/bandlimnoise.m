function [stim,nt] = bandlimnoise (cutoff,tsamp,tmax)
%
%%%%%%%%%%%%%%% Generate band-limited noise with a cutoff frequency %%%%%%%%%%%%%%%
nt=round(tmax/tsamp)+1;
z=round(log2(nt)-0.01);
nt=2^z;
x=randn(1,nt);
xf=fft(x,2^z); 
clear x;
df=1.0/tsamp/2^z;
j1=1;
j2=round(cutoff/df)+1;
nn=length(xf);
xf(1:j1-1)=0; xf(j2+1:nn)=0;
stim=ifft(xf,'symmetric'); clear xf;
stim=stim-mean(stim);
stim=stim/std(stim);
return