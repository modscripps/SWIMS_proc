%function newsig=fftzero(sig,lowbin,highbin)
%Zeros all of the bins in the fft of the signal from lowbin to highbin.
%Remember if w is the sample frequency, then the nth bin contains the
%(n/N)*w frequency component.
%
%Pass highbin=-1 for standard lowpass filter.
%Note: doesn't work for complex signal as I take real part at end
%if sig is a matrix, it does it to each column

function newsig=fftzero(sig,lowbin,highbin)


[N,n]=size(sig);
N
n
if (highbin==-1) 
	highbin=N/2-1;
elseif (highbin> N)
	highbin =N/2-1;	
end

sfft=fft(sig);
sfftsh=fftshift(sfft);

zing=zeros(highbin-lowbin+1,n);
sfftsh(N/2+1+lowbin:N/2+1+highbin,:)=zing;
sfftsh(N/2+1-highbin:N/2+1-lowbin,:)=zing;

newsig=real(ifft(fftshift(sfftsh)));
