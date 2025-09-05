function y = bandpass_iir(x, Fs, f1, f2, order)
Wn = sort([f1 f2])/(Fs/2);
if Wn(1)<=0, Wn(1)=1e-6; end
if Wn(2)>=1, Wn(2)=0.999; end
[b,a] = butter(order, Wn, 'bandpass');   % nécessite Signal Processing Toolbox
y = filtfilt(b,a,x);
end
