% creates flat spectrum chirps and gates

f1 = 100;	% f start
f2 = 10400;	% f stop
fsamp = 44100;	% sampling frequency
dur = 0.012;	% duration of chirp (sec)
alpha = 3.0;	% exponent alpha
lwin = 0.001;	% duration of ramps (sec)

% call to bmchirp with var. alpha, find last non-zero sample

c_up = alchirp(f1,f2,fsamp,dur,alpha);
n = length(c_up);
for i = n :-1:1
	if c_up(i) ~= 0
		sample = i;
		break;
        end;
end;

% produce flat chirp via fft and ifft, only real part is valid
clear i % added (Till Habersetzer 22.10.21) otherwise fft doesn't work

fftc_up = fft(c_up);
n = length(c_up);
fc_up = real(ifft(cos(angle(fftc_up)) + i * sin(angle(fftc_up))));

% cut off zeros for gating, gate

c_up = fc_up(1:sample);
n = length(c_up);
nwin = 2*round(fsamp*lwin);
window = hanning(nwin)';
plateau = 1 + zeros(1,n-nwin);
gate = [ window(1:nwin/2) plateau window(nwin/2+1:nwin) ];
c_up = c_up.* gate;
c_dn = fliplr(c_up);

% write wav files

audiowrite('up.wav',c_up/max(abs(c_up)),fsamp);
audiowrite('down.wav',c_dn/max(abs(c_up)),fsamp);
