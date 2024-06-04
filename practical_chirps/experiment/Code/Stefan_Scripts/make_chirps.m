% function make_chirps (fstart, fstop, fsamp, period, a, duration, scalefactor)
function make_chirps (fstart, fstop, fsamp, period, a, duration, scalefactor)

% make_chirps (fstart, fstop, fsamp, period, a, duration, scalefactor)
%
% make chirps, using function sig = alchirp(f1,f2,fs,len,alpha),
% which is a modified version (varibale alpha) of the optimised chirp 
% routine supplied by Torsten Dau and Oliver Wegner.
%
% S. Uppenkamp, 12/01/99
%
%         f1 : lower start frequency of the chirp (in Hz)
%         f2 : upper stop frequncy of the chirp (in Hz)
%         fs : sampling frequency in Hz (standard: 25000 Hz)
%         len: length of the chirp in s (standard: 0 s)
%              determines period in seconds, pitch is 1/len
%         alpha: in-/de-crease of instantaneous frequency
%         alpha=3.0  ----> bmchirp (to compensate dispersion as
%                          described above) 
%
% chirp parameters and examples
%
%  fstart = 100;   % Hz
%  fstop = 10000;  % Hz

% fsamp = 80000;  % Hz

%  period = 0.02;  % sec
%  a = 2.0;        % cm^^-1

% additional stimulus paramters
%
% duration = 0.5;      % overall duration of sound in sec
% scalefactor = 30000;  % appropriate size for 16 bit D/A

% start, make up one upward chirp 
%
close all;
clear cb stim rstim;
cb=alchirp(fstart,fstop,fsamp,period,a);
figure(1);
plot(cb);

% scale to appropriate size
%           
cb=cb*scalefactor;

% time reverse, keep first non-zero for later
%
n=size(cb);
for i=1:n(2)
  r_cb(i)=cb(n(2)+1-i);
  if r_cb(i) ~= 0
    if r_cb(i-1) == 0
      first = i;
    end;
  end;
end;

% swap zero and non-zero part using "first", we don't want to start with zeros
% first < 3 doesn't work out here, so skip the swap then
%
if first > 2
  r_begin = r_cb(first-1:n(2));
  r_rest = r_cb(1:(first-2));
  r_cb = [r_begin r_rest];
end;

% how many in overall duration?
%
m = duration/period;

% fill stimulus arrays 
%
for i = 1:m
  j = i*n(2);
  stim(j-n(2)+1:j) = cb;
  rstim(j-n(2)+1:j) = r_cb;
end;

% pictures ...
%
figure(2);
plot(stim);
figure(3);
plot(rstim);

% ... and finally two files with PC integer data in it
%
str = sprintf('%1.1f',a);
strr = [str(1) str(3)];
f_up = ['up_' strr '.pcr'];
f_dn = ['dn_' strr '.pcr'];

d = fopen(f_up,'w','l'); 
fwrite(d,stim,'integer*2');
st=fclose(d);
d = fopen(f_dn,'w','l'); 
fwrite(d,rstim,'integer*2');
st = fclose(d);

% full stop.
