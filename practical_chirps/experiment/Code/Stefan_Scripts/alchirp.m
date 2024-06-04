function sig = alchirp(f1,f2,fs,len,alpha)
%BMCHIRP  Chirp which compensates the spatial dispersion on basilar membrane
%
%         sig = alchirp(f1,f2,fs,len,alpha)
%         generates a chirp, which should compensate the spatial
%         dispersion on the basilar membrane according to the
%         linear basilar membrane model of de Boer (1980). A
%         vector sig with the samples is returned. The vector is
%         filled with zeros to match the length of len. If len
%         equals zero the minimum length for the chirp will be
%         calculated.
%        
%         f1 : lower start frequency of the chirp (in Hz)
%         f2 : upper stop frequncy of the chirp (in Hz)
%         fs : sampling frequency in Hz (standard: 25000 Hz)
%         len: length of the chirp in s (standard: 0 s)
%         
%         alpha: in-/de-crease of instantaneous frequency
%         alpha=3.0  ----> bmchirp (to compensate dispersion as
%                          described above) 
%
%         Based upon: E de Boer (1980): "Auditory physics.
%                     Physical principles in hearing theory I",
%                     Phys Rep (62), 87-274
%
%         Frequenz-Orts-Transformation:
%         Based upon: DD. Greenwood(1961): "Critical bandwidth
%                     and the frequency coordinates of the
%                     basilar membrane",
%                     J Acoust Soc Am (33), 1344-1356
%
% author/date: ow,td / 09.96
% last update: ow,td / 10.97
%
% variable alpha:   S. Uppenkamp, 18-09-98, Standard: alpha = 3.0

% All numbered equations refer to the corresponding equations in
% Chapter 5 of the diploma thesis of Oliver Wegner (1997): "Zu-
% sammenhang zwischen Psychoakustik und Hirnstammpotentialen"

% check for f1 and f2
if (exist('f1') ~= 1) | (exist('f2') ~= 1)
  help alchirp
  return
end

% range check
if (f1 < 0) | (f2 > 20457)	% (f > 20457) => (bmpos61(f) < 0)
  error('f1 and f2 must be in the range [0 ; 20457]');
end
if (f1 > f2)
  error('f1 > f2: not allowed');
end

% standard values:
len_std = 0.0;			 % standard value for len (means minimum length)
fs_std = 25000.0;		 % standard value for fs
alpha_std = 3.0;		 % standard value for alpha

% check for len
if (exist('len') ~= 1)
  len = len_std;
elseif isempty(len)
  len = len_std;
end

% check for fs
if (exist('fs') ~= 1)
  fs = fs_std;
elseif isempty(fs)
  fs = fs_std;
end

% check for alpha
if (exist('alpha') ~= 1)
  alpha = alpha_std;
elseif isempty(alpha)
  alpha = alpha_std;
end

len = ceil(len*fs);		 % len in samples (rounded towards next integer)
t0 = bmtime(0,alpha);			 % time for 0 Hz (reference)
t1 = t0 - bmtime(f1,alpha);		 % time for f1
t2 = t0 - bmtime(f2,alpha);		 % time for f2
samples = ceil((t2 - t1) * fs);	 % samples needed (rounded towards next integer)
if (len == 0) len = samples; end % len = 0 => optimal fill
if (samples > len)		 % is len long enough?
  error(sprintf('at least %f s are needed',samples/fs));
end
t = t1 : 1/fs : t2;		 % time-axis

sig = [sin(instphas(t,alpha)-instphas(t1,alpha)) zeros(1,len-length(t))];	
% sig = [(instphas(t,alpha) - instphas(t - 1/fs,alpha))*fs/(2*pi) zeros(1,len-length(t))];	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function time = bmtime(f,alfa)
%BMTIME   traveling time on the basilar membrane
%
%         time = bmtime(f,alfa)
%         returns the time (in s) a pulse needs to get to the
%         position on the basilar membrane which represents the
%         frequency f (in Hz).
%
%         Based upon: E de Boer (1980): "Auditory physics.
%                     Physical principles in hearing theory I",
%                     Phys Rep (62), 87-274
%
% author/date: ow,td / 09.96
% last update: ow,td / 10.97

% Constants:
C0 = 10^9;		% 10^9 g s^(-2) cm^(-2) == 10^4 N cm^(-3)
h = 0.1;		% cm
rho = 1.0;		% g cm^(-3)
alpha = alfa;		% 1/cm
A = sqrt(h*C0*0.5/rho); % cm/s

time = 2/(alpha*A)*(exp(alpha*0.5*bmpos61(f))-1);	% rf. eq. (5.8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = bmpos61(f)
%BMPOS61  position of frequency f on the basilar membrane
%
%         pos = bmpos61(f)
%         returns the position on the basilar membrane (as
%         distance from stapes in cm) which represents the
%         frequency f.
%
%         Based upon: DD. Greenwood(1961): "Critical bandwidth
%                     and the frequency coordinates of the
%                     basilar membrane",
%                     J Acoust Soc Am (33), 1344-1356
%
% author/date: ow,td / 09.96
% last update: ow,td / 10.97

% Constants:
a = 1.67;		% cm
b = 0.006046;		% 1/Hz
c = 1.0;		% []
L = 3.50;		% cm
pos = L - a*log10(b*f+c);	% rf. eq. (5.6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = instphas(t,alfa)
%INSTPHAS instantaneous phase
%
%         phi = instphas(t,alfa)
%         calculates the instantaneous phase at a given time t
%         (in s) for a stimulus that should compensate the
%         spatial dispersion on the basilar membrane. It is used
%         by the function alchirp.
%
% author/date: ow,td / 09.96
% last update: ow,td / 10.97
%
% var. alpha: s. uppenkamp, 17-09-98

% Constants:
C0 = 10^9;		% 10^4 N cm^(-3)
h = 0.1;		% cm
rho = 1.0;		% g cm^(-3)
alpha = alfa;		% 1/cm
A = sqrt(h*C0*0.5/rho); % cm/s
a = 1.67;		% cm
b = 0.006046;		% 1/Hz
c = 1.0;		% []
L = 3.50;		% cm
alpha1 = 2/(alpha*A)*exp(alpha*L*0.5);	% rf. eq. (5.9)
beta1 = -alpha*a/(2*log(10));		% rf. eq. (5.9)
beta2 = 1 + 1/beta1;			% rf. eq. (5.13)
t0 = bmtime(0,alpha);			% time for 0 Hz (reference)
tau = t0 + 2/(alpha*A);			% rf. eq. (5.13)

% range check
% if (t < 0) | (t > 0.0179)
if (t < 0) | (t > 1.0179)	% achtung, range check deaktiviert, um alpha
				% groesser als 3.0 zuzulassen
  error('range error: t is only valid in [0 ; 0.0179]');
end

phi = -2*pi/b*( ((tau-t).^beta2 - tau^beta2)/(alpha1^(1/beta1)*beta2) + c*t );
				% rf. eq. (5.12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	Copyright (C) 1996	Torsten Dau, Oliver Wegner
%%				Carl von Ossietzky University Oldenburg
%%	
%%	Permission to use, copy, and distribute this software/file and its
%%	documentation for any purpose without permission by the author
%%	is strictly forbidden.
%%
%%	Permission to modify the software is granted, but not the right to
%%	distribute the modified code.
%%
%%	This software is provided "as is" without expressed or implied warranty.
%%
%%
%%	AUTHORS
%%
%%		Torsten Dau / Oliver Wegner
%%		Carl von Ossietzky University
%%		Workgroup of Medical Physics
%%		26111 Oldenburg
%%		Germany
%%
%%		e-mail:		torsten@medi.physik.uni-oldenburg.de
%%				ow@medi.physik.uni-oldenburg.de
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

