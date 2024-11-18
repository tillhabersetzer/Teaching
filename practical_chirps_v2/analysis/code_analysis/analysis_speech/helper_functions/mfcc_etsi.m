function [MFCCs] = mfcc_etsi(vSig,fs,Cfg)
% [MFCCs] = mfcc_etsi(vSig,fs,Cfg)
% Options: Cfg.bCalcDeltas = [0|1|2] (2)
%          Cfg.bNormFeat   = [0|1] (0)
% 
disp('Calculating MFCC ETSI features')

Cfg.init = 1;
if ~isfield(Cfg,'iCalcDeltas')
  Cfg.iCalcDeltas = 2;
end

if ~isfield(Cfg,'bNormFeat')
  Cfg.bNormFeat = 0;
end

if ~isfield(Cfg,'bUseEnergy')
  Cfg.bUseEnergy = 0;
end

[MFCCs] = signal2mfcc(vSig, fs, 0, Cfg);

% delta (and double delta) calculation
if Cfg.iCalcDeltas
  MFCCs = calc_deltas(MFCCs,Cfg.iCalcDeltas);
end

% I can also normalize features if you want me to
if Cfg.bNormFeat
  MFCCs = fe_norm_feat(MFCCs);
end
  


%--------------------------------------------------------------------------
function [MFCCs] = signal2mfcc(vSignal, fs, plotflag, Cfg)
%--------------------------------------------------------------------------
% calculate melspectrum from input signal according to ETSI standard
% ETSI ES 201 108 V1.1.2 (2000-04)
%
% MELSPEC(SIGNAL, SAMPLINGFREQ. (HZ), PLOTFLAG)
%
% Feb. 2003 / Bernd Meyer / bernd.meyer@uni-oldenburg.de
%---------------------------------------------------------------------------------

% Set parameters
iVerboseFlag     = 0  ;
% script produces same output as the nokia frontend (delivered with Aurora2
% CDs) when no pre-emphase is used
iUsePreemphase   = 0;
bUseEnergy = Cfg.bUseEnergy;

if fs == 8000
  % default values for 8 kHz input signals
  nfft = 256;                 % number of samples for fft
  WindowLength = 200;         % length of window function
  nW = hamming(WindowLength); % type of window function
  ShiftInt = 80;              % shift interval is 80 samples for 8 kHz
elseif fs == 11000
  % default values for 11 kHz input signals
  nfft = 256;                 % number of samples for fft
  WindowLength = 256;         % length of window function
  nW = hamming(WindowLength); % type of window function
  ShiftInt = 110;             % shift interval is 80 samples for 8 kHz
elseif fs == 16000
  % default values for 16 kHz input signals
  nfft = 512;                 % number of samples for fft
  WindowLength = 400;         % length of window function
  nW = hamming(WindowLength); % type of window function
  ShiftInt = 160;             % shift interval is 80 samples for 8 kHz
else
  error(['sampling frequency of ' num2str(fs) ' Hz is not supported']);
end

%---------------------------------------------------------------------------------

if iVerboseFlag
  disp('Compute melspec with the following parameters: ');
  disp(['nfft = ' num2str(nfft) ', WindowLength = ' num2str(WindowLength) ...
    ' ShiftInt = ' num2str(ShiftInt) ' fs = ' num2str(fs)]);
end

nRows = size(vSignal,1);
nCols = size(vSignal,2);
if (nRows > 1) & (nCols > 1)
  error('only mono input supported')
elseif (nCols > 1)
  vInput = vInput';
end

% starting frequency of the filterbank
f_start = 64;     % Hz

% scale input signal
% use 1.5 as threshold, in case amplitude gets a littler larger due to
% resampling or so.
if abs(max(vSignal)) < 1.5
  vSignal = vSignal * 2^15;
end

% offset compensation
vSignal_of = zeros(size(vSignal));
vSignal_of(1) = 0; % is this correct?
vSignal_of(2:length(vSignal)) = vSignal(2:length(vSignal)) - vSignal(2-1:length(vSignal)-1) + 0.999*vSignal_of(2-1:length(vSignal)-1);
vSignal = vSignal_of;
clear vSignal_of;

% framing: vSignalFrame(:,i) are all samples in ith time frame
i = 0;
for k = 1:ShiftInt:(length(vSignal)-WindowLength+1)
  i = i+1;
  vSignalFrame(:,i) = vSignal(k:k+WindowLength-1);
end
NumberOfFrames = i;

if iUsePreemphase
  for i=1:NumberOfFrames
    vSignalFrame(2:WindowLength,i) = vSignalFrame(2:WindowLength,i) - 0.97*vSignalFrame(2-1:WindowLength-1,i);
  end
end % of if

% windowing & energy calculation
vLogEnergy = zeros(1,NumberOfFrames);
for i=1:NumberOfFrames
  	
  vSignalFrame(:,i) = vSignalFrame(:,i) .* hamming(WindowLength);
end

% fft
for i=1:NumberOfFrames
  spec_frame(:,i) = abs(fft(vSignalFrame(:,i),nfft));
end

% look only at freq bins 1..nfft/2+1
% bin(k) in etsi standard means bin(k+1) here!
bin = zeros(nfft/2+1,NumberOfFrames);

for i=1:NumberOfFrames
  bin(1:nfft/2+1,i) = spec_frame(1:nfft/2+1,i);
end
clear spec_frame;

% fc(i) is the center freq. of the ith mel channel
% fc(1) ~ 124 Hz, fc(23) ~ 3660 Hz

fc = eval_fc(fs,f_start);

% cbin(i) = center freq of the ith mel channel in terms of fft bin indices
% cbin(k) in ETSI paper means cbin(k+1) here!
cbin(25) = nfft/2;
cbin(1)  = round((f_start*nfft)/fs);
cbin(2:24) = round((fc*nfft)/fs);

% evaluate weighted sum of fft-bins in each band; use triangular, half-overlapping windows
[MFCCs] = mfcc_eval(cbin,bin,fs,nfft,f_start,NumberOfFrames);

if bUseEnergy
  MFCCs(end,:) = vLogEnergy;
end
  
  
if iVerboseFlag
  disp(['Processed ' num2str(size(MFCCs,2)) ' frames.']);
end

% plot compressed melspec
if plotflag
  openwindow(plotflag,[650 340 330 180],'Compressed Mel Spectrogram');
  % h=findobj('Type','axes','Tag','MelspecWindow');
  % axes(h);
  imagesc(melspec);
  axis xy; % invert freq axis
  colormap jet;
  cb=colorbar;
  ylabel('Frequency bin');
  xlabel('Time frame');
end

%---------------------------------------------------------------------------------
function [vMFCC] = mfcc_eval(cbin,bin,fs,nfft,f_start,NumberOfFrames)
% evaluate compressed filterbank sum
for j=1:NumberOfFrames % index for time frames
  for k=1:23 % index for frequency bins

    i = cbin(k-1+1):cbin(k+1);
    sum1 = sum((i - cbin(k-1+1)+1)*bin(i+1,j) / (cbin(k+1) - cbin(k-1+1) +1));

    i = (cbin(k+1)+1):cbin(k+1+1);
    sum2 = sum(  [1 - ((i-cbin(k+1)) / (cbin(k+1+1)-cbin(k+1)+1))  ]*bin(i+1,j));

    sum3 = sum1 + sum2;

    % flooring
    if (sum3 < 2*10^(-22))
      % disp('flooring');
      sum3 = 2*10^(-22);
    end

    % non linear transformation
    fbank_comp(k,j) = log(sum3);

  end

end

for j=1:NumberOfFrames % index for time frames
  for k=1:23 % index for frequency bins
    % flooring. filter bank output may not be < 50
    fbank_comp = max(fbank_comp,-50);
    % cepstral coefficients (4.2.11)

    for i1=0:12
      sum4 = 0;
      for j1 = 1:23
        sum4 = (fbank_comp(j1,j) * cos((pi*i1)/23 * (j1 - 0.5))) + sum4;
      end % of for
      vMFCC(i1+1,j) = sum4;
    end % of outer for
  end
end

% Nokia FE puts the c0 coeffient
vMFCC = [vMFCC(2:13,:); vMFCC(1,:)];

%--------------------------------------------------------------------------
function [fc] = eval_fc(fs,f_start)
% calculate center frequencies for all 23 bins
i=1:23;
fc = mel_inv(mel(f_start) + (mel(fs/2)-mel(f_start))*i / (23 + 1));

%--------------------------------------------------------------------------
function [y] = mel(x)
% calculate melfrequency y for input freq. x
y  = 2595*log10(1 + x./700);

%--------------------------------------------------------------------------
function [y] = mel_inv(x)
% calculate inverse melfrequency y for input freq. x
y = (10.^(x./2595)-1)*700;


