function [power,time] = calculate_power(vSig,fs,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data trialwise?
% shift: 10 ms
% window length: 25ms
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create an input parser object
p = inputParser;
% add inputs to the scheme
defaultStyle = 'nolog';
validStyle   = {'log','nolog'};
checkStyle   = @(x) any(validatestring(x,validStyle));

addRequired(p,'vSig',@isnumeric);
addRequired(p,'fs',@isnumeric);
addOptional(p,'Style',defaultStyle,checkStyle);
% parse inputs
parse(p,vSig,fs,varargin{:});

% get access to input
% vSig  = p.Results.vSig;
% fs    = p.Results.fs;
% Style = p.Results.Style;

% settings
%--------------------------------------------------------------------------
winLen = 25/1000;            % 25ms
winLen = round(winLen * fs); % samples
shift  = 10/1000;            % 10ms
shift  = round(shift * fs);  % samples

% size of signal
[num_chan, len_vSig] = size(vSig);

% framing
frames     = 1:shift:(len_vSig-winLen+1); % cut off rest of signal
num_frames = length(frames);
% frame-matrix: channel x time x frames (maybe you need this one later)
% vSignalFrame = zeros(num_chan,winLen,num_frames);

% power-matrix: channel x frames (frame is explizit time)
power = zeros(num_chan,num_frames);
time  = zeros(1,num_frames);

for i = 1:num_frames
  % startsample for epoch
  n = frames(i);
  % vSignalFrame(:,:,n) = vSig(:,n:n+winLen-1);
  power(:,i) = sum(vSig(:,n:n+winLen-1).^2,2);
  time(1,i)  = median((n:n+winLen-1))/fs; % in sec
end

% check log argument
if strcmp('log',p.Results.Style)
    power = log10(power);
end


end