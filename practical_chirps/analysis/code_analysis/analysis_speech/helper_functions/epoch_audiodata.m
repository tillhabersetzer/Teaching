function [audiodata_epoched] = epoch_audiodata(audiodata,fsamp,triallength,downsampling,fsamp_down)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epoch audiodata stored on path_audiodata in trials with triallength.
% It generates continuously following trials with length triallength in
% seconds. The last trial covers the rest of the data until the end.
% It offers also the possibility to downsample the trials with frequency 
% fsamp_down. Then downsampling must be set to 1 otherwise 0.
% - data_epoched is a struct containing the following:
%   - trials: cell matrix with epoched audio data
%   - sampleinfo: includes the sample number of the audiodata in the related
%     trial
%   - fsamp: is the sampling frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% option to downsample data
% downsampling = 1;

%[audiodata,fsamp] = audioread(path_audiodata);
L                 = length(audiodata);
len_epoch         = triallength*fsamp;
N                 = ceil(L/len_epoch); % number of epochs

% check if data is columnvector, than transpose it - would like to have
% rowvector with size 1xtime
if iscolumn(audiodata)
    audiodata = audiodata';
end

trlbegin = (1:len_epoch:L)';
trlend   = [trlbegin(2:end)-1; L]; 
if trlbegin(end) == L
    trlbegin(end) = [];  
end
trl      = [trlbegin trlend];

trials     = cell(1,N);
% sampleinfo = cell(N,1);

for n = 1:N
    trial = audiodata(trl(n,1):trl(n,2));
    if downsampling
        trials{1,n} = resample(trial,fsamp_down,fsamp);
    else
        trials{1,n} = trial;
    end
%     sampleinfo{n,1} = trl(n,1):trl(n,2);
end

audiodata_epoched.trials     = trials;
audiodata_epoched.sampleinfo = trl;
audiodata_epoched.fsamp      = fsamp;
audiodata_epoched.fsamp_down = fsamp_down;

end










