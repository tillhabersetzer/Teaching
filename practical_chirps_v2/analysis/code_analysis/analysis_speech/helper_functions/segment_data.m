function [data_segments,tab_label_segments] = segment_data(data,winlen,shift,fsamp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cuts the data in each trial cell into segments of desired length (winlen)
% and desired shift. It treats is trial separately ans returns also a label
% vector describing which segments relates to which trials
%
% data: cell-array with multiple trials (like trial field in fieldtrip)
% winlen: desired window length for segments in seconds
% shift: shift between following windows
% fsamp: sampling frequency of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_epochs = length(data);

WindowLength = winlen*fsamp; % in samples
ShiftInt     = shift*fsamp;  % in samples

data_segments  = [];
label_segments = zeros(N_epochs,2);
epoch_name     = [];
n              = 1;

for e = 1:N_epochs
    
    label_segments(e,1) = n;
    
    for s = 1:ShiftInt:(length(data{e})-WindowLength+1)
      data_segments{n} = data{e}(:,s:s+WindowLength-1);
      n = n+1;
    end  
    
    add = {['epoch ',num2str(e)]};
    epoch_name = cat(1,epoch_name,add);
    
    label_segments(e,2) = n-1;
end

tab_label_segments = array2table(label_segments,...
                     'VariableNames',{'Start','End'},...
                     'Rownames',epoch_name);
end
