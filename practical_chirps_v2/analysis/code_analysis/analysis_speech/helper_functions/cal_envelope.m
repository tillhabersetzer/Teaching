function [envelope] = cal_envelope(audiodata,fs,b,a,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates envelope of different types:
% - 'envelope'
% - 'onset_envelope'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% envelope
envelope = abs(hilbert(audiodata));
% Zero-phase digital filtering, lowpass 32 Hz
envelope_filtered = filtfilt(b,a,envelope);

switch type
    case 'onset_envelope'
        % derivative (n-1) samples for onset envelopes
        envelope_diff = diff(envelope_filtered)/(1/fs);
        % add last sample artificially
        envelope_diff = [envelope_diff; envelope_diff(end)];
        % half-wave rectified
        envelope_diff(envelope_diff < 0) = 0;
        envelope = envelope_diff;
        
    case 'envelope'
%         envelope = 20*log10(envelope_filtered);
        envelope_filtered(envelope_filtered < 0) = 0;
        envelope_filtered = 20*log10(envelope_filtered+1);
        envelope = envelope_filtered;
    otherwise
        error('specified type is not supported!')
end
