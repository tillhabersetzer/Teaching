function [bitmask, wordlength, bitmaskstr, data] = int2bitmask(data_, wordlength_)
% [bitmask, wordlength] = int2bitmask(data[,wordlength])
%
% int2bitmask - converts a matrix of positive integers into the
% corresponding bitmasks - as logical array 'bitmask' as well as a
% string 'bitmaskstr'
%
% Inputs:
% data       - integer or array of integers 
%             (floating point numbers will be rounded, maximum: flintmax )  
% wordlength - minimum length of bitmask that should be considered.
%              wordlength will be increase automaticalla y when needed.
%              (optional, default: 24)
%
% Outputs:
% bitmask    - is the resulting bitmask as a vector or array (NUM x BITS) of 
%              logicals with NUM being the number of given integers and
%              BITS the number of bits 
% wordlength - might differ from the input 'wordlength' if it was extended
%              automatically to allow the biggest value of 'data'
% bitmaskstr - is the resulting bitmask as a characterarray
% data       - input'data' casted to datatype 'uint64'
%
% see also: EEGTrigWord2info, EEGTrigID2info
%
% (c) Manfred Mauermann 2017
%

%last changes:
%25.08 - completely rewritten by using the Matlab command dec2bin 
%

if nargin < 2,   wordlength_ = []; end
if nargin < 1,   help(mfilename); return ; end
if isempty(wordlength_), wordlength_ = 24; end

if any(data_ > flintmax)
    error('maximum allowed value for data is %d (2^%d+1)',flintmax,log2(flintmax-1));
end
data       = cast(data_,'uint64');
bitmaskstr = dec2bin(data,wordlength_);
bitmask    = logical(bitmaskstr-'0');
wordlength = size(bitmask,2)  ;

if wordlength ~= wordlength_
    warning('int2bitmask: wordlength is extended from %d to %d to allow maximum values',wordlength_,wordlength);
end
    