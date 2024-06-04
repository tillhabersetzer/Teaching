close all
clear all
clc

[up, fs1] = audioread("up_rmseq.wav");
[down, fs2] = audioread("down_rmseq.wav");
[click, fs3] = audioread("click_rmseq.wav");

% [up, fs1] = audioread("up.wav");
% [down, fs2] = audioread("down.wav");
% [click, fs3] = audioread("click.wav");

rms_up    = sqrt(mean(up.^2));
rms_down  = sqrt(mean(down.^2));
rms_click = sqrt(mean(click.^2));

% Ddjust click to same length as up
click_adjusted = [click;zeros(length(up)-length(click),1)];
rms_click_adjusted = sqrt(mean(click_adjusted.^2));


