%%       
%     This code is for testing QRSdetection.m
%     
%     MIT-BIH Arrhythmia DB No. 100 is used as the sample ECG data(100.mat)
%     100.mat includes 2-channel waveforms(variable: M)
%     and time information (variable: TIME)
%     Input arguments of QRSDetection.m are the single channel waveform, 
%     sampling frequency(Fs),
%     time information(TIME) and filter type(MIT-BIH or AHA)
%
%     This source code is for the following article :
%     Jinkwon Kim and Hangsik Shin, "Simple and Robust Realtime QRS
%     Detection Algorithm based on Spatiotemporal Characteristic of the
%     QRS Complex", Plos One, 2016
%     article link : TBD
%
%     Jinkwon Kim and Hangsik Shin (hangsik.shin@gmail.com)
%     @ Healthcare Solution Laboratory, Chonnam National University
%     2015.12.31  
%

close all; clc; clear all;

load('100.mat') % MIT-BIH Record No. 100
lead = 1;       % Channel 1
Fs=360;         % Sampling Frequency
ftype=0;        % Filter for MIT-BIH
[loc, time] = QRSdetection(M(:,lead),Fs,TIME,ftype);

figure('Name','QRS Detection');
plot(TIME, M(:,lead), time, M(loc,lead),'ro')
xlabel 'Time(s)'
ylabel 'Amplitude(V)'

% End of the Code