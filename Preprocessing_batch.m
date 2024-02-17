%! -*- coding:utf-8 -*-

%%
close all; clear; clc;

% 定义数据所在的路径
input_path = 'E:\EEGDATA\PYTHON_EEG\Rawset\';
% Save path
output_path = 'E:\EEGDATA\PYTHON_EEG\Preprocset\';

% list data files
files = dir([input_path '*.set']); 
filenames = {files.name}';

COUNTFILE = size(filenames,1);

Channelz_value_automatic_detection = 3.29;
Number_of_EEG_electrodes = 60; % without reference and ground electrodes

for VP = 1:COUNTFILE

    sname = erase(filenames{VP},'.set');

    % Start EEGLAB
    [ALLEEG, EEG, CURRENTSET] = eeglab;

    %LOAD THE DATA, all channels, all timepoints
    EEG = pop_loadset('filename',filenames{VP},'filepath',input_path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',sname,'gui','off');

    %% Step 1
    % Channel locations
    EEG = pop_chanedit(EEG, 'lookup', 'D:\MATLABtools\eeglab2023.1\plugins\dipfit\standard_BEM\elec\standard_1005.elc');
    EEG = eeg_checkset( EEG );

    % Remove unwanted channels
    EEG = pop_select( EEG,'rmchannel',{'VEO','HEO','M2','CB1', 'CB2'}); 
    EEG = eeg_checkset( EEG );

    % Rereference
    EEG = pop_reref( EEG, [] );
    EEG = eeg_checkset( EEG );

    % save the data before remove the bad channels we have detected before
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filenames{VP},'gui','off');
    % EEG = pop_saveset( EEG, 'filename', filenames{VP}, 'filepath', step1_path);
    EEG = eeg_checkset( EEG );

    %% Step 2
    % Automatic rejection of bad channels
    % THREE DIFFERENT BAD ELECTRODE CHECKS (Probability, kurtosis, Frequency 1-125 Hz): 
    % load the data 3 times for checks, then a last time to carry on with the interpolation of the bad channels
    % look for bad channels
    [~, indelec1] = pop_rejchan(EEG, 'elec',1:Number_of_EEG_electrodes ,'threshold',Channelz_value_automatic_detection,'norm','on','measure','prob'); 	%we look for probability
	[~, indelec2] = pop_rejchan(EEG, 'elec',1:Number_of_EEG_electrodes ,'threshold',Channelz_value_automatic_detection,'norm','on','measure','kurt');	%we look for kurtosis 
    [~, indelec3] = pop_rejchan(EEG, 'elec',1:Number_of_EEG_electrodes ,'threshold',Channelz_value_automatic_detection,'norm','on','measure','spec','freqrange',[1 125] );	%we look for frequency spectra

    % look whether a channel is bad in multiple criteria
    index=sort(unique([indelec1,indelec2,indelec3])); %index is the bad channel array

    for i = 1:size(index,2)
        sub_indexarray(VP,i) = index(1,i);
    end

    % remove channels because of index array
	% Interpolate Channels (Bad Channels)
    EEG = pop_interp(EEG, index, 'spherical');
    EEG = eeg_checkset( EEG );

    % Repair bursts and reject bad portions of data
    EEG = clean_artifacts( EEG, 'Highpass', 'off',...
        'ChannelCriterion', 'off',...
        'LineNoiseCriterion', 'off',...
        'BurstCriterion', 30,...
        'WindowCriterion',0.3);

    % save the data after removed the bad channels
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname', filenames{VP}, 'gui','off'); 
    % EEG = pop_saveset( EEG, 'filename', filenames{VP}, 'filepath', step2_path);
    EEG = eeg_checkset( EEG );

    %% Step 3
    % Filtering
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1); %高通滤波
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 60); %低通滤波
    EEG =  pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1); %凹陷滤波
    EEG = eeg_checkset( EEG );

    % Downsampling to 250 Hz
    EEG = pop_resample(EEG, 250);
    EEG = eeg_checkset( EEG );

    % epoch and baseline correction
    EEG = pop_epoch( EEG, {'11'  '12'  '21' '22'  '31'  '32'}, [-0.2  0.8]);
    EEG = pop_rmbase( EEG, [-200 0] ,[]);
    EEG = eeg_checkset( EEG );

     % save data
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'overwrite','on','gui','off');
    % EEG = pop_saveset( EEG, 'filename', filenames{VP}, 'filepath', step3_path);
    EEG = eeg_checkset( EEG );

    %% Step 4
    %run ICA
     EEG = pop_runica(EEG, 'extended',1,'interupt','on', 'pca', Number_of_EEG_electrodes-size(index,2)); % runICA
     EEG = pop_iclabel(EEG, 'default'); % ICLabel
     EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]); % THESE ARE EXAMPLE THRESHOLDS ! VALIDATE ON YOUR DATA !
     for i = 1:size(EEG.etc.ic_classification.ICLabel.classifications,1)
         if EEG.etc.ic_classification.ICLabel.classifications(i,1)>max(EEG.etc.ic_classification.ICLabel.classifications(i,2:6))% if signal probability is higher than "pure" artifact
             classifyICs(i)=0;
         else
             classifyICs(i)=1;
         end
     end
    EEG.reject.gcompreject=classifyICs; %
    EEG = eeg_checkset( EEG );

    % Save data after ICA
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname',filenames{VP},'overwrite','on','gui','off');
    % EEG = pop_saveset( EEG, 'filename', filenames{VP}, 'filepath', step4_path);
    EEG = eeg_checkset( EEG ); 

    %% Step 5
    % Clean data by rejecting epochs
    [EEG, rejindx] = pop_eegthresh(EEG,1,1:EEG.nbchan ,-100,100,EEG.xmin, EEG.xmax,0,1);
    % Save data
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',sname,'overwrite','on','gui','off');
    EEG = pop_saveset( EEG, 'filename', filenames{VP}, 'filepath', output_path);
    EEG = eeg_checkset( EEG ); 

end

%% End