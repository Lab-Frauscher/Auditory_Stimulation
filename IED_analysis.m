%% IED_allchannels.m
% 
% Title: Interictal Epileptiform Discharges in auditory stimulation
% Author: Barbora Matouskova
% Date: 2025-11-19
% 
% DESCRIPTION:
% This script analyzes intracranial EEG (iEEG) data to quantify the 
% response of Interictal Epileptiform Discharges (IEDs) to auditory 
% stimulation. It loads bipolar data and channel labels (SOZ/Non-SOZ, 
% Mesiotemporal/Neocortical) alongside annotations for auditory stimulation
% and evoked arousals. 
% The analysis focuses on detecting and processing IED/spike activity 
% within a window centered around each stimulation event, calculating 
% changes in IED rate and propagation across defined channel subsets 
% during Pre-Stimulation and Post-Stimulation/Arousal windows.
% 
% KEY PROCESSING STEPS:
% 1. Data segmentation and labeling into SOZ/Non-SOZ and Cortical subsets.
% 2. Spike detection (Janca detector default 60s window) centered on stimulation.
% 3. Post-processing to remove simultaneous spikes.
% 4. Quantification of normalized global IED count across all channels 
%    and each subset in Pre-Stimulation and Post-Stimulation/Arousal Window.
% 5. Calculation of mean spike propagation channel count across all channels 
%    and each subset for both windows.
% 
% -------------------------------------------------------------------------
%
% REQUIRED TOOLBOXES: 
% 
% 1. Signal Processing Toolbox
% 2. Statistics and Machine Learning Toolbox
%
% DEPENDENCIES (External Files/Functions):
% - 'spike_detector_hilbert_v25.m'
% - 'postprocessing_v3.m'
%
% MATLAB Version: R2022b or later
%
% -------------------------------------------------------------------------

% --- START OF SCRIPT ---
clear all;
clc;
close all;

addpath('functions')

% Read Annotations for Stimululation and Arousals
annots = readtable('arousal_annotations.xlsx');
patients = unique(annots.patient);

% Read Channel Subsets
subsets = readtable('channel_labels.xlsx','Sheet','Channel_labels');

arousals_GIED = table('Size', [0, 12],'VariableTypes', repmat({'double'}, 1, 12), ...
    'VariableNames', {'ID','stimulation_start','IED pre','IED post','IED pre mesiotemporal','IED post mesiotemporal','IED pre neocortical','IED post neocortical','IED pre SOZ','IED post SOZ','IED pre non-SOZ','IED post non-SOZ'});
controls_GIED = table('Size', [0, 12],'VariableTypes', repmat({'double'}, 1, 12), ...
    'VariableNames', {'ID','stimulation_start','IED pre','IED post','IED pre mesiotemporal','IED post mesiotemporal','IED pre neocortical','IED post neocortical','IED pre SOZ','IED post SOZ','IED pre non-SOZ','IED post non-SOZ'});

arousals_propagation = table('Size', [0, 12],'VariableTypes', repmat({'double'}, 1, 12), ...
    'VariableNames', {'ID','stimulation_start','IED propagation pre','IED propagation post','IED propagation pre mesiotemporal','IED propagation post mesiotemporal','IED propagation pre neocortical','IED propagation post neocortical','IED propagation pre SOZ','IED propagation post SOZ','IED propagation pre non-SOZ','IED propagation post non-SOZ'});
controls_propagation = table('Size', [0, 12],'VariableTypes', repmat({'double'}, 1, 12), ...
    'VariableNames', {'ID','stimulation_start','IED propagation pre','IED propagation post','IED propagation pre mesiotemporal','IED propagation post mesiotemporal','IED propagation pre neocortical','IED propagation post neocortical','IED propagation pre SOZ','IED propagation post SOZ','IED propagation pre non-SOZ','IED propagation post non-SOZ'});


%% Iterate over Patients
for idx = 1:length(patients)
    pIdx = patients{idx};
    fprintf('Processing %s \n', pIdx);
    
    % Read Bipolar Channels
    File_Name = ['demo_data/',pIdx,'.mat'];
    load(File_Name);
    
    % Annotatios for the Patient
    audio_stim_pos = table2array(annots(strcmp(annots.patient, pIdx),'stim_start'));
    arousal_pos = table2array(annots(strcmp(annots.patient, pIdx),'arousal_start'));

    % Channel Names
    labels = split(bp_ch_names,'-');

    % Channel Subsets
    channels_soz = subsets.Contact(find(strcmp(subsets.PatientID, pIdx) & strcmp(subsets.SOZVsNon_SOZ, 'SOZ')));
    channels_nonsoz = subsets.Contact(find(strcmp(subsets.PatientID, pIdx) & strcmp(subsets.SOZVsNon_SOZ, 'non-SOZ')));
    channels_mes = subsets.Contact(find(strcmp(subsets.PatientID, pIdx) & strcmp(subsets.mesiotemporalVsNeocortical, 'mesiotemporal')));
    channels_neo = subsets.Contact(find(strcmp(subsets.PatientID, pIdx) & strcmp(subsets.mesiotemporalVsNeocortical, 'neocortical')));
    
    % Split Bipolar Channels into Subsets
    index_soz = find(ismember(labels(:,1), channels_soz) | ismember(labels(:,2), channels_soz));
    index_nonsoz = find(ismember(labels(:,1), channels_nonsoz) & ismember(labels(:,2), channels_nonsoz));
    index_mes = find(ismember(labels(:,1), channels_mes) & ismember(labels(:,2), channels_mes));
    index_neo = find(ismember(labels(:,1), channels_neo) & ismember(labels(:,2), channels_neo));

    %% Iterate over Stimulation
    [audio_stim_pos,idx] = sort(audio_stim_pos);
    lastsize = 0;
    for stimIdx = 1:length(audio_stim_pos)
        lastsize = fprintf('Stimulation #%d/%d \n', [stimIdx length(audio_stim_pos)]);

        % Find if There Is Arousal or Not
        offset = audio_stim_pos(stimIdx);
        idx = find(arousal_pos >= offset & arousal_pos <= (offset  + 6.5*fs), 1 );
        
        % Spike Detector Uses 60 Seconds Segment
        duration = 60;

        % Calculate Start of Arousal Window/Post-Stimulus Windows
        if ~isempty(idx)
            % arousals
            arousal_index = (arousal_pos(idx) - offset)/fs;
        else
            % controls
            arousal_index = 5;  
        end

        % Prepare 60 Seconds Segment for Spike Detector
        epoch = offset - duration*fs/2;
        data_epoch = bp_data(:,epoch:epoch+duration*fs)';

        %% Spike detection
        settings = '-bl 10 -bh 60 -h 60 -jl 3.65 -dec 200'; 

        % Janca Spike Detection and Postprocessing for Simultaneous Spikes
        out = spike_detector_hilbert_v25(data_epoch,fs,settings);
        out_pp = postprocessing_v3(out,fs,size(data_epoch,2));

        %% Pre-Stimulus Window
        % Get Spikes Only 3 Seconds Before Stimulation Onset
        ch_index = [];
        sources = [];
        for i = 1:length(out_pp)
            spikes = out_pp{i}(out_pp{i} > duration/2-3 & out_pp{i} < duration/2);
            sources = [sources; spikes'];
            ch_index = [ch_index; i*ones(size(spikes'))];
        end

        % Compute Global Spike Index for Pre-Stimulus Window
        res_pre = computeGIED(sources, ch_index);

        % Compute Global Spike Index for Subsets for Pre-Stimulus Window
        if ~isempty(res_pre)
            respre_soz = res_pre(ismember(res_pre(:,1),index_soz),:);
            respre_nonsoz = res_pre(ismember(res_pre(:,1),index_nonsoz),:);
            respre_mes = res_pre(ismember(res_pre(:,1),index_mes),:);
            respre_neo = res_pre(ismember(res_pre(:,1),index_neo),:);
        else
            respre_soz = [];
            respre_nonsoz = [];
            respre_mes = [];
            respre_neo = [];
        end

        % Compute IED Propagation for Pre-Stimulus Window
        res_pre_propagation = computePIED(sources, ch_index);

        % Compute IED Propagation for Subsets for Pre-Stimulus Window
        if ~isempty(res_pre_propagation)
                respre_soz_propagation = res_pre_propagation(ismember(res_pre_propagation(:,1),index_soz),:);
                respre_nonsoz_propagation = res_pre_propagation(ismember(res_pre_propagation(:,1),index_nonsoz),:);
                respre_mes_propagation = res_pre_propagation(ismember(res_pre_propagation(:,1),index_mes),:);
                respre_neo_propagation = res_pre_propagation(ismember(res_pre_propagation(:,1),index_neo),:);
        else
                respre_soz_propagation = zeros(0,2);
                respre_nonsoz_propagation = zeros(0,2);
                respre_mes_propagation = zeros(0,2);
                respre_neo_propagation = zeros(0,2);
       end

        %% Arousal/Post-Stimulus Window
        % Get Spikes Only 3 Seconds from Arousal/Post-Stimulus Window Onset
        ch_index = [];
        sources = [];
        for i = 1:length(out_pp)
            spikes = out_pp{i}(out_pp{i} > (duration/2+arousal_index) & out_pp{i} < (duration/2+arousal_index+3));
            sources = [sources; spikes'];
            ch_index = [ch_index; i*ones(size(spikes'))];
        end

        % Compute Global Spike Index for Arousal/Post-Stimulus Window
        res_post = computeGIED(sources, ch_index);

        % Compute Global Spike Index for Subsets for Arousal/Post-Stimulus Window
        if ~isempty(res_post)
            respost_soz = res_post(ismember(res_post(:,1),index_soz),:);
            respost_nonsoz = res_post(ismember(res_post(:,1),index_nonsoz),:);
            respost_mes = res_post(ismember(res_post(:,1),index_mes),:);
            respost_neo = res_post(ismember(res_post(:,1),index_neo),:);
        else
            respost_soz = zeros(0,2);
            respost_nonsoz = zeros(0,2);
            respost_mes = zeros(0,2);
            respost_neo = zeros(0,2);
        end

        % Compute IED Propagation for Arousal/Post-Stimulus Window
        res_post_propagation = computePIED(sources, ch_index);

        % Compute IED Propagation for Subsets for Arousal/Post-Stimulus Window
        if ~isempty(res_post_propagation)
                respost_soz_propagation = res_post_propagation(ismember(res_post_propagation(:,1),index_soz),:);
                respost_nonsoz_propagation = res_post_propagation(ismember(res_post_propagation(:,1),index_nonsoz),:);
                respost_mes_propagation = res_post_propagation(ismember(res_post_propagation(:,1),index_mes),:);
                respost_neo_propagation = res_post_propagation(ismember(res_post_propagation(:,1),index_neo),:);
        else
                respost_soz_propagation = zeros(0,2);
                respost_nonsoz_propagation = zeros(0,2);
                respost_mes_propagation = zeros(0,2);
                respost_neo_propagation = zeros(0,2);
       end

        % Global Spike Index Results
        if ~isempty(idx)
            arousals_GIED = [arousals_GIED; [pIdx offset num2cell(size(res_pre,1)/length(labels)) num2cell(size(res_post,1)/length(labels)) num2cell(size(respre_mes,1)/length(index_mes)) num2cell(size(respost_mes,1)/length(index_mes)) num2cell(size(respre_neo,1)/length(index_neo)) num2cell(size(respost_neo,1)/length(index_neo)) num2cell(size(respre_soz,1)/length(index_soz)) num2cell(size(respost_soz,1)/length(index_soz)) num2cell(size(respre_nonsoz,1)/length(index_nonsoz)) num2cell(size(respost_nonsoz,1)/length(index_nonsoz))]];    
        else
            controls_GIED = [controls_GIED; [pIdx offset num2cell(size(res_pre,1)/length(labels)) num2cell(size(res_post,1)/length(labels)) num2cell(size(respre_mes,1)/length(index_mes)) num2cell(size(respost_mes,1)/length(index_mes)) num2cell(size(respre_neo,1)/length(index_neo)) num2cell(size(respost_neo,1)/length(index_neo)) num2cell(size(respre_soz,1)/length(index_soz)) num2cell(size(respost_soz,1)/length(index_soz)) num2cell(size(respre_nonsoz,1)/length(index_nonsoz)) num2cell(size(respost_nonsoz,1)/length(index_nonsoz))]];
        end

        if ~isempty(idx)
            arousals_propagation = [arousals_propagation; [pIdx offset num2cell(mean(res_pre_propagation(:,2))) num2cell(mean(res_post_propagation(:,2))) num2cell(mean(respre_mes_propagation(:,2))) num2cell(mean(respost_mes_propagation(:,2))) num2cell(mean(respre_neo_propagation(:,2))) num2cell(mean(respost_neo_propagation(:,2))) num2cell(mean(respre_soz_propagation(:,2))) num2cell(mean(respost_soz_propagation(:,2))) num2cell(mean(respre_nonsoz_propagation(:,2))) num2cell(mean(respost_nonsoz_propagation(:,2)))]];    
        else
            controls_propagation = [controls_propagation; [pIdx offset num2cell(mean(res_pre_propagation(:,2))) num2cell(mean(res_post_propagation(:,2))) num2cell(mean(respre_mes_propagation(:,2))) num2cell(mean(respost_mes_propagation(:,2))) num2cell(mean(respre_neo_propagation(:,2))) num2cell(mean(respost_neo_propagation(:,2))) num2cell(mean(respre_soz_propagation(:,2))) num2cell(mean(respost_soz_propagation(:,2))) num2cell(mean(respre_nonsoz_propagation(:,2))) num2cell(mean(respost_nonsoz_propagation(:,2)))]];
        end
    end
    fprintf('\n');
end

% Function for Global Spike Index Computation
function global_ied = computeGIED(spikes,channels)
    [sources,idx] = sort(spikes);
    ch_index = channels(idx);
    % Compute events independent from source spike (i.e., 120 ms away
    % from source spike)
    indie_spikes = [];
    ch_spikes = [];
    while ~isempty(sources)
        event = sources(1);
        % Obtain next event > 120ms from the current source 
        next_event = event+0.120;
        indie_spikes = [indie_spikes; event];
        ch_spikes = [ch_spikes; ch_index(1)];
        ch_index = ch_index(sources > next_event);
        sources = sources(sources > next_event);
    end
    global_ied = [ch_spikes indie_spikes];
end

% Function for Spike Propagation Computation
function propagating_ied = computePIED(spike,channels)
    [sources,idx] = sort(spike);
    ch_index = channels(idx);
    % Compute events independent from source spike (i.e., 120 ms away
    % from source spike)
    ied = [];
    ch_spikes = [];
    while ~isempty(sources)
        event = sources(1);
        % Obtain next event > 120ms from the current source 
        last_event = event+0.120;
        num_spikes = length(sources(sources<last_event));
        sources = sources(num_spikes+1:end);
        ch_spikes = [ch_spikes; ch_index(1)];
        ch_index = ch_index(num_spikes+1:end);
        ied = [ied; num_spikes];
    end  
    propagating_ied = [ch_spikes ied];
end

