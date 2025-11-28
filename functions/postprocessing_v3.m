function out_ch=postprocessing_v3(out,fs,num_chans)

% POST PROCESSING JANCA DETECTOR DETECTIONS
% @========================================================================
% Authors: Sayeed, John Thomas
% DATE CREATED: 23-April-2021
% DATE MODIFIED: 22-Spetember-2022
% Wrote a completely new code in a new way. Here, we first remove spikes
% that occur in atleast 50% of the channels with an error of 50 ms as there
% could be a sampling error of 10 samples. 

% FOR RESEARCH PURPOSES ONLY. KINDLY CONTACT THE AUTHOR FOR REUSE AND BUGS.
% ========================================================================@
% Tasks:
% 1. Intput: 
     %out: Output from the Janca detector.
     %fs: Sampling frequency.
     %num_chans: Number of channels from montage.txt.
% 2. Remove the codections occuring simulatneously in atleast 50% of the
     %channels.
% 3. Remove all the detections that occur within 300 milliseonds defined by
     % 'threshold' within asingle channel.
%4. Return 'out_ch': each cell containing detections in a single channel.
%==========================================================================&
threshold=0.120;% Threshold in seconds
threshold_artifact=0.05; % Threshold in seconds

spike=out.pos;
ch_list=out.chan;

% Execute if there are more than one detections
if length(spike)>1
    
    %% Remove codetections on more than 50% of channels   
    % For each spike, find if there are multiple detections within 50 ms on
    % alteast 50% of the channels. If so, remove them.
    
    for spike_no=1:length(spike)
        spike_test=spike(spike_no);
        idx=find(abs(spike-spike_test)<threshold_artifact);
        if length(idx)>num_chans/2
            spike(idx)=0;
        end
    end
     
    ch_list(find(spike==0))=[];
    spike(find(spike==0))=[];
        
    %% Arrange spike output into channels, remove spikes within a window of
    % 300 ms to account for spindles misdetections
    out_ch=[];
    for i=1:num_chans
        spike_ch=sort(spike(find(ch_list==i)));
        for j=1:length(spike_ch)
            spike_test=spike_ch(j);
            idx=find(abs(spike_ch-spike_test)<threshold);
            if length(idx)>1

                idx_all=[];
                for k=1:length(idx)
                    idx_all=[idx_all; find(abs(spike_ch-spike_ch(idx(k)))<threshold)];
                end
                idx_all=unique(idx_all);idx_all2=[];
                for kk=1:length(idx_all)
                    idx_all2=[idx_all2; find(abs(spike_ch-spike_ch(idx_all(kk)))<threshold)];
                end
                 idx_all2=unique(idx_all2);idx_all3=[];
                for kkk=1:length(idx_all2)
                    idx_all3=[idx_all3; find(abs(spike_ch-spike_ch(idx_all2(kkk)))<threshold)];
                end


                spike_ch(unique(idx_all3))=0;
            end
        end
        
        spike_ch(find(spike_ch==0))=[];
        out_ch{i}=spike_ch'; clear spike_ch
    end

else
    out_ch{out.chan}=out.pos;
end