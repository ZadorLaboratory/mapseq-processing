load('barcodematrixL1ZL235_Mseq130_SL.mat');    % Load L1 matrix
load('barcodematrixZL235_Mseq130_SL.mat');      % load L2 matrix
load('spikesZL235_Mseq130_SL.mat');  	% load spike in
target_L2 = 1:40; 						%change target site info
target_L1 = 44; 						%change L1 target site info
num_spikein = zeros(1,length(spikes)); 	%for too many spike-in counts

for i=1:length(spikes)
    num_spikein(i) = length(spikes(i).counts2u);
end

num_target_L2 = sum(sum(barcodematrix(:,target_L2)));          % total num of L2 molecules in L2 targets
num_target_L1 = sum(sum(barcodematrixL1(:,target_L1)));        % total num of L1 molecules in L1 targets
num_spikes_L1 = sum(num_spikein(target_L1));                   % total num of spike in molecules in L1 targets
num_spikes_L2 = sum(num_spikein(target_L2));                   % total num of spike in molecules in L2 targets
num_templateswitching = sum(sum(barcodematrix(:,target_L1)));  % num of L2 molecules detected in L1 targets
c = num_templateswitching / (num_target_L2 * (num_spikes_L1 + num_target_L1));   % template switching coefficient

ratio_ts = 0.5 * c * num_target_L2+c * num_spikes_L2;         % ratio of template swtiching molecules in all L2 targets

