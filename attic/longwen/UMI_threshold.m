threshold_UMI = 0:10;       	% test a variety of UMI threshold_UMI
threshold_injection = 0;       	% threshold for injection sites
idx_injection = [];        		% SSI for injection sites
idx_target = 1:108;          	% SSI for target sites (including negative control)
idx_negative_ctrl = [106,108];	% SSI for ctrl
error_rate_false_positive = zeros(1,length(threshold_UMI));

for i=1:length(threshold_UMI)       % calculate the error rate for each UMI threshold
    num_false_positive = sum( max(barcodematrix(:,idx_injection),[],2)>threshold_injection 
    							& 
    						 max(barcodematrix(:,idx_negative_ctrl),[],2)>threshold_UMI(i) ); 
    % # of neurons with false positive projection: max injection UMI>50 AND max ctrl UMI > threshold
    
    num_total = sum( max(barcodematrix(:,idx_injection),[],2)>threshold_injection 
    				 & 
    				 max(barcodematrix(:,idx_target),[],2)>threshold_UMI(i) ); 
    % # of projection neurons: max injection UMI>50 AND max target UMI > threshold

    error_rate_false_positive(i) = num_false_positive/num_total;
end
