threshold_UMI = 0:10;           % test a variety of UMI threshold_UMI
threshold_injection = 50;       % threshold for injection sites
idx_injection = 1:4;            % SSI for injection sites
idx_target = 5:20;              % SSI for target sites (including negative control)
idx_negative_ctrl = [19,20];    % SSI for ctrl

error_rate_false_positive = zeros(1,length(threshold_UMI));

for i=1:length(threshold_UMI)      
	% calculate the error rate for each UMI threshold
    % # of neurons with false positive projection:  max injection UMI>50 AND max ctrl UMI > threshold
    num_false_positive = 
    	sum( 
    		max( barcodematrix(:,idx_injection) ,[] , 2) > threshold_injection 
    		& 
    		max( barcodematrix(:,idx_negative_ctrl),[],2) > threshold_UMI(i) 
    	); 
    
    % # of projection neurons: max injection UMI > 50 AND max target UMI > threshold    
    num_total = 
    	sum( 
    	    max( barcodematrix(:,idx_injection), [] , 2) > threshold_injection 
    	    & 
    	    max( barcodematrix(:,idx_target),[],2) > threshold_UMI(i) 
    	); 

    error_rate_false_positive(i) = num_false_positive / num_total;
end
