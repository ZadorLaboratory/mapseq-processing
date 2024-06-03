%% this code calculates re-used barcode within each brain and between each pair of brains
load('barcodematrixM253_CR.mat')       % load mapseq data

num_brain = 2;     % 10 brains in total

target_threshold = 5;       % UMI threshold for target site
injection_threshold = 30;       % UMI threshold for injection site

idx_injection = cell(1,num_brain);      % index of injection sites for each brain
idx_injection{1} = 1;
idx_injection{2} = 2;

idx_target = cell(1,num_brain);      % index of target sites for each brain
idx_target{1} = 6:39;
idx_target{2} = 40:70;


num_neuron = zeros(1,num_brain);    % now calculate # of projection neurons in each brain. criteria: max UMI in injection > injection_threshold AND max UMI in target > target_threshold
for i =1:num_brain
    num_neuron(i) = sum( max(barcodematrix(:,idx_injection{i}),[],2)>injection_threshold & max(barcodematrix(:,idx_target{i}),[],2)>target_threshold  );
end


load('barcodematrixL2.mat')  % now load viral library data
viralbarcodecount = sum(barcodematrixL2,2);       % barcode count of each barcode in the viral library

samplingpool=zeros(sum(viralbarcodecount),1,'int32');          % now simulate a viral library sampling pool; if barcode n has x counts in the library, then there are x copies of n in this sampling pool 
j=1;
for i=1:length(viralbarcodecount)
    samplingpool(j:j-1+viralbarcodecount(i))=i*ones(viralbarcodecount(i),1);
    j=j+viralbarcodecount(i);
end

% Now simulate the sampling process, it calculates the re-used barcode counts within each brain and between any two brains

num_reusedbarcodes = zeros(num_brain);  % element (i,j) in this matrix stores the number of re-used barcodes between brain i and brain j (when i=j, it stores within-brain re-used barcodes)

for i=1:(num_brain)
    for j=1:(num_brain)
        if i==j             % re-used barcodes within a single brain
            samplesize=num_neuron(i);
            a=randi([1,length(samplingpool)],samplesize,1,'int32');
            b=samplingpool(a);
            c=unique(b);
            num_reusedbarcodes(i,i) = length(b) - length(c);
        elseif i~=j         % re-used barcodes between two brains
            samplesize1=num_neuron(i);
            samplesize2=num_neuron(j);
            a1=randi([1,length(samplingpool)],samplesize1,1,'int32');
            b1=samplingpool(a1);
            c1=unique(b1);
            a2=randi([1,length(samplingpool)],samplesize2,1,'int32');
            b2=samplingpool(a2);
            c2=unique(b2);
            num_reusedbarcodes(i,j) = length(c1)+length(c2)-length(unique([c1;c2]));
            num_reusedbarcodes(j,i) = length(c1)+length(c2)-length(unique([c1;c2]));
        end
    end
end
