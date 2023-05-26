
function [B,Bseq,Bnorm] = normBCmat(barcodematrix,refbarcodes,spikes,sourcethresh,projthresh,sourcesite,projsite)

%Filter barcodes so that all remaining barcodes have counts >sourccethresh
%at sourcesite, and max(counts) at other sites > projthresh. then normalize
%barcodes by spikes. sourcethresh must be equal or larger than projthresh.

if ~exist('sourcesite','var')
    sourcesite = [];
    if ~exist('projsite','var')
        projsite = 1:size(barcodematrix,2);
    end
else
    if ~exist('projsite','var')
        projsite = ~ismember(1:size(barcodematrix,2),sourcesite);
    end
end


%convert spike counts
x=zeros(1,length(spikes));
for i=1:length(spikes)
    x(i)=size(spikes(i).counts2u,1);
end


%filter barcodes
if ~isempty(sourcesite)
    B=barcodematrix(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
    Bseq=refbarcodes(max(barcodematrix(:,projsite),[],2)>projthresh & max(barcodematrix(:,sourcesite),[],2)>sourcethresh,:);
else
    B=barcodematrix(max(barcodematrix,[],2)>projthresh,:);
    Bseq=refbarcodes(max(barcodematrix,[],2)>projthresh,:);
end


%normalize barcode matrix
Bnorm=B./repmat(x,size(B,1),1);
save('filtBCmat.mat','B','Bseq','Bnorm');
