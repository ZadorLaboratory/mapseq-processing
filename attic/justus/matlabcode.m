#%find connected components
#graph=[];
#for i=1:length(bcnfilt)
#    %find the connected graph components
#    [graph(i).S, graph(i).G]= graphconncomp( clustermatrix1(i).C, 'Directed', 'false'); 
#end
#% save('graph1.mat','graph'); 
# https://stackoverflow.com/questions/53277121/dulmage-mendelsohn-matrix-decomposition-in-python

#
# /grid/zador/data_nlsas_norepl/MAPseq/processing/MATLAB common/MATLAB common/charlienash-nricp-5d1cb79/dependencies/geom3d-2017.12.01/geom3d/meshes3d
#
#function [S,C] = conncomp(G)
#% CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
#% toobox. G is an n by n adjacency matrix, then this identifies the S
#% connected components C. This is also an order of magnitude faster.
#%
#% [S,C] = conncomp(G)
#%
#% Inputs:
#%   G  n by n adjacency matrix
#% Outputs:
#%   S  scalar number of connected components
#%   C

#% Transpose to match graphconncomp
#G = G';

#[p,~,r] = dmperm(G+speye(size(G)));
#S = numel(r)-1;
#C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
#C(p) = C;
#end

#%collapse barcodes to most abundant member of the connected graph component
#  
#for i=1:length(bcnfilt)
#    x=1:graph(i).S;
#    [tf,loc]=ismember(x,graph(i).G,'R2012a');
#    collapsedreads=data(i).reads(loc,:);
#    collapsedcounts=accumarray(graph(i).G',data(i).counts);%'
#    [corrected(i).counts2u,ix]=sort(collapsedcounts,'descend');
#    corrected(i).reads2u=collapsedreads(ix,:);
#end