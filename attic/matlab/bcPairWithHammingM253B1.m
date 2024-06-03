% bcPairWithinHammingD
% finding pairs of BC within max hamming disatance
% Speed estimation: single CPU (i7-8750H), 448385 bc, 1680 s
%
% 03042024 LY

% Settings
% Input matrix
bc = B1seq;

% (LY: I just like to represent gtac as 1234, personal preference)
gtac = [71 84 65 67];

% Max hamming distance
maxHamming = 3;

% Block size for each loop for hamming distance computation
% Please adjust the block number to fit your CPU and RAM
blockSize = 500;

%% Find pairs < maxHamming
% 1. Split data into blocks, calculate hamming distance using blocks
% 2. Because the goal is to find the pairs with < maxHamming, so exclude
% the previous blocks to cut down ~50% computation

% Only include the first 30 nt
bc = bc(:,1:30);

% Change gtac to 1-4 (LY: I am more used to this...)
[~,bc] = ismember(bc,gtac);

% Use uint8 to decrease the size
% bc2 is bc, but swap the 1st and 3rd dimesnion, so rows now is the 3rd-D
bc = uint8(bc);
bc2 = permute(bc,[3 2 1]);

% Number of BC
n = size(bc,1);

% Block edges
block = 0:blockSize:n;
block = [block,n];
block = unique(block);

I = {};
for i = 2:numel(block) % parfor, but single CPU has similar speed
    tic;

    % Prepare blocks of barcodes for comparison ---------------------------

    % Current block 
    iBlock = (block(i-1)+1):block(i);
    iBlock = reshape(iBlock,[],1);

    % To speed up, exclude the IDs before the current block
    iBlock2 = iBlock(1):n;
    iBlock2 = reshape(iBlock2,[],1);

    iBC = bc(iBlock,:);
    iBC2 = bc2(:,:,iBlock2);

    % Hamming distance ----------------------------------------------------

    % logical: nt is different
    iD = iBC ~= iBC2;
    % sum it as hamming distance
    % (nt is the 2nd dimension)
    iD = sum(iD,2);
    % (The max of uint8 is 255)
    iD = uint8(iD);
    % Squeeze for using the find function
    iD = squeeze(iD);

    TF = iD < maxHamming;

    % Get IDs -------------------------------------------------------------
    % row: its ID; col: its partners ID
    [row,col] = find(TF);
    % Convert back to the original ID
    row = iBlock(row);
    col = iBlock2(col);

    % Exclud itself
    TF = row ~= col;    

    I{i,1} = [row(TF),col(TF)];

    disp(['running block ',num2str(i),'; ',num2str(toc)]);
end

% Get unique pairs
I = vertcat(I{:});
I = sort(I,2);
I = unique(I,'rows');

%% Barcode matrix comparison

% BC comparison
bcCmp = {};
for i = 1:size(I,1)
    iI = I(i,:)';
    bcCmp{i,1} = B1seq(iI,:);
end

% Data comparison ---------------------------------------------------------
% Getting the corresponding rows of B1
bcMatCmp = {};
for i = 1:size(I,1)
    iI = I(i,:)';
    bcMatCmp{i,1} = B1(iI,:);
end

% Reporting counts
countCmp = cellfun(@(X) sum(X,2),bcMatCmp,'UniformOutput',false);