function [bundleLibs,avg_M,PPidx,EM_pix,IDX_comp] = extract_bundles_by_angle(data,M0,tresh,flag_Npx,weightEuc)
% ========================================================================= 
% Extract endmember bundles from the HSI by selecting pixels similar to
% reference endmembers (i.e. pure)
% 
% INPUT:  data      - HS image (nr * nc * L)
%         M0        - Matrix containing reference signatures of each endmember (L * P)
%         tresh     - Pixels closer than thresh to the reference endmembers are
%                     considered pure (can be scalar or 1 * P vector)
%         flag_Npx  - If true, "thresh" will contain the number of purest
%                     pixels of  each class to extract
%         weightEuc - also considers the Eclidean distance, weighted by
%                     this factor (i.e. dist = cosdist + weightEuc*eucdist)
% 
% optional arg pairs: 'M true',Mth - true endmember matrix to order the
%                                    estimates against
% 
% OUTPUT: bundleLibs - extracted endmembers for each class
%         avg_M      - Average endmembers from the extracted pixels
%         PPidx      - Indices of pure pixels of each class (between 1...N)
%         EM_pix     - Pure pixels of each class in cube format (remaining positions are zeros)
%         IDX_comp   - Indicator function of pure pixels positions of each class
% 
% 
% Author: Ricardo Borsoi
% 04/10/2018
% =========================================================================

if nargin < 4
    flag_Npx = false;
end
if nargin < 5
    weightEuc = 0;
end


[L,P] = size(M0);
[nr,nc,L] = size(data);
N = nr*nc;

% check consistency of input arguments
if ndims(data) ~= 3, error('HS image must be ordered as a cube!'), end
if size(data,3) ~= L, error('Image and endmember matrix have different number of bands!'), end
if numel(tresh) > 1 && numel(tresh) ~= P, error('Set of thresholds or numbers of pure pixels do not match the number of endmembers!'), end

if numel(tresh) == 1
    tresh = tresh * ones(1,P);
end

% Reorder data as a matrix with column-wise lexicographic ordering
data_r = reshape(data,nr*nc,L)';


% compute the distances from each pixel to M0
cosdists = 1 - M0' * data_r ./ sqrt(repmat(sum(data_r.^2,1),[P,1]) .* repmat(sum(M0.^2,1)',[1,N]));
% eucdists = sqrt(permute(sum(abs(repmat(data_r,[1 1 P]) - permute(repmat(M0,[1 1 N]),[1 3 2])).^2, 1), [3 2 1])/L);
eucdists = (permute(sum(abs(repmat(data_r,[1 1 P]) - permute(repmat(M0,[1 1 N]),[1 3 2])).^2, 1), [3 2 1])/L);


% weight both distances
cosdists = cosdists + weightEuc * eucdists;

% Get indices of pixels similar to pure pixels (do this per EM)
p_pixels = zeros(size(cosdists));

% Extract pure pixels (detect which pixels are pure)
if flag_Npx == false
    % extract pure pixels by threshold
    for i=1:P
        p_pixels(i,:) = (cosdists(i,:) < tresh(i));
    end
else
    % extract Nk purest pixels
    for i=1:P
        [~,Idsort] = sort(cosdists(i,:));
        p_pixels(i, Idsort(1:tresh(i))) = 1;
    end
end




% Extract pure pixels from the image
EM_pix = cell(1,P);      % pure pixels
IDX_comp = cell(1,P);    % pure pixels indices (indicator function)
PPidx = cell(1,P);       % linear indices of pure pixels
avg_M = zeros(size(M0)); % average endmember signatures
bundleLibs = cell(1,P);  % Extracted signatures per class

for i=1:P
    EM_pix{i} = zeros(L,N);
    IDX_comp{i} = zeros(L,N);

    % get indices
    idx = find(p_pixels(i,:));
    PPidx{i} = idx;
    
    % assign values
    EM_pix{i}(:,idx) = data_r(:,idx);
    IDX_comp{i}(:,idx) = ones(L,length(idx));
    
    % Store bundles
    bundleLibs{i} = data_r(:,idx);
    
    % Compute mean EM signature
    avg_M(:,i) = mean(data_r(:,idx),2);
end


% reorder as cube
for i=1:P
    EM_pix{i}   = reshape(EM_pix{i}',nr,nc,L);
    IDX_comp{i} = reshape(IDX_comp{i}',nr,nc,L);
end







