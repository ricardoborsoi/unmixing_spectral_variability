function networks = loadNetworks(baseFilename,P)
% -------------------------------------------------------------------------
% Load keras network weights
% 
% Author: Ricardo Borsoi
% last revision: 02/09/2019
% -------------------------------------------------------------------------

% this command get the information in the file
% hinfo = hdf5info('vae_EM_idx1.h5');
% 
%                       enc/dec    net_layer     bias/weights
% hinfo.GroupHierarchy. Groups(1) .Groups(4).    Datasets(1)

actFunStr = 'ReLU';


networks = cell(1,P);
for em=1:P
    % hinfo = hdf5info(['vae_EM_idx' num2str(em) '.h5']);
    hinfo = hdf5info([baseFilename num2str(em) '.h5']);
    numLayers = length(hinfo.GroupHierarchy.Groups(1).Groups);
    
    networks{em} = cell(1,numLayers);
    for layer=1:numLayers
        networks{em}{layer}.b = hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(layer).Datasets(1)); % bias
        networks{em}{layer}.W = hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(layer).Datasets(2)); % weights
        if layer < numLayers
            switch actFunStr
                case {'relu','ReLU'}
                    networks{em}{layer}.actFun = @(x)ReLU(x);
                case {'softplus'}
                    networks{em}{layer}.actFun = @(x)softplus(x);
                otherwise
                    error('Unkwnon activation function selected!')
            end
        else
            % Use a sigmoid at the output layer
            networks{em}{layer}.actFun = @(x)sigmoid(x);
        end
    end
end

end


% activation fucntions ----------------------------------------------------
function out = ReLU(in)
    out = max(in,0);
end

function out = softplus(in)
    out = log(exp(in) + 1);
end

function out = sigmoid(in)
    out = 1./ (1 + exp(-in));
end


