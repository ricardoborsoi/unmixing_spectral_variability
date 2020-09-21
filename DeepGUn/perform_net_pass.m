function [outp]=perform_net_pass(network,input)
% compute predictions from the network based on a latent vector
% 
% Author: Ricardo Borsoi
% last revision: 02/09/2019

netLen = length(network);

% Use activation functions stored in the layers ------------------------
for layer = 1:netLen
    temp = network{layer}.W * input + network{layer}.b;
    input = network{layer}.actFun(temp);
end
outp = input;

end


% activation functions ----------------------------------------------------
function out = ReLU(in)
    out = max(in,0);
end

function out = softplus(in)
    out = log(exp(in) + 1);
end



