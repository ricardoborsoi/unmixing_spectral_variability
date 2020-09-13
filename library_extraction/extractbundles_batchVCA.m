function [Lib]=extractbundles_batchVCA(Y, M0, bundle_nbr, percent)


P = size(M0,2);

[groups, bundle] = batchvca(Y, P, bundle_nbr, percent);
[groups,idx] = sort(groups);
bundle = bundle(:,idx);

Lib_temp = cell(P,1);
for i=1:P
    Lib_temp{i} = bundle(:,(i-1)*bundle_nbr+1 : i*bundle_nbr);
end

% order the signatures to M0
M_avg = zeros(size(M0));
for i=1:P
    M_avg(:,i) = mean(Lib_temp{i},2);
end


[~, idxPerm] = sort_endmembers_to_ref(M0,M_avg);

Lib = cell(P,1);
for i=1:P
    Lib{i} = Lib_temp{idxPerm(i)};
end



