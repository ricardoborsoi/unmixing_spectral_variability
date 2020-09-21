function [mtxout] = row2col_lexico_order(inmtx, nr, nc)
% -------------------------------------------------------------------------
% convert row-wise to column-wise lexicographic ordering
% 
% Author: Ricardo Borsoi, 2018
% -------------------------------------------------------------------------



N = nr*nc;


if length(size(inmtx)) == 2
    K = size(inmtx,1);
    
    % Convert image matrix back to cube
    im_cube = permute(reshape(inmtx',nc,nr,K),[2 1 3]); 
    
    % Convert cube to matrix with column lexicographic ordering 
    mtxout = reshape(im_cube,N,K)'; 
    
    
    
elseif length(size(inmtx)) == 3
    L = size(inmtx,1);
    P = size(inmtx,2);
    mtxtmp = zeros(L,P,nr,nc);
    
    % Convert 3D endmember matrix to 4D EM matrix indexed spatially
    n = 0;
    for i=1:nr
        for j=1:nc
            n = n+1;
            mtxtmp(:,:,i,j) = inmtx(:,:,n);
        end
    end
    
    % Convert 3D 4D EM matrix to 3D matrix with column lexicographic ordering
    mtxout = zeros(L,P,N);
    n = 0;
    for j=1:nc
        for i=1:nr
            n = n+1;
            mtxout(:,:,n) = mtxtmp(:,:,i,j);
        end
    end
    
    
else
    error('Unknown input!')
end











