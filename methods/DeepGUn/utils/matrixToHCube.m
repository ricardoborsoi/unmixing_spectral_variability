function [ HCube ] = matrixToHCube( Y, nRow, nCol,invRC)
%function [ HCube ] = matrixToHCube( Y, nRow, nCol,invRC)
%   
%   invRC - inverts column-wise (default) to row-wise conversion.

if nargin < 4
    invRC = 0;
end

[L,N] = size(Y);

if (N ~= nRow*nCol)
    error('Size missmatch! Check number of columns and rows.')
end




HCube = zeros(nRow, nCol, L);
count =1;
if invRC==0
    for i=1:nRow, 
        for j=1:nCol,
            HCube(i,j,:) = Y(:,count);
            count = count +1;
        end
    end
else
    for j=1:nCol,
        for i=1:nRow, 
            HCube(i,j,:) = Y(:,count);
            count = count +1;
        end
    end    
end

end

