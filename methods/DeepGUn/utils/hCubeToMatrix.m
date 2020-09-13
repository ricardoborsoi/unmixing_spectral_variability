function [ hypMatrix ] = hCubeToMatrix(hCube,invertRowCol)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    invertRowCol = 0;
end

[n1,n2,L] = size(hCube);


hypMatrix = zeros(L,n1*n2);

count = 1;

if ~invertRowCol
    for i=1:n1,
        for j=1:n2
            hypMatrix(:,count) = squeeze(hCube(i,j,:));
            count = count + 1;
        end
    end
else
    for j=1:n2
        for i=1:n1
            hypMatrix(:,count) = squeeze(hCube(i,j,:));
            count = count + 1;
        end
    end
end

