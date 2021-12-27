function [X2] = tikhonovRT_singLam(Y, A, singLambda)

% Uses a single (median) lambda value at all times
% TIKHONOVRT    Tikhonov Regularization, uses Regularization toolbox by 
% Per Christian Hansen

%
% Inputs:
%       Y: Torso measurements
%       A: Forward transfer matrix
%       sigLambda: regularization parameter
%
% Output:
%       X2: the inverse problem solution
%
% Usage:
%       [X2] = tikhonovRT_singLam(Y, A, singLambda)
%
% Author:
%       Assoc. Prof. Yesim Serinagaoglu Dogrusoz <yserin@metu.edu.tr>



nFrames = size(Y,2);
[U,s,V] = csvd(A);
%[U2,sm,XX,V2] = cgsvd(A,L);

X2 = zeros(size(A,2), nFrames);

for fr = 1:nFrames,
    X2(:,fr) = tikhonov(U,s,V,Y(:,fr),singLambda);
end;
