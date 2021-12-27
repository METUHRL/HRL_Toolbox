function [X,lambda] = tikhonovRT(Y, A)
% TIKHONOVRT    Tikhonov Regularization, uses Regularization toolbox by
% Per Christian Hansen
%
% Inputs:
%       Y: Torso measurements
%       A: Forward transfer matrix
%
% Outputs:
%       X: the inverse problem solution 
%       lambda: regularization parameter lambda vector
%
% Usage:
%       [X,lambda] = tikhonovRT(Y, A)
%
% Author: 
%       Assoc. Prof. Yesim Serinagaoglu Dogrusoz <yserin@metu.edu.tr>


nFrames = size(Y,2);
[U,s,V] = csvd(A);

lambda = zeros(nFrames, 1);

X = zeros(size(A,2), nFrames);

for fr = 1:nFrames,
   [lambda(fr),wi,wi1,wi2] = l_curve(U, s, Y(:,fr));
   if lambda(fr)>2
       if fr == 1,
           lambda(fr) = 0.05;
       else
           lambda(fr)=lambda(fr-1);
       end
   end
   if lambda(fr)<0.0005
       if fr == 1,
           lambda(fr) = 0.05;
       else
           lambda(fr)=lambda(fr-1);
       end
   end
   X(:,fr) = tikhonov(U,s,V,Y(:,fr),lambda(fr));

end;
