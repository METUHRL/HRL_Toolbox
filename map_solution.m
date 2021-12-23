function [Xinv, Cerr] = map_solution(Y, A, meanX, Cx, Cn)

% MAP_SOLUTION  
% Find the inverse solution using the Bayesian MAP approach when the data
% and the noise are Gaussian distribution. 
%
% Usage:
%   [Xinv, Cerr] = map_solution(Y, A, meanX, Cx, Cn)
%
% Inputs:
%   Y       Torso measurements
%   A       Forward transfer matrix
%   meanX   (nHeart x 1) mean vector of epicardial potentials
%   Cx      Covariance matrix of the epicardial potentials
%   Cn      Covariance matrix of the noise
%
% Outputs:
%   Xinv    Solution (MAP estimate of epicardial potentials)
%   Cerr    Error covariance matrix (YSD thesis, p.60, Eq. (2.35))
%
%
% When the prior is a Gaussian pdf, MAP estimation is the same as the MMSE
% estimat5ion. Therefore, this code uses the MMSE approach to the solution,
% i.e., the solution is (YSD thesis, p.54, Eq. (2.20)): 
%
%       Xinv = meanX + (Cx*A') * inv(A*Cx*A'+Cn) * (y - A*meanX)
%
% Author: Assoc. Prof. Yesim Serinagaoglu Dogrusoz <yserin@metu.edu.tr>

%[nTorso, nHeart] = size(A);
nFrames = size(Y,2);

CxAtr = Cx*A';
K = A*CxAtr + Cn;
B = CxAtr/K;

meanX = repmat(meanX, 1, nFrames); 
Xinv = meanX + B*(Y - A*meanX);

% the error covariance matrix:
Cerr = Cx - B*A*Cx;
