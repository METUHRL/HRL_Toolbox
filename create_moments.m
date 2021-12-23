function [Cx, mnx] = create_moments(qrs)

% CREATE_MOMENTS    calculate the mean vector and the covariance matrix of a given data
% 
% Usage:
%   [Cx, mnx] = create_moments(qrs)
%
% Inputs:
%   qrs     (nLeads x nFrames) data for which we calculate the mean and the
%           covariance. If it is not provided, a dialog box appears and the
%           users can make their selections from previously saved .mat
%           file. But the variable name saved inside the data file should
%           be qrs
%
% Outputs:
%   Cx      Covariance matrix
%   mnx     Mean vector (mean is along the time dimension)
%
% Author: Assoc. Prof. Yesim Serinagaoglu Dogrusoz <yserin@metu.edu.tr>

if nargin < 1
    [TSFileName, rpathname] = uigetfile('*.mat', 'Load training dataset');
    load([rpathname TSFileName]); % data should be saved with the variable name 'qrs'
end;
fprintf('\n Training dataset loaded \n');

nFrames = size(qrs,2);

mnx = mean(qrs,2);
mnxmat = repmat(mnx, 1, nFrames);
Cx = (qrs - mnxmat) * (qrs - mnxmat)' / nFrames; % Maybe we should use (nFrames-1) 
