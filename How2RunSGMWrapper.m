function [] = How2RunSGMWrapper()
% Compute disparity map with semi global matching algorithm.
% "Heiko Hirschmuller. Stereo processing by semiglobal matching and mutual information. 
%   Pattern Analysis and Machine Intelligence, IEEE Transactions on, 30(2):328–341, 2008."


%% Initialize parameters

% Block size for raw cost aggergation. simple block matching
params.block_size = 5; 

% [minimum disparity, maximum disparity
params.disparity_range = [-16 16]; 

%Penlize disparity diffrent than 1 for neigber 
params.P1 = (1/16)*8*params.block_size.^2; 

% Penlize for more than 1 disparity diffrent from the neighber
params.P2 = (1/16)*32*params.block_size.^2; 

% Number of SGM path direction
params.directions_num = 8; 

%% Load images

% left image in gray scale. used as the reference image
im_left  = rgb2gray(imread('image_left.png'));

% right image in gray scale. used as the support image
im_right = rgb2gray(imread('image_right.png'));



%% Run Semi Global Matching algorithm

% the dispairty map from left view point with values between 
% params.disparity_range(1) %and params.disparity_range(2)
disparity_map = SGMWrapper( im_left, im_right, params );


end

