function [] = How2RunSGMWrapper()
% Compute disparity map with semi global matching algorithm.
% "Heiko Hirschmuller. Stereo processing by semiglobal matching and mutual information. 
%   Pattern Analysis and Machine Intelligence, IEEE Transactions on, 30(2):328�341, 2008."


%% Initialize parameters

% Block size for raw cost aggergation. simple block matching
params.block_size = 5; 

% [minimum disparity, maximum disparity
params.disparity_range = [-16 16]; 

% Penalise disparity different than 1 for neighbor 
params.P1 = (1/16)*8*params.block_size.^2; 

% Penalise for more than 1 disparity different from the neighbor
params.P2 = (1/16)*32*params.block_size.^2; 

% Number of SGM path direction
params.directions_num = 8; 

%% Load images

% left image in gray scale. used as the reference image
im_left  = rgb2gray(imread('image_left.png'));

% right image in gray scale. used as the support image
im_right = rgb2gray(imread('image_right.png'));



%% Run Semi Global Matching algorithm
% the disparity map from left view point with values between 
% params.disparity_range(1) %and params.disparity_range(2)
disparity_map = SGMWrapper( im_left, im_right, params );
imwrite(uint8(255*(disparity_map-min(disparity_map(:)))/(max(disparity_map(:)) - min(disparity_map(:)) )),...
    'disparity_map.tiff');

end

