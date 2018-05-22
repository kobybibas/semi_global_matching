function [ disparity_map ] = SGMWrapper( im_left, im_right, params )
% Compute disparity map with semi global matching algorithm.
% "Heiko Hirschmuller. Stereo processing by semiglobal matching and mutual information. 
%   Pattern Analysis and Machine Intelligence, IEEE Transactions on, 30(2):328–341, 2008."
%
% Inputs:   im_left - left image in gray scale. used as the reference image
%           im_right - right image in gray scale. used as the support image
%           params - diffrent parameters for thr algorithm: 
%               block_size, disparity_range, P1, P2, directions_num
% Output:   disparity_map - the disparity map from left view point
%               with values between params.disparity_range(1) 
%               and params.disparity_range(2)


%% Check inputs parameters. put default param if not exist

% Block size for raw cost aggregation. simple block matching
if ~isfield(params, 'block_size')
    params.block_size = 5; 
end

% [minimum disparity, maximum disparity] - disparity range
if ~isfield(params, 'disparity_range')
    params.disparity_range = [-16 16]; 
end

% Penalise for 1 disparity different from the neighbor
if ~isfield(params, 'P1')
    params.P1 = (1/16)*8*params.block_size.^2; 
end

% Penalise for more than 1 disparity different from the neighbor
if ~isfield(params, 'P2')
    params.P2 = (1/16)*32*params.block_size.^2; 
end

% Number of SGM path direction. currently only 8 is supported.
if ~isfield(params, 'directions_num')
    params.directions_num = 8; 
end

%% Raw cost aggregation 
% Compute cost using simple block matching with Sum-Square-Diff distance  
block_size      =   params.block_size;
disparity_range =   params.disparity_range;
C = costCalculation(im_left, im_right, block_size, disparity_range);

%% Cost aggregation (Dynamic programming)
% calculates cost for each pixel and for each disparity along some path.
% Penalise disparity differences between neighbors.
P1              =   params.P1;
P2              =   params.P2;
directions_num  =   params.directions_num;
Lr = costAggregation(C, P1, P2, directions_num);

%% Choose disparity from all paths
% for each pixel for each disparity summing all costs of the different path.
% choosing disparity with minimum cost.
disparity_range = params.disparity_range;
disparity_map = disparitySelection(Lr, disparity_range);

end




function disparity_map = disparitySelection(Lr , disparity_range)
% Generate disparity map from aggregated path.
% Current implementation is summing all paths and choose the minimum cost. 
assert(disparity_range(2) >  disparity_range(1));

% Sum all paths. eq. (14) in the paper
S = sum(Lr, 4);
[~, disparity_map] = min(S, [], 3);

%Normalize based on disparity range
disparity_map = disparity_map + disparity_range(1) -1 ; %-1 because Matlab idx start from 1.
end

function Lr = costAggregation(C, P1, P2, directions_num)
% Inputs:   C - cost array, dim: ImageHeight x ImageWidth x DisparityRange
%           P1 - Penalise for 1 disparity different from the neighbor
%           P2 - Penalise for more than 1 disparity different from the neighbor
%           directions_num - number of direction for which to calculate the patch cost
% Outputs:  Lr - for each pixel and for each disparity the calculated minimal costs along the paths. 
%               dim: ImageHeight x ImageWidth x DisparityRange


%currently only 8 directions are implemented
assert(directions_num == 8); 
assert(P1 > 0); 
assert(P2 > 0);


[size_y ,size_x , size_d] = size(C);
max_b = size_x + size_y;

%Initialize output
Lr = zeros([size(C),directions_num], 'double');

%Iterate over direction
for dir = 1:directions_num
    Li = zeros(size(C), 'double'); %temp buffer
    
    %Iterate over all paths in specific direction
    for b = -1*max_b : max_b
        
        %Extract path from Cost array
        inds = directionBasedRemap(dir, b, size_x, size_y, size_d);
        slice = reshape(C(inds), [size(inds,1)/size_d, size_d]); %%dimensions
        
        %If path exist evaluate the path cost
        if all(size(slice) ~= 0)
            
            %Evaluate cost over the path
            graded_slice = evaluatePathCost(slice',P1,P2)';
            Li(inds) = graded_slice(:);
            
            %Assign to output buffer
            Lr(:,:,:,dir) = Li;
        end        
    end
end

end

function [inds] = directionBasedRemap(direction, b, size_x, size_y, size_d)
% Parameterization of linear line a1*y = a2*x + b 
% Each direction means different parameters for a1, a2 and b
% Inputs:   direction - direction of the path. value between 1 and 8
%           b - parametrization of the line.
%           [size_x, size_y, size_d] - dimensions of C array (Cost array)
% Outputs: inds - indices of array C which represent the path

%--------------------------------------
% direction | a1 |  a2 |  flip_order
% -------------------------------------
%   1       | 1  |  0  |      0       |
%   2       | 1  |  1  |      0       |
%   3       | 0  |  1  |      0       |
%   4       | 1  |  -1 |      1       |
%   5       | 1  |  0  |      1       |
%   6       | 1  |  1  |      1       |
%   7       | 0  |  1  |      1       |
%   8       | 1  |  -1 |      0       |

% Note: probably there is a better way to do it. need to be improved 


% direction mapping
switch direction
    case 1
        a1=1; a2=0; flip_order=0;
    case 2
        a1=1; a2=1; flip_order=0;
    case 3
        a1=0; a2=1; flip_order=0;
    case 4
        a1=1; a2=-1; flip_order=1;
    case 5
        a1=1; a2=0; flip_order=1;
    case 6
        a1=1; a2=1; flip_order=1;
    case 7
        a1=0; a2=1; flip_order=1;
    case 8
        a1=1; a2=-1; flip_order=0; 
end

if (a1 ~= 0)
    x_inds = 1:size_x;
    y_inds = (a2*x_inds+b)*a1; % a1 is either 1 or -1
    inds_in = find(y_inds >= 1 & y_inds <= size_y);
else
    y_inds = 1:size_y;
    x_inds = -1*b*a2*ones(size(y_inds));
    inds_in = find(x_inds >= 1 & x_inds <= size_x);
end

x_inds = x_inds(inds_in);
y_inds = y_inds(inds_in);


slice_length = size(x_inds, 2);
x_inds = repmat(x_inds, 1, size_d)';
y_inds = repmat(y_inds, 1, size_d)';
z_inds = repmat(1:size_d, slice_length, 1);
z_inds = z_inds(:);

% Flip based on direction
if flip_order
    x_inds = fliplr(x_inds')';
    y_inds = fliplr(y_inds')';
end

% Merge indices to inds vector- the remap vector
inds = sub2ind([size_y size_x size_d], y_inds, x_inds, z_inds);
end % directionBasedRemap



function Lr_slice = evaluatePathCost(C_slice, P1, P2)
% Evaluate the path cost: for each position and for each disparity
% calculate the cost using dynamic programming.
% Inputs:   C_slice - path generated from C array. from which path cost will be calculated.
%           P1 - Penalise for 1 disparity different from the neighbor
%           P2 - Penalise for more than 1 disparity diffrent from the neighbor
% Outputs: Lr_slice - The cost the path C_slice. same dim as C_slice


nLabels = size(C_slice, 1);
nCols = size(C_slice, 2);

% create const_grades:
% the constant part of the cost for disparity from YY to move to disp from XX
xx = 1:nLabels;
[YY, XX] = meshgrid(xx,xx);
const_grades = zeros(nLabels, 'double');
const_grades(abs(XX - YY) == 1) = P1;
const_grades(abs(XX - YY) > 1) = P2;

% init L_slice
Lr_slice = zeros(size(C_slice));
Lr_slice(:, 1) = C_slice(:,1);

% iterate over slice direction
for col = 2:nCols
    L_slice_prev = Lr_slice(:,col-1);
    
    % calculate C and M. Eq. (12) in the paper
    C = C_slice(:,col);
    M = min(repmat(L_slice_prev',nLabels,1) + const_grades, [], 2);
    
    % save the result in L_slice. Eq. (13) in the paper
    Lr_slice(:,col) = C + M - min(L_slice_prev);

end

end % evaluatePathCost


function [C] = costCalculation(im_reference, im_support, block_size, disparity_range)
% Calculate cost array using simple block matching algorithm
% Inputs:   im_left - left image, gray scale, will be converted to double.
%               size: WxHx1
%           im_right - right image, gray scale, will be converted to double
%               size: WxHx1
%           block_size - window size of the matching process
%           disparity_range - vector of the kind [minimum disparity maximum disparity]
% Output:   C - cost calculation result. size WxHxD where D is disparity range


%Check inputs arguments
assert(all(size(im_reference) == size(im_support)));
assert(block_size > 1); 
assert(disparity_range(2) > disparity_range(1)); 


if (~isa(im_reference,'double') || ~isa(im_support,'double'))
    im_reference = im2double(im_reference);
    im_support  = im2double(im_support);
end

% output image:
C = zeros(size(im_support,1), size(im_reference,2), disparity_range(2) - disparity_range(1) + 1, 'double');


% allocate im_diff. im_diff is in left coords:
im_diff = zeros(size(im_reference, 1), size(im_reference,2), 'double');

% iterate over all disparities
for disp = disparity_range(1) : disparity_range(2)
    

    % calc borders for images    
    borders_left = [1, size(im_reference,2)] - disp *  [(disp < 0), disp > 0];
    borders_right = borders_left + [disp, disp];
    
    % im_diff = (im_left - im_reference).^2, inside borders
    im_diff(:, borders_left(1) : borders_left(2)) = ...
        im_support(:, borders_left(1) : borders_left(2)) - ...
        im_reference(:, borders_right(1) : borders_right(2));
    im_diff = abs(im_diff);
    
    % replicate areas outside borders for im_diff:
    im_diff(:, 1:borders_left(1)) = im_diff(:, borders_left(1)) * ...
        ones(1, borders_left(1));
    im_diff(:, borders_left(2):end) = im_diff(:,borders_left(2)) * ...
        ones(1, size(im_diff,2) - borders_left(2) + 1);
    
    % for summing, do convolution with ones:
    filt = ones(block_size);
    C(:, :, disp - disparity_range(1) + 1) = ...
        conv2(im_diff, filt, 'same');
    
end

% copy last supported value to boundaries of the output:
frame_copy = (block_size - 1) / 2;

% replicate areas outside borders for im_diff:
% replicate on x
C(:, 1:frame_copy,:) = repmat(C(:, frame_copy+1,:),1,frame_copy);
C(:,end-frame_copy+1:end,:) = repmat(C(:,end-frame_copy, :), 1, frame_copy);
% replicate on y
C(1:frame_copy, :, :) = repmat(C(frame_copy+1,:,:),frame_copy,1);
C(end-frame_copy+1:end,:,:) = repmat(C(end-frame_copy,:,:),frame_copy,1);

% normalize output:
min_val = min(C(:));
max_val = max(C(:));

C = 255.0*(C - min_val)/(max_val-min_val);

end