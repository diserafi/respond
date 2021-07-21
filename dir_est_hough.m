% -------------------------------------------------------------------------
% Copyright (C) 2021 by D. di Serafino, G. Landi, M. Viola.
%
%                           COPYRIGHT NOTIFICATION
%
% This file is part of ResPoND.
%
% ResPoND is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ResPoND is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ResPoND. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

function [theta] = dir_est_hough(b)

%==========================================================================
%
% Authors:
%   Daniela di Serafino (daniela.diserafino [at] unina.it)
%   Germana Landi       (germana.landi [at] unibo.it )
%   Marco Viola         (marco.viola [at] unicampania.it)
%
% Version: 1.1
% Last Update: 18 July 2021
%
%==========================================================================
%
% This function estimates the direction of the main striped pattern in a
% directional image. After extracting the edges by applying the Sobel
% filter, the direction is estimated by using the Hough transform to
% detect lines in the edge image.
%
% See Section 3 in [1] for further details.
%
%==========================================================================
%
% REFERENCES:
% [1] D. di Serafino, G. Landi and M. Viola,
%     "Directional TGV-Based Image Restoration under Poisson Noise",
%     Journal of Imaging, 7 (6), p. 99, 2021, DOI: 10.3390/jimaging7060099,
%     open access.
%
%==========================================================================
%
% INPUT ARGUMENTS
% 
% b      = double matrix, containing the (possibly blurry and noisy)
%          observed image.  
%
% OUTPUT ARGUMENTS
% 
% theta  = double, estimate of the angle (in radiants) determining the
%          striped pattern direction. 
%
%==========================================================================

[m,n] = size(b);

% Detect the edges via the Sobel filter
edgeb = edge(b);

% Create the largest possible circular mask centered in the image center
h = min(m,n);
mask = fspecial('disk',floor(0.5*(h-1)));
mask = padarray(mask, floor((size(b)-size(mask))/2));
if size(mask,1)<m || size(mask,2)<n
    mask(m,n) = 0;
end

% Apply mask to the edge image
edgeb = edgeb.*sign(mask);

% Compute the Hough transform 
[H,eta_vect] = hough(edgeb);

% Compute scores for each eta and determine the eta with maximum score
sumH2 = sum(H.^2,1); % max(H,[],1); % 
[~,imax] = max(sumH2);
etamax = eta_vect(imax);

% Determine the estimate theta
if etamax >= 0
    theta = 90-etamax;
else
    theta = -90-etamax;
end
theta = theta/180*pi;

end