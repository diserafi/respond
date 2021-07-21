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

function plot_line_rad(theta,N,co,ls)

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
% Given an angle theta, this function plots a line over the current matlab
% figure, starting from the image center and using the provided values for
% color and line style.
%
%==========================================================================
%
% INPUT ARGUMENTS
% 
% theta = double, angle in radiants;
% N     = integer, smallest image size;
% co    = line color, specified as a MATLAB ColorSpec (Color Specification),
%         e.g., for green it can have the values 'green', 'g', or [0 1 0];
% ls    = char, linestyle specified as one of the MATLAB linestyles, e.g.,
%         '--' for dashed line or ':' for dotted line.
%
%==========================================================================
% 
% This function is a slight modification of the function 'plot_line' by
% Rasmus Dalgas Kongskov, included in the 'DTGV-Reg' package, available
% from http://www2.compute.dtu.dk/~pcha/HDtomo.
%
%==========================================================================

if nargin<4 || isempty(ls)
    ls = '-';
end

if theta>(pi/2)
    theta = pi - theta;
else
    theta = -theta;
end

hold on;
line([N/2 N/2*(1 + cos(theta)) ],[N/2 N/2*(1 + sin(theta))],...
                'LineWidth',4,'LineStyle',ls,'Color',co);
hold off

end 