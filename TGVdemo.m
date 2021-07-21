%==========================================================================
%        EXAMPLE OF USE OF RESPOND FOR TGV-BASED IMAGE RESTORATION
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
% This script illustrates the use of the of the ResPoND package for
% TGV-based restoration of Poisson-noise corrupted images.
%
% The script loads a phantom image and creates a test problem with either
% Gaussian or Out-Of-Focus blur, adds Poisson noise, and then calls the
% 'respond' function for solving the image restoration problem. Various
% error metrics are computed to assess the quality of the restoration.
%
%==========================================================================

clear all
close all

%% Choose the blur type
blur_ind = 1; % 1 - Out-Of-Focus blur with radius r = 5
              % 2 - Gaussian blur with standard deviation sigma = 2

fprintf(1,'\nStarting TGVdemo\n\n');

%% Create the test problem
% Load phantom
load smooth_phantom
u_true = u_true/max(u_true(:)); 

% Show original image
figure; imagesc(u_true); colormap gray; axis equal; axis off; title('Original','FontSize',15);

% Build blurring operator
[m,n] = size(u_true);
switch blur_ind
    case 1
        % Out-Of-Focus blur - radius 5
        h = fspecial('disk', 5);
        psf = zeros(m,n); center = [ceil(m/2), ceil(n/2)];
        psf(center(1)-5:center(1)+5,center(2)-5:center(2)+5)=h;
        psf = psf/sum(psf(:));
    case 2
        % Gaussian blur - std 2
        psf = fspecial('gaussian',[m,n],2);
end
S =  fft2(fftshift(psf));

% Blur image and add background
b = real(ifft2(S.*fft2(u_true)));
bg = 1e-10; 
b = b+bg;

% Add Poisson noise
rng('default')
mag   = 7.5e7; % noise level
scale = mag/sum(sum(b))*1e-12;
b     = imnoise(b*scale,'poisson')/scale;
snr   = 10*log10(mag/sqrt(mag));

% Show blurry and noisy image
figure; imagesc(b); colormap gray; axis equal; axis off;
titlestring = sprintf('Blurry and noisy (SNR = %.2f)',snr);
title(titlestring,'FontSize',15);

% Define regularization parameters
switch blur_ind
    case 1
        lambda = 15;
    case 2
        lambda = 16; 
end

% TGV parameters 
alpha0 = .1;
alpha1 = (1-alpha0);

%% Set ADMM parameters
% Starting point
u0 = b;
% Stopping conditions
maxit = 200; 
tol = 1e-4;
% Penalty parameter
rho = 10;

%% Call 'respond' to solve the restoration problem
[u,it,output] = respond(S,b,bg,lambda,alpha0,alpha1,0,1,u0,rho,tol,maxit,u_true);

%% Results
% Show restored image
figure; imagesc(u); colormap(gray); axis equal; axis off; title('Restored','FontSize',15);

% Plot RMSE history
figure; plot(output.rmse_vec,'LineWidth',2); title('RMSE history','FontSize',15);
xlabel('iter','FontSize',15);
ylabel('RMSE','FontSize',15);

% Compute ISNR
ISNR = 10*log10((norm(u_true-b,'fro').^2)/(norm(u_true-u,'fro').^2));

% Compute MSSIM
K = [0.01 0.03];
window = fspecial('gaussian', 11, 1.5);
MSSIM = ssim_index(u_true, u, K, window, max(max(u_true)));

% Display metrics
fprintf('RMSE        ISNR        MSSIM       It   Time\n');
fprintf('%6.4e  %6.4e  %6.4e  %d  %4.2f \n',...
    output.rmse_vec(end),ISNR,MSSIM,it,output.total_time);
