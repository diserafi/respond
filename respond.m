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

function [u,iter,output] = respond(A,b,gamma,lambda,alpha0,alpha1,theta,a,u0,rho,tol,maxit,u_true)

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
% This is the main function of the ResPoND (Restoration of Poisson-Noisy
% Directional images) package. Given the noisy and blurry observed image
% 'b', the linear bluring operator 'A', the background noise 'gamma', and
% the regularization parameter 'lambda', this function restores 'b' by
% minimizing the KL-DTGV2 model described in [1], i.e., 'respond' solves
% the nonsmooth constrained optimization problem
% 
%  (P)             min  lambda * KL(A*u + gamma, b) + DTGV_2(u)
%                  s.t. x >= 0.
% 
% The objective function is the sum of a data-fidelity term, with weight
% lambda, consisting of the generalized Kullback-Leibler (KL) divergence of
% the blurred image A*u + gamma from the observed image b (see eq (2)
% in [1]) and a regularization term consisting of the discrete second-order
% Directional Total Generalized Variation (DTGV_2) (see eq (3) in [1]). The
% regularization term is the sum of two terms, accounting respectively for
% first- and second-order information on the image, whose relative weight
% is governed by the parameters alpha0 (first-order term) and alpha1
% (second-order term). For a detailed description of the discrete DTGV_2
% regularization term see Section 1 in [1] and the references therein.
% 
% The DTGV_2 regularization term depends on an estimate 'theta' of the
% angle (expressed in radiants) determining the direction of the main
% striped pattern, and on a scalar 'a', regulating the relative weight
% between the differences along the diretion determined by theta and the
% differences along the orthogonal direction. If an estimate of theta is
% not provided, the software automatically determines it by calling the
% 'dir_est_hough' function included in the package. Observe that, by
% setting (theta,a) = (0,1), the DTGV_2 becomes the standard TGV_2
% regularization. Hence, the present software can also be used for the
% restoration of general Poissonian images by minimizing the KL-TGV_2 model. 
% 
% The minimization problem (P) is solved by a specialized version of the
% Alternating Direction Method of Multipliers (ADMM) (Algorithm 2 in [1]).
% In detail, the algorithm is based on a reformulation of the original
% problem which allows one to rewrite it equivalently as
%
%               min  F_1(x) + F_2(z)
%               s.t. H*x + G*z = 0,
%
% in which F_1 and F_2 are closed, proper and convex functions, and the
% matrices H and G have full (column) rank. The reformulated problem is
% solved by an ADMM scheme where the subproblem in x coincides with the
% solution of a linear system that can be performed directly with a small
% computational cost by exploiting the Discrete Fourier Transform.
% The subproblem in z can be split in smaller subproblems, which can be
% solved in closed form.
%
% The ADMM algorithm stops when
%
%    ( ||u_k - u_{k-1}|| / ||u_{k-1}|| < tol ) or ( k = maxit ),
%
% where u_k is the approximation to the restored image at iteration k, and
% 'tol' and 'maxit' are a tolerance and a maximum number of iterations,
% respectively, given in input by the user (default values are considered
% if they are not provided).
%
% See Section 4 in [1] for further details.
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
% A      = double matrix, point spread function (PSF) representing the
%          blurring operator (periodic boundary conditions are assumed);
% b      = double matrix, containing the observed image, affected by blur
%          and Poisson noise; 
% gamma  = double, background noise; note that gamma is a scalar, i.e., we
%          assume the background noise to be the same for each pixel of b;
% lambda = double, weight of the Kullback-Leibler divergence;
% alpha0 = [optional] double, weight of the first order term in the DTGV/TGV
%          regularization term [default 0.5];
% alpha1 = [optional] double, weight of the second order term in the DTGV/TGV
%          regularization term [default 1-alpha0];
% theta  = [optional] double, angle (in radiants) determining the dierction
%          of the striped texture in the directional image; if left empty,
%          the parameter is estimated automatically via the 'dir_est_hough'
%          function;
% a      = [optional] double, weight associated with the direction theta+pi/2
%          in the DTGV regularization term [default 0.25];
% u0     = double matrix, starting guess; if empty then u0 is set to the
%          blurry and noisy image b;
% rho    = [optional] double, weight of the penalty term in the augmented
%          Lagrangian function, see eq. (12) in [1] [default 10]; 
% tol    = double, tolerance for the stopping criterion [default 1e-4];
% maxit  = integer, maximum number of iterations [default 200];
% u_true = [optional] double matrix, true image; if available then the
%          relative error w.r.t. the true image is stored at each iteration.
%
% OUTPUT ARGUMENTS
% 
% u      = double matrix, restored image;
% iter   = integer, number of ADMM iterations performed;
% output = struct containing further information on the execution, i.e.:
%          total_time = elapsed time;
%          conv_hist  = double vector, history of the relative
%                       distance (Frobenius norm) between iterates, i.e.,
%                           conv_hist(k) = ||u_k - u_{k-1}|| / ||u_{k-1}||,
%          rmse_vec   = [optional] double vector, if u_true is provided then
%                       it stores the RMSE at each iteration;
%          best_rec   = [optional] double matrix, if u_true is provided then
%                       it contains the "best reconstruction", i.e., the
%                       iterate associated with the minimum RMSE;
%          min_rmse_it = [optional] integer, if u_true is provided then
%                       it contains the iteration at which 'best_rec' has
%                       been obtained;
%          min_rmse   = [optional] double, if u_true is provided then it
%                       contains the RMSE at 'best_rec'.
%
%==========================================================================

if nargin<4
   error('You must specify at least the values of A, b, gamma, and lambda,\n for restoration with the TGV_2 regularization.'); 
end

%% Initializing parameters and checking for input
if ~exist('alpha0','var') || isempty(alpha0)
    alpha0 = 0.5;
end
if ~exist('alpha1','var') || isempty(alpha1)
    alpha1 = 1-alpha0;
end
if ~exist('theta','var') || isempty(theta)
    theta = dir_est_hough(b);  % automatic estimation of the angle
end
if ~exist('a','var') || isempty(a)
    a = 0.25;
end

if ~exist('u0','var') || isempty(u0)
    u0 = b;
end
if ~exist('rho','var') || isempty(rho)
    rho = 10;
end
if ~exist('tol','var') || isempty(tol)
    tol = 1e-4;
end
if ~exist('maxit','var') || isempty(maxit)
    maxit = 200;
end
if ~exist('u_true','var') || isempty(u_true)
    comp_rmse = false;
else
    comp_rmse = true; 
    sqrtN = sqrt(numel(u_true));
end

%% Starting time
tic

%% Definition of the derivative operators with periodic boundary conditions
[m,n] = size(u0);

c     = cos(theta);
s     = sin(theta); 
% Eigenvalues of the difference matrices
D1     = fft2([(c-s) -c; s 0],m,n);
conjD1 = conj(D1);
D2     = fft2(a*[(s+c) -s; -c 0],m,n);
conjD2 = conj(D2);
D1TD1  = abs(D1).^2;
D2TD2  = abs(D2).^2;

% Define difference operators via fft2
du1 = @(uin) real(ifft2(D1.*fft2(uin)));
du2 = @(uin) real(ifft2(D2.*fft2(uin)));

%% Initialization
% Unknowns
u    = u0;
w1   = du1(u);          % differences along main direction
w2   = du2(u);          % differences along secondary direction
z1   = zeros(m,n);      % constrained to be equal to A*u
z2_1 = zeros(m,n);      % constrained to be equal to D1*u-w1
z2_2 = zeros(m,n);      % constrained to be equal to D2*u-w2
z3_1 = zeros(m,n);      % constrained to be equal to (E*w)_1
z3_2 = zeros(m,n);      % constrained to be equal to (E*w)_2
z3_3 = zeros(m,n);      % constrained to be equal to (E*w)_3
z3_4 = zeros(m,n);      % constrained to be equal to (E*w)_4
z4   = zeros(m,n);      % constrained to be equal to u

% Scaled Lagrange multipliers
mu1   = zeros(m,n);     % size of z1
mu2_1 = zeros(m,n);     % size of z2_1
mu2_2 = zeros(m,n);     % size of z2_2
mu3_1 = zeros(m,n);     % size of z3_1
mu3_2 = zeros(m,n);     % size of z3_2
mu3_3 = zeros(m,n);     % size of z3_3
mu3_4 = zeros(m,n);     % size of z3_4
mu4   = zeros(m,n);     % size of z4

% Other variables
stop_crit  = 0;
iter       = 0;
if comp_rmse
    rmse_vec   = zeros(maxit,1);
    best_rec   = u0;
    min_rmse    = Inf;
    min_rmse_it = 0;
end
conv_hist = zeros(maxit,1);

%% Precomputing the matrix factorization for the update of u and w

conjA = conj(A);

% (1,1) block of the system matrix
Gamma = 1+abs(A).^2+D1TD1+D2TD2;

% Blocks of matrix Phi
Phi11 = 1+D1TD1+0.5*D2TD2;
Phi12 = 0.5*conjD2.*D1;
Phi21 = 0.5*conjD1.*D2;
Phi22 = 1+0.5*D1TD1+D2TD2;

% Computing Psi = inv(Phi) via the block-inversion formula
tmp1  = Phi11-Phi12.*Phi21./Phi22;
Psi11 = 1./(tmp1);
Psi12 = (-Phi12./Phi22)./(tmp1);
tmp2  = Phi22-Phi21.*Phi12./Phi11;
Psi21 = (-Phi21./Phi11)./(tmp2);
Psi22 = 1./(tmp2);

% Computing Xi^-1
invXi = Gamma-(conjD1.*Psi11.*D1+conjD2.*Psi21.*D1+conjD1.*Psi12.*D2+conjD2.*Psi22.*D2);
invXi = 1./invXi;

% Computing Omega
Omega11 = Phi11-D1./Gamma.*conjD1;
Omega12 = Phi12-D1./Gamma.*conjD2;
Omega21 = Phi21-D2./Gamma.*conjD1;
Omega22 = Phi22-D2./Gamma.*conjD2;

% Computing Y = inv(Omega) via the block-inversion formula
tmp1 = Omega11-Omega12.*Omega21./Omega22;
Y11  = 1./(tmp1);
Y12  = (-Omega12./Omega22)./(tmp1);
tmp2 = Omega22-Omega21.*Omega12./Omega11;
Y21  = (-Omega21./Omega11)./(tmp2);
Y22  = 1./(tmp2);

% Computing the blocks of Delta'*Psi
DtPsi1 = conjD1.*Psi11+conjD2.*Psi21;
DtPsi2 = conjD1.*Psi12+conjD2.*Psi22;

% Computing the blocks of Delta*inv(Gamma)
DinvG1 = D1./Gamma; 
DinvG2 = D2./Gamma;

%% Main loop starts
% To have a potentially good initialization of u and w, we avoid their
% first update. We move the update of u and w at the end of each iteration
% in the loop

while (~stop_crit) && (iter < maxit)

    %% Update of z = [z1,z2,z3,z4]
    % Build v_z
    vz1   = mu1   + real(ifft2(A.*fft2(u)));
    vz2_1 = mu2_1 + du1(u) - w1;
    vz2_2 = mu2_2 + du2(u) - w2;
    vz3_1 = mu3_1 + du1(w1);
    vz3_2 = mu3_2 + 0.5*(du2(w1)+du1(w2));
    vz3_3 = vz3_2;
    vz3_4 = mu3_4 + du2(w2);
    vz4   = mu4   + u;

    % Update of z1: exact 1-dimensional minimization
    rhotil = rho/lambda;
    B_eq = vz1 - 1/rhotil - gamma;
    C_eq = (rhotil*gamma*vz1 + b - gamma)/rhotil;
    z1 = 0.5*(B_eq + sqrt(B_eq.^2 + 4*C_eq));
    
    % Update of z2: shrinkage operator
    alpha0til = alpha0/rho;

    vz2n = sqrt(vz2_1.^2 + vz2_2.^2);
    vz2n(vz2n==0) = 1;
    vz2n = max(vz2n - alpha0til, 0)./vz2n;
    z2_1 = vz2_1.*vz2n;    
    z2_2 = vz2_2.*vz2n;

    % Update of z3: shrinkage operator
    alpha1til = alpha1/rho;

    vz3n = sqrt(vz3_1.^2 + vz3_2.^2 + vz3_3.^2 + vz3_4.^2);
    vz3n(vz3n==0) = 1;
    vz3n = max(vz3n - alpha1til, 0)./vz3n;
    z3_1 = vz3_1.*vz3n;
    z3_2 = vz3_2.*vz3n;
    z3_3 = vz3_3.*vz3n;
    z3_4 = vz3_4.*vz3n;
    
    % Update of z4: projection
    z4 = max(0, vz4);
    
    
    %% Update of the Lagrange multipliers mu
    mu1   = vz1   - z1;
    mu2_1 = vz2_1 - z2_1;
    mu2_2 = vz2_2 - z2_2;
    mu3_1 = vz3_1 - z3_1;
    mu3_2 = vz3_2 - z3_2;
    mu3_3 = mu3_2;
    mu3_4 = vz3_4 - z3_4;
    mu4   = vz4   - z4;
        
    iter = iter+1;
    uold = u;
    
    %% Update of x=[u,w]: exact solution of linear system via FFT
    
    % Build v_x
    vx1   = z1   - mu1;
    vx2_1 = z2_1 - mu2_1;
    vx2_2 = z2_2 - mu2_2;
    vx3_1 = z3_1 - mu3_1;
    vx3_2 = z3_2 - mu3_2;
    vx3_3 = z3_3 - mu3_3;
    vx3_4 = z3_4 - mu3_4;
    vx4   = z4   - mu4;
    
    % Compute the rhs of the linear system
    rhs_1 = conjA.*fft2(vx1) + conjD1.*fft2(vx2_1) + conjD2.*fft2(vx2_2) + fft2(vx4);
    rhs_2 = - fft2(vx2_1) + conjD1.*fft2(vx3_1) + 0.5*conjD2.*fft2(vx3_2+vx3_3);
    rhs_3 = - fft2(vx2_2) + conjD2.*fft2(vx3_4) + 0.5*conjD1.*fft2(vx3_2+vx3_3);

    % Compute the solution to the linear system by matrix-vector products
    y1 = rhs_1+DtPsi1.*rhs_2+DtPsi2.*rhs_3;
    y2 = DinvG1.*rhs_1+rhs_2;
    y3 = DinvG2.*rhs_1+rhs_3;

    u  = invXi.*y1;        u  = real(ifft2(u));
    w1 = Y11.*y2+Y12.*y3;  w1 = real(ifft2(w1));
    w2 = Y21.*y2+Y22.*y3;  w2 = real(ifft2(w2));

    %% Compute and store the error
    rel_diff = norm(u-uold,'fro')/norm(u,'fro');
    conv_hist(iter) = rel_diff;
    if comp_rmse
        rmse = norm(u-u_true,'fro')/sqrtN;
        if rmse < min_rmse
            min_rmse = rmse;
            best_rec = u;
            min_rmse_it = iter;
        end
        rmse_vec(iter) = rmse;
    end

    % Check if the stopping criterion is satisfied
    stop_crit = (rel_diff<tol);
end

output.total_time = toc;
output.conv_hist = conv_hist(1:iter);
if comp_rmse
    output.rmse_vec = rmse_vec(1:iter);
    output.min_rmse = min_rmse;
    output.min_rmse_it = min_rmse_it;
    output.best_rec = best_rec;
end

end
