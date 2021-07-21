# ResPoND

## A MATLAB package for the restoration of directional images corrupted by Poisson noise.

## Authors
Daniela di Serafino, University of Naples Federico II, Napoli, Italy, daniela.diserafino[at]unina.it
Germana Landi, University of Bologna, Bologna, Italy, germana.landi[at]unibo.it
Marco Viola, University of Campania "Luigi Vanvitelli", Caserta, Italy, marco.viola[at]unicampania.it

## Last Update
Version 1.1 - July 18, 2021

## Description
ResPoND (Restoration of Poisson-Noisy Directional images) is a MATLAB 
package for deblurring directional images corrupted by Poisson noise.

The main function in the package, named `respond`, performs the restoration
by minimizing the KL-DTGV_2 model presented in [1]. Given the noisy and
blurry observed image `b`, the linear bluring operator `A`, the background
noise `gamma`, and the regularization parameter `lambda`, the function
finds an approximate solution to the nonsmooth constrained minimization
problem

              min  lambda*KL(A*u + gamma, b) + DTGV_2(u)
              s.t. x >= 0,

where KL is the so-called Kullback-Leibler divergence of the blurred image
`A*u + gamma` from the observed image `b` and DTGV_2 is the discrete
second-order Total Generalized Variation regularization term. The
minimization problem is solved by a specialized version of a two-block
Alternating Direction Method of Multipliers (ADMM) (see Algorithm 2 in [1]).

The `respond` function is also suited for the deblurring of general
Poissonian images by the minimization of the KL-TGV_2 model (see the
function documentation for further details).

### References
[1] D. di Serafino, G. Landi and M. Viola,
*Directional TGV-Based Image Restoration under Poisson Noise*,
Journal of Imaging, volume 7(6), 2021, p. 99, DOI: 10.3390/jimaging7060099
(open access).

## Software requirements
ResPoND runs under MATLAB. It has been tested under MATLAB R2020a.

## Contents of the package
Here is the list of ResPoND files in alphabetical order:

- `respond.m`       : main function;
- `dir_est_hough.m` : function estimating the main direction of
                      directional images;
- `plot_line_rad.m` : function plotting a line along a specified direction
                      at the center of the current image.
- `ssim_index.m`    : implementation of the algorithm for calculating the
                      Structural SIMilarity (SSIM) index between two images
                      by Zhou Wang.

See the documentation inside each file for further details.

## Example of use
- `DTGVdemo.m`         : example of use for the case of directional images;
- `fibre_phantom.mat`  : phantom directional image;
- `TGVdemo.m`          : example of use for the case of general images;
- `smooth_phantom.mat` : phantom smooth image.


## License
[![GNU GPL v3.0](http://www.gnu.org/graphics/gplv3-127x51.png)](http://www.gnu.org/licenses/gpl.html)
