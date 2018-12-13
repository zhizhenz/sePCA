function [ mean_image ] = denoise_mean(fn, mean_coeff, L0, R)
% Description
% This function gives mean and denoised images in real space
%
% Input: 
%	U: eigenvectors of C^{(k)}, saved in cell array
% 	fn: inverse Fourier transform of the Fourier- Bessel function in real space on a Cartesian grid.
%	mean_coeff: mean Fourier-Bessel expansion coefficients for k=0
%	Coeff: sPCA expansion coefficients (Wienter filtered).
%	L0: original image size
%	n_max: number of images you want to denoise
%	 
% Output:
%	mean_image: mean image in real space on a Cartesian grid of size L0 x L0
%	denoised: denoised images in real space on a Cartesian grid L0 x L0
% Update 10/15 Zhizhen Zhao

L = 2*R;
%Computes eigen images, need output from IFT_FB.m.

%Original image size
origin = floor(L0/2) + 1;

tmp = fn{1};
tmp = reshape(tmp, L^2, size(tmp, 3));
mean_Im = reshape(tmp*mean_coeff, L, L);
mean_Im = real(mean_Im);
mean_image = zeros(L0);
mean_image(origin-R:origin+R-1, origin-R:origin+R-1) = mean_Im;

