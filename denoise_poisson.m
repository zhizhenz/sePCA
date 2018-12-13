function [ denoised ] = denoise_poisson(coeff, L0, R, n_max)
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

c = 0.5
[ fn ] = IFT_FB(R, c);
max_ang_freqs = size(coeff, 1) - 1; %Should be the maximum angular frequency
L = 2*R;

%Original image size
origin = floor(L0/2) + 1;
tmp1 = fn{1};
tmp1 = reshape(tmp1, L^2, size(tmp1, 3));
denoised = zeros(L0, L0, n_max);

tmp2 = tmp1*coeff{1}(:, 1:n_max);
tmp2 = reshape(tmp2, L, L, n_max);
I = real(tmp2);

for k = 1:max_ang_freqs
    if size(coeff{k+1},1)~=0 
        tmp = fn{k+1};
        tmp = reshape(tmp, L^2, size(tmp, 3));
        tmp2_pos = tmp*coeff{k+1}(:, 1:n_max);
	tmp_2 = 2*real(tmp2_pos);
        I = I + reshape(tmp_2, L, L, n_max);
    end;
end;
denoised = zeros(L0, L0, n_max);
denoised(origin-R:origin+R-1, origin-R:origin+R-1, :) = real(I);
