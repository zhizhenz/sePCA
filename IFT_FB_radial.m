function [ fn ] = IFT_FB_radial(r, c, R)

%%% Description
%The function computes the inverse Fourier transform of Fourier-Bessel Basis
%   This evaluate Bessel radial functions for an image of size
%   Table bessel.mat gives the R_ns (3rd column) for n (1st column) and s (2nd column)
% Input:  
%       1. c: band limit
%       2. R: compact support radius in real domain
% Output: 
%       1. fn: IFT of FB basis, stored in cell structure. cell number corresponds to the angular index.     
% Zhizhen Zhao 04/2015

load bessel.mat
bessel = bessel(bessel(:, 4)<= 2*pi*c*R, :);
%bessel = bessel(bessel(:, 4)<= pi*c*R, :);
k_max = max(bessel(:, 1));
fn = cell(k_max+1,1);
for i = 1:k_max+1
    bessel_k = bessel(bessel(:, 1) == i-1, :);
    l = size(bessel_k, 1);
    tmp = zeros(length(r), l);
    for lr = 1:l
    	Jk = besselj(i-1, 2*pi*c*r(:));
    	Rkq = bessel_k(lr, 3);
    	f_r = 2*c*sqrt(pi)*(-sqrt(-1))^(i-1)*(-1)^lr*Rkq*Jk./((2*pi*c*r).^2 - Rkq^2);
    	tmp(:, lr) = f_r;
    end;
    fn{i} = tmp;
end;

