function [ coeffw, Binv, BD, c, R, mean_est ] = FBwhiten(data)
%Description
%Get Fourier expansion coefficients after whitening     
%    Input:
%	 data: image dataset
%	 c: bandlimit
%	 R: radius of the signal in the real space
%    Output:
%	 coeffw: FB expansion coefficents after whitening
%        Binv: matrices for recoloring

n = size(data, 3); %number of images
L = size(data, 1); %size of the images
[ x, y ] = meshgrid(-floor(L/2):ceil(L/2)-1, -floor(L/2):ceil(L/2)-1); %x-y grid points
r = sqrt(x.^2 + y.^2);

R = est_support(data, 0.99); %Support size
c = 0.5; %Bandlimit is 0.5
mean_img = mean(data, 3);
n_r = ceil(4*c*R);
[ basis, sample_points ] = precomp_fb( n_r, R, c );
[ coeff ] = FBmain(mean_img, R, basis, sample_points, 1); %use 8 parallel pools 
mean_coeff = mean(coeff{1}, 2);
[ fn ] = IFT_FB(R, c);
[ mean_est ] = denoise_mean(fn, mean_coeff, L, R);
mean_est(r>R) = 0;

data = reshape(data, L^2, n);
mean_est = mean_est(:);
mean_est(mean_est<0) = 0; %the reconstructed mean_est might be negative
g = 1./sqrt(mean_est); %for prewhitening
g(isinf(g)) = 0;
whiten = bsxfun(@times, data, g);
whiten = reshape(whiten, L, L, n);

[ coeffw ] = FBmain(whiten, R, basis, sample_points, 8);
clear whiten
%compute the recoloring matrices Binv
max_ang_freqs = size(coeffw, 1) - 1;
Binv = cell(max_ang_freqs + 1, 1);
[ r2, w2 ] = lgwt(length(r), 0, R);
[ fn ] = IFT_FB_radial(r2, c2, R);
w2 = r2.*w2;

load bessel.mat
bessel = bessel(bessel(:, 4)<= 2*pi*c*R, :);
bessel_k = bessel(bessel(:, 1) == 0, :);
l = size(bessel_k, 1);
ift_fb = zeros(length(r2), l);
for lr = 1:l
    Jk = besselj(0, 2*pi*c*r2);
    Rkq = bessel_k(lr, 3);
    ift_fb(:, lr) = 2*c*sqrt(pi)*(-1)^lr*Rkq*Jk./((2*pi*c*r2).^2 - Rkq^2);
end;
f = real(ift_fb*mean_coeff);
f(f<0) = 0;
g2 = sqrt(f);
g2 = g2.*w2;
g3 = f.*w2;

BD = cell(max_ang_freqs+1, 1);
for k = 1:max_ang_freqs + 1
    tmp2 = bsxfun(@times, fn{k}, g2);
    tmp3 = bsxfun(@times, fn{k}, g3);
    Binv{k} = 2*pi*fn{k}'*tmp2;
    BD{k} = 2*pi*fn{k}'*tmp3;
end;
