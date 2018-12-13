function [ coeff ] = FBmain(data, R, basis, sample_points, num_pool)
%Description:
% Computes Fourier-Bessel expansion coefficients 'coeff', and filtered steerable PCA expansion coefficients 'sPCA_coeff'.
%Input:
%	data: image data
%	R: compact support radius
%       basis: FB basis
%       sample_points: quadrature points and weights for FB transform
%       num_pool: number of pools for parallel pools
%Output
%	coeff: Fourier-Bessel expansion coefficients
%Update: 10/15 Zhizhen Zhao

n = size(data, 3);
data2 = cell(num_pool, 1); %First set data into batches in cell structure
nb = floor(n/num_pool);
remain = n - nb*num_pool;
for i = 1:remain
    data2{i} = data(:, :, (nb+1)*(i-1)+1: (nb+1)*i);
end;
count = (nb+1)*remain;
for i = remain+1:num_pool
    data2{i} = data(:, :, count + (i-remain-1)*nb+1: count + (i-remain)*nb);
end; 
clear data;

[ coeff ]= FBcoeff_nfft(data2, R, basis, sample_points, num_pool);

end
