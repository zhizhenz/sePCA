%need to add nfft matlab path

%read data
[ data ] = cxi2datamatrix( 'condor_20000.cxi' ); 
%replace 'condor_20000.cxi' with your data filename

% choose average pixel intensity to scale to 0.01, noise type is Poisson
average_intensity = 0.01;
data = average_intensity*data/mean(data(:)); %rescale data
nim = size(data, 2);
L = size(data, 1);
high = max(data(:));
low = 0;
c = 0.1; %set bandlimit of the image

for iter = 1:10 %10 repititions
    iter
    noisy = poissrnd(data);
    for i = 1:length(nim)
        ni = nim(i);
        tic_start_coeff = tic;
        [ coeffw, Binv, BD, c, R, mean_est ] = FBwhiten(noisy(:, :, 1:nim(i)), c); %coeffw is the steerable PCA expansion coefficients of whitened image
        toc_coeff = toc(tic_start_coeff);
        
        tic_start_sePCA = tic;
        [ pc, eigval, coeff, scale, index, num_pc ] = sePCA(coeffw, nim(i), Binv, BD, c, R, L, 1); %use steerable ePCA to denoise the coefficients 
        toc_sePCA = toc(tic_start_sePCA); 
        
        tic_start_denoise = tic;
        [ denoised ] = denoise_poisson(coeff, edge, c, R, nim(i)); %use denoised coeff to reconstruct the images. 
        toc_denoise = toc(tic_start_denoise);
        err = denoised - data(:, :, 1:nim(i));
        mse = mean(err(:).^2, 1);
        filename = sprintf('denoise_nim%d_iter%d.mat', nim(i), iter);
        denoised = denoised(:, :, 1:6);
        save(filename, 'pc', 'eigval', 'index', 'num_pc', 'c', 'denoised', 'mse', 'toc_coeff', 'toc_sePCA', 'toc_denoise');
    end;
end;
