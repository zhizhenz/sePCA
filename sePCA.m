function [ pc, eigval, coeff_WF, scale2, index, numpc ] = sePCA( coeffw, num_im, Binv, BD, c, R, L, scale )
%Description
% 
%
% Zhizhen Zhao Jan 2017

%shrinkage
%n = size(coeffw{1}, 2);
n = num_im;
max_freqs = size(coeffw, 1) - 1;
V = cell(max_freqs+1, 1);
D = cell(max_freqs+1, 1);
coeff_WF = cell(max_freqs + 1, 1);
scale2 = cell(max_freqs + 1, 1);
numpc = zeros(max_freqs+1, 1);
index = cell(max_freqs+1, 1);

for i = 1:max_freqs+1
    tmp = coeffw{i}(:, 1:n);
    if i == 1
        mean_coeff = mean(tmp, 2);
        tmp = bsxfun(@minus, tmp, mean(tmp, 2));
        p = size(tmp, 1);
        gamma = p/n;
    else
        p = size(tmp, 1);
        gamma = p/(2*n);
    end;

    C = 1/n*real(tmp*tmp');
    [ eig_vec, eig_val ] = eig(C); 
    eig_val = diag(real(eig_val));
    [eig_val, id] = sort(eig_val, 'descend');
    eig_vec = eig_vec(:, id);	
    %[ r, right_edge ] = permute_rankEst_sPCA(eig_val, tmp.', perm_alpha, 100);
    %r = r
    %r_list(i) = r;
    %right_edge_list(i) = right_edge;
    %right_edge_list(i) = (1 + sqrt(gamma))^2;
    %right_edge = right_edge
    shr_eig = op_norm_shrink(eig_val, gamma); %Shrink eigenvalues
    r = length(find(shr_eig>0));
    %r = min(r, r2)
    %r = r2
    %r = r2;
    shr_eig = shr_eig(1:r);
    eig_vec = eig_vec(:, 1:r);
    if r > 0 
    	S = eig_vec*(bsxfun(@times, eig_vec, shr_eig.'))'; %Covariance after eigenvalue shrinkage
    	S_recolored = Binv{i}*S*Binv{i}'; %Recolor the covariance matrix
    	[ V_tmp, D_tmp ] = eig(S_recolored);
    	D_tmp = real(diag(D_tmp)); 
        [ D_tmp , id ] = sort(D_tmp, 'descend');
        V_tmp = V_tmp(:, id(1:r));
        D_tmp = D_tmp(1:r);
        if scale == 1 %scaling
            Yvar = trace(BD{i});
            [~, c2] = standard_spiked_forward(shr_eig, gamma);
            s2 = 1 - c2;
            tau = (Yvar*shr_eig)./(D_tmp*p);
            alpha = zeros(r,1);
            for j = 1:r
                if c2( j ) > 0
                    alpha( j ) = max((1-s2( j )*tau( j ))/c2( j ), 0);
                else
                    alpha( j ) = 1;
                end
            end
            id_tmp = find(alpha>0);
            alpha = alpha(id_tmp);
            D_tmp = alpha.*D_tmp(id_tmp);
            V_tmp = V_tmp(:, id_tmp);
            r = length(id_tmp);
            scale2{i} = alpha;
        end
        if r > 0
            V{i} = V_tmp;
            D{i} = D_tmp;
            numpc(i) = r;
            index{i} = [(i)*ones(length(D_tmp), 1), [1:length(D_tmp)]'];
            %Wiener Filter the expansion coefficients
            coeff_recolor = Binv{i}*coeffw{i};
            S2 = V_tmp*diag(D_tmp)*V_tmp';
            coeff_WF{i} = (S2 / (S2 + BD{i}))*coeff_recolor;
            if i == 1
                mean_c = (BD{i} / (BD{i} + S2))*(Binv{i}*mean_coeff);
                coeff_WF{i} = bsxfun(@plus, coeff_WF{i}, mean_c);
            end;
         end;
    end;
end;
total_pc = sum(numpc)
D = cell2mat(D);
index = cell2mat(index);
[ D, D_id ] = sort(D, 'descend');
index = index(D_id, :);
[ fn ] = IFT_FB(R, c); %This generate inverse Fourier transform of Fourier-Bessel basis

cutoff = total_pc; %Get the top 15 eigenvalues
orig = floor(L/2) + 1;
pc = zeros(L, L, 2*cutoff);
eigval = zeros(2*cutoff, 1);
count = 1;
for i = 1:cutoff
    tmp = fn{index(i, 1)};
    size(tmp)
    tmp = reshape(tmp, (2*R)^2, size(tmp, 3));
    size(V{index(i, 1)})
    tmp = tmp*V{index(i, 1)}(:, index(i, 2));
    if index(i, 1) == 1
        pc(orig-R:orig+R-1, orig-R:orig+R-1, count) = reshape(real(tmp), 2*R, 2*R); %take the real part of the principal components
    else 
        pc(orig-R:orig+R-1, orig-R:orig+R-1, count) = reshape(sqrt(2)*real(tmp), 2*R, 2*R); %
    end;
    eigval(count) = D(i);
    count = count + 1;
    if index(i, 1) ~= 1
        pc(orig-R:orig+R-1, orig-R:orig+R-1, count) = reshape(sqrt(2)*imag(tmp), 2*R, 2*R); %take the imaginary part of the principal components
        eigval(count) = D(i);
        count = count + 1;
    end
end
pc = pc(:, :, 1:count-1);
eigval = eigval(1:count-1);
numpc = length(eigval);
