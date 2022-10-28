function [eta,tmp_eta ]= compute_eta(A,B)
    [num_user, ~, ~] = size(A);
    tmp_eta = zeros(num_user,num_user);
    for kk =1:num_user
        for jj=1:num_user
             a_kj = A(jj,:,kk).';
             b_kj = B(jj,kk);
             tmp_eta(kk,jj) = abs(b_kj)/norm(a_kj,1);
        end
    end
    indx = 1:num_user+1:num_user*num_user;
    tmp_eta0 = tmp_eta;
    tmp_eta0(indx) = [];
    eta = max(tmp_eta0);
end