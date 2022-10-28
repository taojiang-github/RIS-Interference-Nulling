function v = mamnopt_func(A,p,alpha,sigma2_dB,v0)
        [num_user, num_elements_irs, ~] = size(A);
        R1 = zeros(num_user,num_elements_irs,num_elements_irs);
        R2 = zeros(num_user,num_elements_irs,num_elements_irs);
        for kk =1:num_user
            tmp=0;
            for jj=1:num_user
                a_kj = A(jj,:,kk)';
                tmp=tmp+a_kj*a_kj'*p(jj);
            end
            R1(kk,:,:) = tmp;
            a_kk = A(kk,:,kk)';
            tmp = tmp-a_kk*a_kk'*p(kk);
            R2(kk,:,:) = tmp;
        end
        problem.M = complexcirclefactory(num_elements_irs);
        problem.cost  = @(x) -compute_rate(A, x, p, alpha, sigma2_dB);
        problem.egrad = @(x) -compute_egrad(R1,R2,num_user, x, alpha, sigma2_dB);      % notice the 'e' in 'egrad' for Euclidean
        options.verbosity=0;
        options.stopfun = @mystopfun;
        options.maxiter = 5000;
%         v =exp(1j.*2*pi*rand(num_elements_irs,1));
%         v = v./abs(v);
        [v, ~, ~, ~] = conjugategradient(problem,v0,options);
%         [v, ~, ~, ~] = trustregions(problem,v,options);
end