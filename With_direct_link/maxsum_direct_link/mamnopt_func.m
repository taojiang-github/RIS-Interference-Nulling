function v = mamnopt_func(A,B,p,alpha,sigma2_dB,v0)
        [num_user, num_elements_irs, ~] = size(A);
        R1 = zeros(num_user,num_elements_irs,num_elements_irs);
        R2 = zeros(num_user,num_elements_irs,num_elements_irs);
        c1 = zeros(num_user,num_elements_irs,1);
        c2 = zeros(num_user,num_elements_irs,1);
        d1 = zeros(num_user,1);
        d2 = zeros(num_user,1);
        for kk =1:num_user
            tmp=0;
            tmp1=0;
            tmp2=0;
            for jj=1:num_user
                a_kj = A(jj,:,kk)';
                b_kj = B(jj,kk)';
                tmp=tmp+a_kj*a_kj'*p(jj);
                tmp1=tmp1+a_kj*b_kj'*p(jj);
                tmp2 = tmp2+abs(b_kj)^2;
            end
            R1(kk,:,:) = tmp;
            c1(kk,:)=tmp1;
            d1(kk)=tmp2;
            a_kk = A(kk,:,kk)';
            b_kk = B(kk,kk)';
            tmp = tmp-a_kk*a_kk'*p(kk);
            tmp1 = tmp1-a_kk*b_kk';
            tmp2=tmp2-abs(b_kk)^2;
            R2(kk,:,:) = tmp;
            c2(kk,:)=tmp1;
            d2(kk)=tmp2;
        end
        problem.M = complexcirclefactory(num_elements_irs);
        problem.cost  = @(x) -compute_rate(A, B, x, p, alpha, sigma2_dB);
        problem.egrad = @(x) -compute_egrad(R1,R2,c1,c2,d1,d2, num_user, x, alpha, sigma2_dB);      % notice the 'e' in 'egrad' for Euclidean
        options.verbosity=0;
        options.maxiter=5000;
% %         options.stopfun = @mystopfun;
        
%         v =exp(1j.*2*pi*rand(num_elements_irs,1));
%         v = v./abs(v);
        [v, cost, info, options] = conjugategradient(problem,v0,options);
%         [v, ~, ~, ~] = trustregions(problem,v,options);
end