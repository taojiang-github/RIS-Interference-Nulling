function v =  update_v(v0,A,y,gamma,max_iter)
[num_user, num_elements_irs, ~] = size(A);
R = 0;
b = 0;
for kk =1:num_user
    tmp=0;
    for jj=1:num_user
        a_kj = A(jj,:,kk)';
        tmp=tmp+a_kj*a_kj';
    end
    R = R + abs(y(kk))^2*tmp;
    a_kk = A(kk,:,kk)';
    b = b+sqrt(1+gamma(kk))*y(kk)'*a_kk;
end
problem.M = complexcirclefactory(num_elements_irs);
problem.cost  = @(x) compute_obj(R,b,x);
problem.egrad = @(x) compute_egrad(R,b,x);      % notice the 'e' in 'egrad' for Euclidean
options.verbosity=0;
options.stopfun = @mystopfun;
options.maxiter = 50000;
% v0 =exp(1j.*2*pi*rand(num_elements_irs,1));
[v, ~, ~, ~] = conjugategradient(problem,v0,options);

end

function y = compute_obj(R,b,v)
y = real(v'*R*v)-2*real(b'*v);
end

function y = compute_egrad(R,b,v)
y = 2*R*v-2*b;
end