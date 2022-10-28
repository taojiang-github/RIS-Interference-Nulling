function v =  update_v(v,A,y,gamma,max_iter)
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

theta = angle(v);
obj = nan(max_iter+1,1);
obj(1) = compute_obj(R,b,theta);
alpha = 0.3;
beta = 0.8;
for ii=1:max_iter
    grad_theta = 2*real(-1j*conj(v).*(R*v-b));
    theta_old = theta;
    lr = 1;
    while compute_obj(R,b,theta_old-lr*grad_theta)>(compute_obj(R,b,theta_old)-alpha*lr*norm(grad_theta)^2)
        lr = beta*lr;
        if lr<1e-20
            break;
        end
    end
    theta = theta_old-lr*grad_theta;
    v = exp(1j*theta);
    obj(ii+1) = compute_obj(R,b,theta);
    if ii>1 && abs(obj(ii)-obj(ii-1))<1e-4
%         plot(obj)
        break;
    end
end
end

function y = compute_obj(R,b,theta)
v = exp(1j*theta);
y = real(v'*R*v)-2*real(b'*v);
end