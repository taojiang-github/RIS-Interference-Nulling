function [signal_power,interfence_power, sir] = compute_SIR(A,B,v)
    [num_user, num_irs_elements, ~] = size(A);
    cond = abs(sum(abs(v))-num_irs_elements)<1e-3;
    assert(cond,'phase shifts v not valid')
    
    signal_power = nan(num_user,1);
    interfence_power = nan(num_user,1);
    sir = nan(num_user,1);
    for kk=1:num_user
        a_kk = A(kk,:,kk).';
        b_kk = B(kk,kk);
        signal_power(kk) = abs(a_kk.'*v + b_kk)^2;
        denominator = 0;
        for jj=1:num_user
            if jj ~=kk
                a_kj = A(jj,:,kk).';
                b_kj = B(jj,kk);
                denominator = denominator+abs(a_kj.'*v+ b_kj)^2;
            end
        end
        interfence_power(kk) = denominator;
        sir(kk) = signal_power(kk)/interfence_power(kk);
    end
end