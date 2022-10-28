function [channel_irs_user,location_user_set]=generate_channel_Sparse(num_irs_elements, num_user, irs_Nh, ...
    num_samples,num_path)
channel_irs_user = nan(num_samples,num_irs_elements,num_user);
location_user_set = nan(num_samples,num_user,3);
phi_lim = [-pi/2,pi/2];
for ii =1:num_samples
    i1 = mod(0:num_irs_elements-1,irs_Nh);
    i2 = floor((0:num_irs_elements-1)/irs_Nh);
    for kk=1:num_user
        a_irs_user = 0;
        for jj=1:num_path
            pathloss = randn/sqrt(2)+1j*randn/sqrt(2);
            phi_1 = phi_lim(1)+(phi_lim(2)-phi_lim(1)) * rand;
            phi_2 = phi_lim(1)+(phi_lim(2)-phi_lim(1)) * rand;
            tmp1 = cos(phi_2)*sin(phi_1);
            tmp2 = sin(phi_2);
            a_irs_user = a_irs_user+ pathloss*exp(1j*pi*(i1.*tmp1+i2.*tmp2)).';
        end
        channel_irs_user(ii,:,kk) =  a_irs_user/sqrt(num_path);
    end
end
end

