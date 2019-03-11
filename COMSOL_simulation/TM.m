% clc
clear

%% parameters
thickness =4*10^(-6);%8 micron
cutout = 56;%modes to discard
E_in_amp = 1;

lambda = 6*10^(-7);%600 nm
output_cut = 90*10^(-6);%90 micron
num_points = 499;
theta = linspace(-1.4741,1.4741,299);
num_modes = floor(output_cut/lambda*2-1);
modes = num_modes-cutout;
tmfp = 372.6*10^(-9);
mean_eig_theory =  tmfp/thickness;

%% read data
[filename, directory_name] = uigetfile('*.dat', 'Select Input');
fullname = fullfile(directory_name, filename);
data_input = load(fullname);
data_input = sortrows(data_input,2);
[filename, directory_name] = uigetfile('*.dat', 'Select Output');
fullname = fullfile(directory_name, filename);
data_output = load(fullname);
data_output = sortrows(data_output,2);
[dim_output,temp] = size(data_output);
num_ang = (temp-1)/2;
if num_ang ~=num_modes
    error('number of modes are not set up correct')
end
if dim_output ~= num_points
    error('number of spatial data points not right')
end
data_output_ang = zeros(dim_output,num_ang);
% data_input_ang = zeros(num_ang,num_ang);
data_input_ang = zeros(dim_output,num_ang);
% for m = 1: num_ang
%     data_input_ang(m,m) = E_in_amp;
% end

for i = 1: num_ang
    for j = 1 : dim_output
        data_input_ang(j,i) = data_input(j,i*2)+ data_input(j,i*2+1) .*sqrt(-1);
    end
end

for i = 1: num_ang
    for j = 1 : dim_output
        data_output_ang(j,i) = data_output(j,i*2)+ data_output(j,i*2+1) .*sqrt(-1);
    end
end

output_x_space = (data_output(:,2)-20)*10^(-6);
output_k_space = 2*pi./output_x_space;


%% k space
k_0 = 1/lambda*2*pi;
k_in = k_0*cos(theta);
k_out_unit = 2*pi/output_cut;
k_out_limit = pi/output_cut*dim_output;

%% Input Fourier transform
data_in_f = zeros(dim_output,num_ang);

for i = 1:num_ang
    data_in_f(:,i) = fftshift(fft(data_input_ang(:,i)));
end
data_in_f_corr = linspace(-k_out_limit,k_out_limit,dim_output);

for cutoff1 = 1:dim_output
    if data_in_f_corr(cutoff1)>=-k_0
        break
    end
end

for cutoff2= 1:dim_output
    if data_in_f_corr(cutoff2)>= k_0
        break
    end
end

data_in_f_prop = data_in_f(cutoff1:cutoff2-1,:);

input_power_k = zeros(1,num_ang);
for i = 1:num_ang
    input_power_k(i) = sum(abs(data_in_f_prop(:,i)));
end

%% output Fourier transform
data_out_f = zeros(dim_output,num_ang);

for i = 1:num_ang
    data_out_f(:,i) = fftshift(fft(data_output_ang(:,i)));
%     data_out_f(:,i) = fft(data_output_ang(:,i));
end
data_out_f_corr = linspace(-k_out_limit,k_out_limit,dim_output);
for cutoff1 = 1:dim_output
    if data_out_f_corr(cutoff1)>=-k_0
        break
    end
end
for cutoff2= 1:dim_output
    if data_out_f_corr(cutoff2)>= k_0
        break
    end
end
num_out_k = cutoff2-cutoff1+1;
data_out_f_prop = data_out_f(cutoff1:cutoff2-1,:);
[tm_dim,~] = size(data_out_f_prop);
k_space = data_out_f_corr;
output_power_k = zeros(1,num_ang);
for i = 1:num_ang
    output_power_k(i) = sum(abs(data_out_f_prop(:,i)));
end

%% TM for k space
tm_temp = zeros(tm_dim,tm_dim);

for i = 1: num_ang
tm_temp(:,i) = data_out_f_prop(:,i)/min(input_power_k);
end

tm = tm_temp(cutout/2+1:end-(cutout/2),cutout/2+1:end-(cutout/2));
[U,S,V] = svd(tm);
singular = zeros(num_ang-cutout,1);
for i = 1: num_ang-cutout
        singular(i) =  S(i,i);
end

eigenvalues = singular.^2;
trans_simu = mean(eigenvalues);
tmfp_simu = thickness*trans_simu;


figure(1)
plot(1:tm_dim-cutout,eigenvalues)
axis([0 300 0 max(eigenvalues)])

projected_amp = eigenvalues./(max(eigenvalues));

figure(2)
plot(1:tm_dim-cutout,projected_amp)
axis([0 300 0 max(projected_amp)])

figure(3)
imagesc(abs(tm))
colorbar

figure(4)
imagesc(angle(tm))
colorbar

[N,edges] = histcounts(eigenvalues,50);


figure(5)
% histogram(eigenvalues,40,'Normalization','probability')
histogram(eigenvalues,30)
% histogram(singular/mean(singular.^2),20)
% set(gca, 'Yscale',  'log')
xlim([0,1])

%% Eigenvalues Distribution
k = tmfp_simu/thickness/2;
tran_uniform = linspace(0,1,modes);
distribution = k./(tran_uniform.*sqrt(1-tran_uniform));
for non_zero =1 : modes;
    if eigenvalues(non_zero)< 0.0011;
        break
    end
end
non_zero_modes = non_zero-1;
zero_modes = num_modes-non_zero_modes;
one_modes = non_zero_modes - round(sum(distribution(2:end-1)));
if one_modes <0
    distribution(end) = nan;
    distribution(1) = num_modes-sum(distribution(2:end-1));
else
    distribution(end) = one_modes;
    distribution(1) = zero_modes;
end


bins = 60;
% distribution(2:end-1) = distribution(2:end-1)*non_zero_modes/bins;
distribution = distribution*modes/bins;

figure(6)
histogram(eigenvalues,bins)
hold on
plot(tran_uniform,distribution)
hold off
xlim([0,1])
ylim([0,60])

input = ones(243,1);
input1 = input;
input2 = input;
input2(23) = -1.5;
input2(43) = -0.5;
input2(46) = -1;
input2(53) = -1;
input2(33) = -1;
input2(13) = -1;
input2(52) = -1;
input2(5) = -1;
sum(abs(real(conj(V*input1).*((singular/mean(singular)).^2.*(V*input1)))))
sum(abs(real(conj(V*input1).*((singular/mean(singular)).^2.*(V*input2)))))
sum(abs(real(conj(V*input1).*((singular/mean(singular)).^2.*(V*input2)))))/sum(abs(real(conj(V*input1).*((singular/mean(singular)).^2.*(V*input1)))))
sum(abs(real(conj(V*input1).*(V*input2))))/sum(abs(real(conj(V*input1).*(V*input1))))