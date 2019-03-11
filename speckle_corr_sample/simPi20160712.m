
close all
clc
clear
drawplots = 1;
reconPi = 0;
phaseRetrieval = 0;
% load('pi255by255.mat');
load('perfectPi.mat')
Sample = 2;
% qzgf = qzgf.*qCorr2(ingf);
[imn1, imn2]    = size(qzgf);
imn11           = Sample * imn1;    % size of the image with zero paddings
imn22           = Sample * imn2;
Pi_padded            = zeros(imn11,imn22);          % image with zero paddings
Phi_measured = zeros(imn11,imn22);
Phi_measured_proc = zeros(imn11,imn22);
start1          = round(1 + (imn11-imn1)/2);  % position of the true image
end1            = start1 + imn1 - 1;
start2          = round(1 + (imn22-imn2)/2);
end2            = start2 + imn2 - 1;



qzgf = (qzgf-min(min(qzgf)))/(max(max(qzgf))-min(min(qzgf)));
qzgf(qzgf<0)=0;
qzgf = sqrt(qzgf)/max(max(sqrt(qzgf)));
Measured = qzgf;
%Measured = qz;
imagesc(qzgf)
Phi_measured(start1:end2,start2:end2)    = Measured;
%qzgf_proccessed = zeros(imn1,imn2);
qzgf_proccessed = Measured;
qzgf_proccessed(qzgf_proccessed<0.13) = 0.05;
% qzgf_proccessed = qzgf_proccessed -0.13;
% qzgf_proccessed = (qzgf_proccessed-min(min(qzgf_proccessed)))/(max(max(qzgf_proccessed))-min(min(qzgf_proccessed)));
imagesc(qzgf_proccessed);
axis image
sizeWinX = 75;
sizeWinY = 75;
windFunc1X = zeros(imn1);
windFunc1X((imn1+1)/2-sizeWinX:(imn1+1)/2+sizeWinX) = tukeywin(2*sizeWinX+1,.35);
windFunc1X = windFunc1X/max(windFunc1X);
windFunc2X = zeros(imn2);
windFunc2X((imn2+1)/2-sizeWinY:(imn2+1)/2+sizeWinY) = tukeywin(2*sizeWinY+1,.45);
windFunc2X = windFunc2X/max(windFunc2X);
windFunc = windFunc1X*windFunc2X';
Phi_measured_proc(start1:end2,start2:end2)    = qzgf_proccessed.*windFunc;
%Phi_measured_proc(start1:end2,start2:end2)    = qzgf_proccessed;
%Phi_measured_proc = Phi_measured_proc.*(bohmanwin(imn11)*bohmanwin(imn22)');
imagesc(Phi_measured_proc)
axis image
% figure(1)
% imagesc(Pi_pad);
axis image
% title('Pi simulated, larger than actual experiment');
% Pi_spec = fftshift(fft2(fftshift(Pi_pad)));
% figure(2)
% subplot(1,2,1)
% imagesc(abs(Pi_spec));
% axis image
% title('magnitude of spectrum of simulated Pi')
% subplot(1,2,2)
% imagesc(abs(angle(Pi_spec)));
% axis image
% title('phase of spectrum of simulated Pi')


Phi = ifftshift(ifft2(ifftshift(abs(Pi_spec).^2)));
Phi = Phi/max(max(abs(Phi)));
% figure(3)
% subplot(1,2,1)
% imagesc(abs(Phi))
% axis image
% title('magnitude of Phi from simulated Pi')
% subplot(1,2,2)
% imagesc(abs(angle(Phi)))
% axis image

figure(4)
subplot(1,2,1)
imagesc(Phi_measured)
axis image
title('magnitude of Phi measured')
subplot(1,2,2)
imagesc(Phi_measured_proc)
axis image
title('magnitude of Phi measured and processed')

figure(5)
subplot(1,2,1)
imagesc(Phi_measured_proc)
axis image
title('magnitude of Phi measured and processed')
subplot(1,2,2)
imagesc(abs(Phi))
axis image
title('magnitude of Phi')

PhiSyn = abs(Phi);
SpecSq = fftshift(fft2(fftshift(PhiSyn))); % imaginary part really small
Pi_spec_sim_abs = sqrt(abs(SpecSq));
% figure(6)
% subplot(1,2,1)
% imagesc(abs(Pi_spec));
% axis image
% title('spectrum of simulated Pi')
% subplot(1,2,2)
% imagesc(Pi_spec_sim_abs)
% axis image
% title('if simulated Phi has only zero phase, invert')



sizeWin = 180;
windFunc1XPad = zeros(imn11);
% windFunc1XPad((imn11+1)/2-sizeWin:(imn11+1)/2+sizeWin) = tukeywin(2*sizeWin+1,1);
windFunc1XPad((imn11+1)/2-sizeWin:(imn11+1)/2+sizeWin) = hamming(2*sizeWin+1);
windFunc1XPad = windFunc1XPad/max(windFunc1XPad);
windFunc2XPad = zeros(imn22);
% windFunc2XPad((imn22+1)/2-sizeWin:(imn22+1)/2+sizeWin) = tukeywin(2*sizeWin+1,1);
windFunc2XPad((imn22+1)/2-sizeWin:(imn22+1)/2+sizeWin) = hamming(2*sizeWin+1);
windFunc2XPad = windFunc2XPad/max(windFunc2XPad);
windFuncPad = windFunc1XPad*windFunc2XPad';

SpecSqPR = fftshift(fft2(fftshift(Phi_measured_proc)));
Pi_spec_sim_abs_PR = sqrt(abs(SpecSqPR));
Pi_spec_sim_abs_PR = abs(Pi_spec_sim_abs_PR)/(max(max(abs(Pi_spec_sim_abs_PR))));
%For_PR  = Pi_spec_sim_abs_PR.*windFuncPad;
%For_PR  = Pi_spec_sim_abs_PR.*(hann(imn11)*hann(imn22)');
%For_PR  = Pi_spec_sim_abs_PR.*(hamming(imn11)*hamming(imn22)');
For_PR  = Pi_spec_sim_abs_PR.*(blackman(imn11)*blackman(imn22)');
% For_PR = zeros(imn11,imn22);
% For_PR(start1:end2,start2:end2)   = Pi_spec_sim_abs_PR(start1:end2,start2:end2) .*(hamming(imn1)*hamming(imn2)');
figure(7)
subplot(1,2,1)
imagesc(For_PR)
axis image
title('spectrum of real Pi')
subplot(1,2,2)
imagesc(abs(Pi_spec));
axis image
title('spectrum of simulated Pi')


%[mask, RfacR, RfacF, RFD] = OSS_samplecode(For_PR,[90,90],2000,0.9,1,Pi_pad);
%[mask, RfacR, RfacF, RFD] = OSS_samplecode_triag(For_PR,[90,90],2000,0.9,1,Pi_pad);
%[mask, RfacR, RfacF, RFD] = OSS_beta(For_PR,[100,100],2000,1,Pi_pad);


RFD = (RFD- min(min(RFD)))/(max(max(RFD))-min(min(RFD)));
temp = RFD;
temp(temp<0.2) = 0.0;
temp = temp.*mask;
imagesc(temp)
axis image




