% clear
close all

load('perfectPi.mat')
threshold = 0.34;
Pi_unprocessed = perfectPi;

Sample = 2;
[imn1, imn2]    = size(qzgf);
imn11           = Sample * imn1;    % size of the image with zero paddings
imn22           = Sample * imn2;

start1          = round(1 + (imn11-150)/2);  % position of the true image
end1            = start1 + 150 - 1;
start2          = round(1 + (imn22-150)/2);
end2            = start2 + 150 - 1;

Pi_proccessed = zeros(imn11,imn22);
Pi_proccessed(start1:end2,start2:end2) = perfectPi(start1:end2,start2:end2);

Pi_proccessed(Pi_proccessed<threshold) = 0;


figure(1)
imagesc(qx(1,:)/1000,qy(:,1)/1000, Pi_unprocessed(start1:end1,start2:end2));
% title('Reconstructer Pi')
ax = gca;
ax.XTick = -3:1:3;
ax.YTick = -3:1:3;
xlabel('x (mm)')
ylabel('y (mm)')
axis image
colorbar

figure(2)
imagesc(qx(1,:)/1000,qy(:,1)/1000, Pi_proccessed(start1:end1,start2:end2));
%title('Reconstructer Pi after thresholding')
ax = gca;
ax.XTick = -3:1:3;
ax.YTick = -3:1:3;
xlabel('x (mm)')
ylabel('y (mm)')
axis image
colorbar


figure(3)
imagesc(qx(1,:)/1000,qy(:,1)/1000,(qzgf-min(min(qzgf)))/(max(max(qzgf))-min(min(qzgf))))
ax = gca;
ax.XTick = -3:1:3;
ax.YTick = -3:1:3;
xlabel('\Deltax (mm)')
ylabel('\Deltay (mm)')
axis image
colorbar
axis image


sx = linspace(min(min(qx)), max(max(qx)), 700);
sy = linspace(min(min(qy)), max(max(qy)), 700);
piSample = im2double(sample(151:850,81:780));
figure(4)
imagesc(qx(1,:)/1000,qy(:,1)/1000,piSample/max(max(piSample)))
axis image
ax = gca;
ax.XTick = -3:1:3;
ax.YTick = -3:1:3;
xlabel('x (mm)')
ylabel('y (mm)')
colorbar
colormap gray

figure(5)
imagesc(itx/1000,ity/1000,ingf./max(max(ingf)));
ax = gca;
ax.XTick = -3:1:3;
ax.YTick = -3:1:3;
xlabel('\Deltax (mm)')
ylabel('\Deltay (mm)')
axis image
hold on
plot3(ix/1000, iy/1000, intensity, 'o');
hold off

