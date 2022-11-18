%%% Image processing Electrochemical fluorescence microscopy
%%% Thomas Cochard - Harvard University - Aziz lab - 11/18/2022
%%% This version works with the exemple available on the github repository
%%% Thomcochard/Ectrochemical-Fluorescence-Microscopy

clear all; close all;
lw=2; Ms=20;

%% Load directory

dir='C:\Users\thcoc\Dropbox (Harvard University)\BES porous electrodes proposal 2018-11\Researchers files\Cochard\3_Experiments\20211209 - AQDS - commercial\Analysis\Code publication'; % Directory with the set of images
addpath(genpath(dir));

%% Inputs

C=10;                  % AQDS Concentiontration (mM)
E=[0 50 100 175];      % Applied potential (mV)

%% Load Video

I=loadtiff('Video.tif'); I=im2double(I);                  % Read video files
I_Ref=loadtiff('Ref.tif'); I_Ref=im2double(I_Ref);        % Read 100 SOC image

%% Removing impurities

p=99;

for k=1:size(I,3)
A = reshape(I(:,:,k),[],1); 
Y_vid=prctile(A,[(100-p)/2 (100+p)/2],'all');
lo_vid = find(A<Y_vid(1)); hi_vid = find(A>Y_vid(2));
A(lo_vid)=nan; A(hi_vid)=nan; 
I_C(:,:,k)= reshape(A ,size(I_Ref));
end

B = reshape(I_Ref,[],1);
Y_img=prctile(B,[(100-p)/2 (100+p)/2],'all');
lo_img = find(B<Y_img(1)); hi_img = find(B>Y_img(2));
B(lo_img)=nan; B(hi_img)=nan;
I_C_Ref= reshape(B,size(I_Ref));

figure(1)
subplot 221
imagesc(I_Ref(:,:),'AlphaData',~isnan(I_Ref(:,:))); axis image; colormap jet ;colorbar; %caxis([0 1]);
title('Raw reference image')
subplot 222
h=histogram(I_Ref(:,:),100,'Normalization','probability'); axis square
subplot 223
imagesc(I_C_Ref(:,:),'AlphaData',~isnan(I_C_Ref(:,:))); axis image; colormap jet ;colorbar; %caxis([0 1]);
title('Impurity corrected reference image')
subplot 224
h=histogram(I_C_Ref(:,:),100,'Normalization','probability');  axis square

figure(2)
for k = 1: size(I,3)
subplot (2,4,2*k-1)
imagesc(I_C(:,:,k),'AlphaData',~isnan(I_C(:,:,k))); axis image; colormap jet ;colorbar; %caxis([0 1]);
title(['Impurity corrected E=' num2str(E(k)) 'mV'])
subplot (2,4,2*k)
h=histogram(I_C(:,:,k),100,'Normalization','probability');  axis square; box on;
end


%% Normalization: dividing by reference image

I_Hat(:,:,:)=(I_C(:,:,:)-I_C(:,:,1))./(I_C_Ref(:,:));

figure(3)
for k = 1:size(I,3)
subplot (2,4,2*k-1)
imagesc(I_Hat(:,:,k),'AlphaData',~isnan(I_Hat(:,:,end))); axis image; colormap jet ;colorbar; caxis([0 1]);
title(['Normalized E=' num2str(E(k)) 'mV'])
subplot (2,4,2*k)
hold on
h=histogram(I_Hat(:,:,k),'BinEdges',0:0.01:1); axis square; box on;
xlabel('$$\hat{I}$$','Interpreter','Latex')
ylabel('# pixel')
xlim([0 1])
box on
end

 
 %% Removing fibers

Thresh_Fiber=0.2; % Intensity threshold value for estimated for the fiber
A = reshape(I_Hat(:,:,end),[],1); 
Fiber_Idx=find(A<Thresh_Fiber);

for k=1:size(I,3)
B(:,k)= reshape(I_Hat(:,:,k),[],1); 
B(Fiber_Idx,k)=NaN;
I_Hat_C(:,:,k)=reshape(B(:,k),size(I_Ref));
end

figure(4)
for k = 1: size(I,3)
subplot (2,4,2*k-1)
imagesc(I_Hat_C(:,:,k),'AlphaData',~isnan(I_Hat_C(:,:,k))); axis image; colormap jet ;colorbar; caxis([0 1]);
title(['Fiber corrected E=' num2str(E(k)) 'mV'])
subplot (2,4,2*k)
hold on
h=histogram(I_Hat_C(:,:,k),'BinEdges',0:0.01:1); axis square; box on;
xlabel('[H_2AQDS]/[AQDS_0]')
ylabel('# pixel')
end


 %% Average

for k=1:size(I,3)
I_Hat_C_Av(k)=mean(I_Hat_C(:,:,k),'all','omitnan');
end
 

figure(5)
hold on
plot(E,I_Hat_C_Av,'x-','LineWidth',lw)
xlabel('Overpotential (mV)')
ylabel('Final corrected Intensity (a.u)')

