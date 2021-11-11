%% this script prompts user to select a .mat file created during the simulation of spheroids uptake,
% and extract azimuthal average at every time point.
% select and read file
clear all
[FileName,PathName,FilterIndex] = uigetfile('*.mat');
fullname=[PathName  FileName];
load(fullname);
%% define increment in R to be sampled
sampleR=10; % how often to sample in um final radii rangeR
incR=10; %increment to sample radial profile in pixels...bigger value speeds up the calcualtion per time pt, but samples profile sparsely
ch1=squeeze(cubeT(:,:,:,1));
%% XY plane select speroid centre
imagesc(squeeze(nanmean(ch1,3)));colormap(pink)
title('Please select the spheroid centre in XY plane')
[x,y]=ginput(1); 
close
dely=nanmin(floor(y),floor(size(ch1,1)-y));
delx=nanmin(floor(x),floor(size(ch1,2)-x));
% YZ plane select speroid centre
imagesc(squeeze(nanmean(ch1,1)));colormap(pink)
title('Please select the spheroid centre in YZ plane')
[z,y]=ginput(1);
close
delz=nanmin(floor(z),floor(size(ch1,3)-z));
%% pixel and deltaz size
voxeldepth=pixelsize;
ratZtoX=voxeldepth/pixelsize;
%% adjust maximum radius
delz=delz*ratZtoX;
del=min(delx,dely);
del=min(del,delz);
[masknow]=ellipsoid2(sizeX,sizeY,sizeZ,0,0,0,largeR,largeR,largeR);
%%  or each time point do azimuth average
tic
clear tab tab2
delete(gcp('nocreate'))
parpool
for tpt=1:size(cubeT,4)
I=squeeze(cubeT(:,:,:,tpt));
masknow=single(masknow);
masknow(find(masknow==0))=NaN;
I=I.*masknow;
% circular averaging
[X,Y,Z] = meshgrid(1:size(I,2),1:size(I,1),1:ratZtoX:size(I,3)*ratZtoX);
X=X-round(x); 
Y=Y-round(y);
Z=Z-round(z)*ratZtoX;
r = (X.^2 + Y.^2 + Z.^2); clear X Y Z
r=round(r);
runique=unique(r(:));
aveProf=[];
rind=1:incR:length(runique);
parfor i=1:length(rind)
curind=find(r==runique(rind(i)));
aveProf(i)=nanmean(I(curind));
%[i tpt]
end 
tab(:,tpt)=aveProf;
end
%Shutdown parpool
delete(gcp('nocreate'))
close
toc

%% export data into excel sheet
clear T1
T1=table(runique(rind).^0.5);
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet','R');
clear T1
T1=table(tab);
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet','AveIntProfile');

%% save workspace
clear cubeT ch1 I X Y Z im reader r rcut mask indices jj kk ii 
save([PathName  'AA3D_' FileName(1:end-4) '.mat'])
%% basic plot
close all
mycol=jet(size(tab,2));
for i=1:size(tab,2)
plot((runique(rind).^0.5),tab(:,i),'Color',mycol(i,:),'LineWidth',2)
hold on
end
set(gca,'FontSize',20)
xlabel('R (pixel)','FontSize',20)
ylabel('Azimuthal average at R (A.U.)','FontSize',20)
set(gcf,'Color',[1 1 1])
saveas(1,[PathName FileName(1:end-4) '.png'])