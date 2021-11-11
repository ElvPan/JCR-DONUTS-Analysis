% This script loads 3d+T spheroid data and for every time point averages
% the voxel intentisities at given radius from spheroid centre...Data is
% then exported 
clear all
% prompt for file selection...can be modifed from 'czi' to 
[FileName,PathName,FilterIndex] = uigetfile('*.czi;*.nd2;*.oif;*.lif','Select an OME data file (czi, lif, oif, nd2)');
fullname=[PathName  FileName];
MetaData = GetOMEData(fullname);
reader = bfGetReader(fullname);
sampleR=10; % how often to sample in um final radii rangeR
rescale=0.3; %reduce data size, 1 = original size.
incR=10; %increment to sample radial profile in pixels...bigger value speeds up the calculation per time pt, but samples profile sparsely

% loading 1st time point z-stack to be used to define some parameters below
tpt=1;
h = waitbar(0,'Processing Data ...');
totalframes =  MetaData.SizeZ ;
clear data
data=single(zeros(MetaData.SizeY,MetaData.SizeX,MetaData.SizeZ));
channel=1; %change if membrane channel is ''not=1
for zload=1: MetaData.SizeZ 
      iPlane = reader.getIndex(zload - 1, channel -1, tpt - 1) + 1;
      im=single(bfGetPlane(reader, iPlane));
      data(:,:,zload) = im;
      % update waitbar
      wstr = {'Reading Images: ', num2str(zload), ' of ', num2str(totalframes) };
      waitbar(zload / totalframes, h, strjoin(wstr))
end
close(h)
%rescale in 3d t=1 data
data=imresize3(data,rescale);
%% selecting the box around spheroid alone
imagesc(squeeze(data(:,:,round(size(data,3)/2))))
title('Please select the rectangle containing spheroid')
[im,rect]=imcrop;
close
clear ch1
for i=1:size(data,3)
    ch1(:,:,i)=imcrop(data(:,:,i),rect);
end
% selecting spheroid centre in XY and YZ plane
%% XY plane
imagesc(squeeze(nansum(ch1,3)));colormap(pink)
title('Please select the spheroid centre in XY plane')
[x,y]=ginput(1); 
close
dely=nanmin(floor(y),floor(size(ch1,1)-y));
delx=nanmin(floor(x),floor(size(ch1,2)-x));
% YZ plane
imagesc(squeeze(nansum(ch1,1)));colormap(pink)
title('Please select the spheroid centre in YZ plane')
[z,y]=ginput(1);
close
delz=nanmin(floor(z),floor(size(ch1,3)-z));
%% extracting pixel and deltaz size from meta data...prompts user to change if needed
pixelsize=MetaData.ScaleX; % pixel size in x and y
voxeldepth=MetaData.ScaleZ; %delta z in micormeters
prompt = {'Pixel size (\mum)','deltaZ (\mum)'};
dlg_title = 'voxel params';
num_lines = 1;
defaultans = {num2str(pixelsize),num2str(voxeldepth)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
pixelsize=str2double(answer{1}); 
voxeldepth=str2double(answer{2});
ratZtoX=voxeldepth/pixelsize;
%% adjust maximum radius of spheroid
delz=delz*ratZtoX;
del=min(delx,dely);
del=min(del,delz);
%% prompt for which channel to analyze, channel 2 (default for nanoparticle) is the one analysis is always applied on
prompt = {'Which channel should the analysis proceed on?'};
dlg_title = 'Channel for analysis...';
num_lines = 1;
defaultans = {num2str(2)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
channToAnalyse=str2double(answer{1}); 
%% Proceed with analysis...loading for every time point and every channel, the z-stack

tic
clear tabAve tab2Ave tabSum tab2Sum
delete(gcp('nocreate'))
parpool
for tpt=1:MetaData.SizeT
 clear data
h = waitbar(0,'Processing Data ...');
totalframes = MetaData.SizeZ * MetaData.SizeC;
framecounter=1;
data=single(zeros(MetaData.SizeY,MetaData.SizeX,MetaData.SizeZ, MetaData.SizeC));
    for zload=1:MetaData.SizeZ 
      for channel = 1: MetaData.SizeC
       iPlane = reader.getIndex(zload - 1, channel -1, tpt - 1) + 1;
       im=single(bfGetPlane(reader, iPlane));
       data(:,:,zload,channel) = im;
       wstr = {'Reading Images: ', num2str(framecounter), ' of ', num2str(totalframes) ' for t=' num2str(tpt) ' of total ' num2str(MetaData.SizeT)};
       waitbar(framecounter / totalframes, h, strjoin(wstr))
       framecounter=framecounter+1;
      end
    end
    close(h)
    %resizing both volumetric data to accommodate for available RAM
ch1a=squeeze(data(:,:,:,1));
ch1a=imresize3(ch1a,rescale);
Ia=squeeze(data(:,:,:,channToAnalyse));  % channel on which analysis is applied
Ia=imresize3(Ia,rescale);
% croping both data
clear I ch1
for i=1:size(ch1a,3)
    ch1(:,:,i)=imcrop(ch1a(:,:,i),rect);
    I(:,:,i)=imcrop(Ia(:,:,i),rect);
end
% define masks using average intensity of ch1
maskcells=ch1>0.5*mean(ch1(:));
se=strel('sphere',round(del/10));
maskcells=imclose(maskcells,se);
maskcells=single(maskcells);
maskcells(find(maskcells==0))=NaN;
I=I.*maskcells;
% circular averaging
[X,Y,Z] = meshgrid(1:size(I,2),1:size(I,1),1:ratZtoX:size(I,3)*ratZtoX);
X=X-round(x); 
Y=Y-round(y);
Z=Z-round(z)*ratZtoX;
r = (X.^2 + Y.^2 + Z.^2); clear X Y Z
r=round(r);
runique=unique(r(:));
% for every radius find the average or sum of intensity accross all angles in
% 3D...store in matrix and do for every time point
aveProf=[];
sumProf=[];
rind=1:incR:length(runique);
parfor i=1:length(rind) % parralel for loop to speed up
curind=find(r==runique(rind(i)));
sumProf(i)=nansum(I(curind));
aveProf(i)=nanmean(I(curind));
end
% resample (moving average) data to smooth it out and store for exporting

rangeR=0:sampleR/pixelsize:max((runique(1:incR:length(runique)).^0.5)*(1/rescale));%(0:sampleR/pixelsize:(del*(1/rescale)));
tabAve(:,tpt)=aveProf;
aveProf2=movmean(aveProf,20,'omitnan');
yNew=spline((runique(1:incR:length(runique)).^0.5)*(1/rescale),aveProf2,rangeR); %Add aveProf2 for moving average (smoothing of data)
tab2Ave(:,tpt)=yNew;

tabSum(:,tpt)=sumProf;
sumProf2=movmean(sumProf,20,'omitnan');
yNew=spline((runique(1:incR:length(runique)).^0.5)*(1/rescale),sumProf2,rangeR); %Add aveProf2 for moving average (smoothing of data)
tab2Sum(:,tpt)=yNew;
end
%Shutdown parpool
delete(gcp('nocreate'))
close
toc

%% export data
clear T1
T1=table((runique(1:incR:length(runique)).^0.5)*(1/rescale));
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet','R');
clear T1
T1=table(rangeR);
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet','rangeR');
%export profiles for analyzed channel
clear T1
T1=table(tabSum);
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet',['SumCh' num2str(channToAnalyse)]);
clear T1
T1=table(tab2Sum);
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet',['SumCh' num2str(channToAnalyse) 'resamp']);
clear T1
T1=table(tabAve);
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet',['AveCh' num2str(channToAnalyse)]);
clear T1
T1=table(tab2Ave);
writetable(T1,[PathName 'AA3D_' FileName(1:end-4) '.xlsx'],'Sheet',['AveCh' num2str(channToAnalyse) 'resamp']);

%% save workspace, excluding large matrices which are too big to save
clear data ch1 I X Y Z im reader r rcut mask indices jj kk ii 
save([PathName 'workspace_' FileName(1:end-4) '.mat'])