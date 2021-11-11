tic
%% prompt for selecting files, OME format (czi, lif, oif, nd2...)
clear all
[FileNames,PathName,FilterIndex] = uigetfile('*.czi;*.nd2;*.oif;*.lif','Select One or More Files','MultiSelect', 'on');
%%  reading raw data using OME scripts and performing analysis for every file selected
for file=1:length(FileNames)
    % in case a single file is selected...
    if length(FileNames)==1
    FileName=FileNames;
    else
    FileName=FileNames{file};
    end
   % clearvars -except FileNames PathName file incr
fullname=[PathName  FileName];
MetaData = GetOMEData(fullname);
pixelsize=MetaData.ScaleX; % pixel size in x and y
% prompt user to enter which is the first and last z-plane to be loaded
% from 3D+T data sets
prompt = {'Enter first slice of the range you want to load:','Enter last slice of the range you want to load:'};
titleinput = 'Input';
dims = [1 35];
definput = {num2str(1),num2str(MetaData.SizeZ)};
answer = inputdlg(prompt,titleinput,dims,definput);
firstLoad=str2double(answer{1});
lastLoad=str2double(answer{2});
reader = bfGetReader(fullname);
h = waitbar(0,'Processing Data ...');
totalframes = MetaData.SizeC * MetaData.SizeT*length(firstLoad:lastLoad);
clear maxXY maxXZ maxYZ
data=zeros(MetaData.SizeY,MetaData.SizeX,length(firstLoad:lastLoad),MetaData.SizeT,MetaData.SizeC);
framecounter = 0; 
zload2=1;
for zload=firstLoad:lastLoad
for tpt=1:MetaData.SizeT    
for channel = 1: MetaData.SizeC
       iPlane = reader.getIndex(zload - 1, channel -1, tpt - 1) + 1;       
        % get frame for current series
      data(:,:,zload2,tpt,channel) = bfGetPlane(reader, iPlane);
      framecounter = framecounter + 1;
      % update waitbar
      wstr = {'Reading Images: ', num2str(framecounter), ' of ', num2str(totalframes) };
      waitbar(framecounter / totalframes, h, strjoin(wstr))
     
end
end
 zload2=zload2+1;
end
close(h)
% end reading in z+t stacks
%% this section prompts user to select spheroid centre and then for every time point, 
% calculates the circular (azimuthal) average at varying radius from the
% spheroid centre...It is done for every z-plane separately, as in 2D
% averaging, and data is saved for every plane.
imagesc(squeeze(mean(mean(mean(data,3),4),5)));colormap(pink)
title('Please select the spheroid centre')
[x,y]=ginput(1);
close
% dtermine the maximum isotropic radius available for averaging, according
% to the spheroid centre selected and data extent in XY plane
dely=min(floor(y),floor(size(data,1)-y));
delx=min(floor(x),floor(size(data,2)-x));
del=min(delx,dely);
rangeR=0:(1/pixelsize):del;
miny=max(1,floor(y-del));
minx=max(1,floor(x-del));
clear tab
h = waitbar(0,'Processing Data ...');
totalframes = size(data,3)* MetaData.SizeT;
% for loops from first to last z plane to be loaded and over time
framecounter=0;
for zload=1:size(data,3)
for tpt=1:MetaData.SizeT
    % limiting data in XY plane according to the selected spheroid centre
    % and maximum extent available isotropically...
I=squeeze(data(miny:miny+2*del+1,minx:minx+2*del+1,zload,tpt,:));    
% defining the radii range
[X,Y] = meshgrid(-size(I,1)/2+0.5:size(I,1)/2-0.5,-size(I,1)/2+0.5:size(I,1)/2-0.5);
r = (X.^2 + Y.^2);
indices=find(r<(del^2));
[ii,jj]=ind2sub([size(r,1),size(r,2)],indices);
rcut=r(indices);
roundTol = 1e-11;
runique = consolidator(r(:),[],'mean',roundTol);
clear corrslice kBin
% for each plane and for each time point, extract the average intensity at
% given radius, using the function consolidator, which groups and averages
% entries with same value of r (runique defined above)
for i=1:size(I,3)
    for j=1:length(indices)
   corrslice(j,1)=squeeze(I(ii(j),jj(j),i));
    end
 [xg,yg] = consolidator(rcut,corrslice,@mean,roundTol);
 yg2=moving_average(yg,50);
 yg3=spline(xg.^0.5,yg2,rangeR);
    kBin(:,i) = yg3;   
end  

 
tab(:,zload,tpt,1)=rangeR;
tab(:,zload,tpt,2:size(kBin,2)+1)=kBin;

framecounter = framecounter + 1;
      % update waitbar
      wstr = {'Calculating Azimuthal Averages: ', num2str(framecounter), ' of ', num2str(totalframes) };
      waitbar(framecounter / totalframes, h, strjoin(wstr))
end
end
close(h)
% end of averaging

% storing the results into table and exporting into excel file...for each
% z-plane and time channel, save an excel sheet with averaged intensity
% and corresponding radii
rangeZ=firstLoad:lastLoad;
for zload=1:size(data,3)
for chan=1:size(tab,4)

tab2=squeeze(tab(:,zload,:,chan));
clear T1
T1=table(tab2);
writetable(T1,[PathName 'AA2Dperum_AzimuthAverage_' FileName(1:end-4) '.xlsx'],'Sheet',['z' num2str(rangeZ(zload)) 'ch' num2str(chan)]);

end
end

end % ends loop of analysis over every file selected
toc