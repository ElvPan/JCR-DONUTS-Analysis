clear all
%% read spheroid data of nuclei density
[FileName,PathName,FilterIndex] = uigetfile('*.czi;*.nd2;*.oif;*.lif','Select an OME data file (czi, lif, oif, nd2)');
fullname=[PathName  FileName];
MetaData = GetOMEData(fullname);

%% intialise reader of OME
reader = bfGetReader(fullname);
% find rectangle to crop 
iPlane = reader.getIndex(round(MetaData.SizeZ/2)- 1, 0, 0) + 1;
im=single(bfGetPlane(reader, iPlane));
 imagesc(im)
 title('Please select the rectangle containing spheroid')
 [im,rect]=imcrop;
%im=imcrop(im,rect);
close
tpt=1;
h = waitbar(0,'Processing Data ...');
totalframes =  MetaData.SizeZ ;
clear data
data=single(zeros(size(im,1),size(im,2),MetaData.SizeZ));
channel=1;
for zload=1: MetaData.SizeZ 
      iPlane = reader.getIndex(zload - 1, channel -1, tpt - 1) + 1;
      im=single(bfGetPlane(reader, iPlane));
      im=imcrop(im,rect);
      data(:,:,zload) = im;
      % update waitbar
      wstr = {'Reading Images: ', num2str(zload), ' of ', num2str(totalframes) };
      waitbar(zload / totalframes, h, strjoin(wstr))
end
close(h)
%% selecting a contour around spheroid in XY projection
imagesc(squeeze(max(data,[],3)))
title('scaled XY projection')
clear radius
for i=1:10
   h=imellipse;
   BW = createMask(h);
   radius(i)=(length(find(BW))/pi)^0.5;
end
close
averadi=mean(radius);
%% apply LoG filter to smooth data
[invconv,h2]=smoothSeries3D(data,[averadi/3 averadi/3 averadi/3],'log','n');
%% find nuclei using feature3D
featSize=averadi;
 inconv=feature3dMB(invconv,[featSize featSize featSize],[2*featSize 2*featSize 2*featSize],size(invconv),[1 1 1],2*featSize-1,2,0.3);
plot3((inconv(:,2)),(inconv(:,1)),(inconv(:,3)),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');

%xlabel('X (\mum)','FontSize',30)
%xlabel('X (\mum)','FontSize',20)
%ylabel('Y (\mum)','FontSize',20)
%zlabel('Z (\mum)','FontSize',20)
%set(gca,'FontSize',20)
%set(gcf,'Color',[1 1 1])
%% Nearest Neighbor distance calculation
D= squareform(pdist([inconv(:,1)*MetaData.ScaleX inconv(:,2)*MetaData.ScaleY inconv(:,3)*MetaData.ScaleZ], 'euclidean'));
minD=[];
for i=1:size(D,1)
distcur=min(D(find(D(:,i)>0),i));
minD(i)=distcur;
end
%% Plot histogram of  and save the figure
close all
nBins=round(length(minD)/20);
hist(log10(minD),nBins)
[yb,xb]=hist(log10(minD),nBins);
ylabel('Occurence (A.U.)','FontSize',25)
xlabel('log_{10} d_{N.N.} (\mum)','FontSize',25)
set(gca,'FontSize',25)
set(gcf,'Color',[1 1 1])
text(1,1.15*max(yb),['<d_{NN}>\pm\sigma_{d_{NN}}=' num2str(mean(minD)) '\pm' num2str(std(minD)) '\mum'],'FontSize',25)
%xlim([0 2])
ylim([0 1.2*max(yb)])
saveas(1,[PathName 'HistNNDist_' FileName(1:end-4) '.png'])
close

%% apply Delaunay, find convex Hull, spheroid volume and nucei density
DT = delaunayTriangulation(double([inconv(:,1)*MetaData.ScaleX inconv(:,2)*MetaData.ScaleY inconv(:,3)*MetaData.ScaleZ]));
[C,v2] = convexHull(DT);
volumeinmmcube=v2*(10^-9);
densityNuclei=size(inconv,1)/volumeinmmcube;
%% define center of spheroid and find distance of each nuclei to the center
xc=mean(inconv(:,1)*MetaData.ScaleX);
yc=mean(inconv(:,2)*MetaData.ScaleY);
zc=mean(inconv(:,3)*MetaData.ScaleZ);

distToC=((inconv(:,1)*MetaData.ScaleX-xc).^2+(inconv(:,2)*MetaData.ScaleY-yc).^2+(inconv(:,3)*MetaData.ScaleZ-zc).^2).^0.5;
%% export data into excel file
inconv(:,8)=distToC;
inconv(:,1:2)=inconv(:,1:2)*MetaData.ScaleX;
inconv(:,3)=inconv(:,3)*MetaData.ScaleZ;
T=array2table(inconv,'VariableNames',{'X','Y','Z','TotInt','Radius','PeakInt','FracVoxAboveThresh','DistToCenter'});
writetable(T,[PathName FileName(1:end-4) '.xlsx'],'Sheet','NucleiStats');
T=table(minD','VariableNames',{'MinDist'});
writetable(T,[PathName FileName(1:end-4) '.xlsx'],'Sheet','NNdist');
T=table([volumeinmmcube;densityNuclei],'RowNames',{'VolumeSpheroInmm','densityNucleiPermm'});
writetable(T,[PathName FileName(1:end-4) '.xlsx'],'WriteRowNames' ,true,'Sheet','DensityNuclei');

%% export boxplot of NN distance vs distance from spheroid centre
% binsize=20;
% numbins=round(max(distToC)/binsize);
% clear binsind MatricesToPlot
% for i=1:numbins;
%     lb=1+(i-1)*binsize;
%     ub=i*binsize;
%     ind1=find(distToC>lb);
%     ind2=find(distToC<ub);
%     ind=intersect(ind1,ind2);
%     binsind(ind)=i;
%     MatricesToPlot{i}=minD(ind);
%     ticklabel{i}=i*binsize;
% end
% 
% bh=notBoxPlot(MatricesToPlot)
% d = [bh.data];
% 
% set(d,'marker','o', 'markerfacecolor',[0.3 0.3 0.3], 'color', [0,0,0])
% set(d, 'markersize', 0.5);
% set(gca,'XTickLabel',ticklabel)
% ylabel('d_{NN} (\mum)','FontSize',25)
% xlabel('distance from center (\mum)','FontSize',25)
% set(gca,'FontSize',25)
% set(gcf,'Color',[1 1 1])
% 
% saveas(1,[PathName 'BoxplotNNDist_' FileName(1:end-4) '.png'])
% close
%% output some images before clearing data
close all
del=30;
for center=[del+1:200:size(data,3)-del-1]
Z=inconv(:,3)/MetaData.ScaleZ;
X=inconv(:,2)/MetaData.ScaleX;
Y=inconv(:,1)/MetaData.ScaleX;
ind1=find(Z<center+del);
ind2=find(Z>center-del);
ind=intersect(ind1,ind2);
imagesc(squeeze(mean(data(:,:,center-10:center+10),3)));colormap(gray)
hold on
text(size(data,2)-200,50,['#Nuc=' num2str(length(ind))],'Color',[1 1 1],'FontSize',20)
plot(X(ind),Y(ind),'oc')
set(gcf,'Position',[267    61   974   918])
saveas(1,[PathName 'Slice' num2str(center) FileName(1:end-4) '.png']) 
close 
end
%%
clear data invconv mask reader h2 
save([PathName 'workspace_' FileName(1:end-4) '.mat'])