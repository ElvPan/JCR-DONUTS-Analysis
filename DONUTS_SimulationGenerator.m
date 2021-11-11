%%
clear all
path=uigetdir([],'Select directory to which simulation data will be saved');
path=[path '\'];
%% prompt some inputs for simulations
prompt = {'Image size in x (pixels)','Image size in y (pixels)','Nmber of Z slices(pixels)','# of temporal stacks','Pixel size (um)','Frame time (s)','PSF lateral e^-2 radius (pixels)','PSF axial e^-2 radus (pixels)','Ratio Spheroid to full box','Cell density fraction','Cell radius to spheroid ratio','# of particles simulated','Save every Nth','Diffusion coefficient inside cells (um^2/s)','Diffusion coefficient outside cells (um^2/s)','Probability of particles entering cells (0 to 1)','Probability of particles exiting cells (0 to 1)','What file type (mp4,avi) would you like to export?','What frame rate do you want to save video (fps)?'};
dlg_title = 'Enter Simulation Parameters';
defaultans = {num2str(300),num2str(300),num2str(300),num2str(1000),num2str(0.1),num2str(1),num2str(2),num2str(4),num2str(0.5),num2str(0.2),num2str(0.1),num2str(100000),num2str(10),num2str(0.001),num2str(0.1),num2str(0.9),num2str(0.1),'mp4','10'};
answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
numPops=str2double(answer{1});
sizeX=str2double(answer{1}); 
sizeY=str2double(answer{2});
sizeZ=str2double(answer{3});
sizeT=str2double(answer{4});
pixelsize=str2double(answer{5});
timesize=str2double(answer{6});
PSFSizeXY=str2double(answer{7});
PSFZ=str2double(answer{8});
SpheroidBoxFraction=str2double(answer{9});
CellDensityFrac=str2double(answer{10});
CellRadiusToSpheroidRatio=str2double(answer{11});
numParticles=str2double(answer{12});
saveNth=str2double(answer{13});
Din=str2double(answer{14});
Dout=str2double(answer{15});
Pin=str2double(answer{16});
Pout=str2double(answer{17});
MovType=answer(18);
FrameRate=str2double(answer(19));
%% PSF extent and radius of spheroid defined from inputs
omega=[PSFSizeXY PSFSizeXY PSFZ];
largeR=SpheroidBoxFraction*sqrt(round(sizeX/2)^2+round(sizeY/2)^2+round(sizeX/2)^2);
% intialise all particles positions
xCoor=ceil(sizeX*rand(1,numParticles));
yCoor=ceil(sizeY*rand(1,numParticles));
zCoor=ceil(sizeZ*rand(1,numParticles));
% define mask in 3D which will be spheroid volume
dim=3;
d=inline('sqrt(sum(p.^2,2))-1','p');
[p,t]=distmeshnd(d,@huniform,CellDensityFrac,[-ones(1,dim);ones(1,dim)],[]);
p=p*largeR;
close
ptRads=CellRadiusToSpheroidRatio*largeR*ones(1,size(p,1));
mask=zeros(sizeX,sizeY,sizeZ);
for j=1:size(p,1)
       [masknow]=ellipsoid2(sizeX,sizeY,sizeZ,p(j,1),p(j,2),p(j,3),ptRads(j),ptRads(j),ptRads(j));
       mask=mask+double(masknow);
end   
% define spheroid boundary by populating surface with boundary points...to
% be used to define which particles are within/outside of spheroid...
N_new=0;
clear X Y Z
for Phi = 0:pi/10:2*pi
    for Theta=0:pi/10:pi
     N_new = N_new + 1;
X(N_new) =round(sizeX/2)+largeR*sin(Theta)*cos(Phi);
Y(N_new) =round(sizeY/2)+largeR*sin(Theta)*sin(Phi);
Z(N_new) =round(sizeZ/2)+largeR*cos(Theta);
    end
end
tess = convhulln([X' Y' Z']);
%% simulate particle movement via diffusion

%define which particles are inside and which are out of spheroid
in = inhull([xCoor' yCoor' zCoor'],[X' Y' Z'],tess);
% initialixe simulation propereties such as initial particle positions
    population.xCoor = xCoor(~in);
    population.yCoor = yCoor(~in);
    population.zCoor = zCoor(~in);
    population.xCoorDisplay=xCoor(~in);
    population.yCoorDisplay=yCoor(~in);
    population.zCoorDisplay=zCoor(~in);
    population.Din=Din;
    population.Dout=Dout;
    population.Pin=Pin;
    population.Pout=Pout;
imsize=[sizeX sizeY sizeZ];
% initial 3D spheroid data generated
cubeT=single(zeros(sizeX,sizeY,sizeZ,round(sizeT/saveNth)));
A=accumarray(ceil([  population.xCoorDisplay;population.yCoorDisplay;population.zCoorDisplay]'),ones(length(population.zCoorDisplay),1), [sizeX, sizeY, sizeZ]);
B = imgaussfilt3(A,omega);
cubeT(:,:,:,1)=single(B);
%move particles around, taking into account various simulation parameters,
%and update positions at each time step...
j=2;
for k=1:sizeT  
[population] = simul8trMovementDomain3D(population,timesize,pixelsize,sizeX,sizeY,sizeZ,mask);
if mod(k,saveNth)==0
A=accumarray(ceil([  population.xCoorDisplay;population.yCoorDisplay;population.zCoorDisplay]'),ones(length(population.zCoorDisplay),1), [sizeX, sizeY, sizeZ]);
B = imgaussfilt3(A,omega);
cubeT(:,:,:,j)=single(B);
j=j+1;
[k j]
end
end
% save 3D+T data for further processing
clear masknow in A B 
save([path 'Dout' num2str(Dout) 'Din' num2str(Din) 'Pout' num2str(Pout) 'Pin' num2str(Pin) 'Rspheroid' num2str(largeR) 'CellDensityFrac' num2str(CellDensityFrac) 'CellRadiusToSpheroidRatio' num2str(CellRadiusToSpheroidRatio) '.mat'],'-v7.3')
%% generate a movie of simulated data and save
name=['Dout' num2str(Dout) 'Din' num2str(Din) 'Pout' num2str(Pout) 'Pin' num2str(Pin) 'Rspheroid' num2str(largeR) 'CellDensityFrac' num2str(CellDensityFrac) 'CellRadiusToSpheroidRatio' num2str(CellRadiusToSpheroidRatio)];

if exist('writerObj')
close(writerObj);
end
% define movie file
if strcmp(MovType,'mp4')
writerObj = VideoWriter([path name '.mp4'],'MPEG-4');
elseif strcmp(MovType,'avi')
writerObj = VideoWriter([path name '.avi'],'Uncompressed AVI');
else
    disp('Acceptable movie type is mp4 or avi')
end
writerObj.FrameRate =FrameRate;
open(writerObj);
close all
f = waitbar(0,'Exporting movie of simulated data...');
for i=1:size(cubeT,4)
V1=squeeze(cubeT(:,:,:,i));
V2 = imadjustn(V1,[0 1],[]);
V2=V2.^0.3;
 h1 = figure;set(h1, 'Visible', 'off');
slice(V2,size(V2,2)/2,size(V2,1)/2,size(V2,3)/2)
colormap gray
shading flat
axis equal
grid off

xlabel('X (pixel)','FontSize',20)
ylabel('Y (pixel)','FontSize',20)
zlabel('Z (pixel)','FontSize',20)
hold on
plot3(1:size(V1,1),(size(V1,2)/2)*ones(1,size(V1,1)),(size(V1,3)/2)*ones(1,size(V1,1)),'-','Color',[1 0 0])
hold on
plot3((size(V1,1)/2)*ones(1,size(V1,2)),1:size(V1,2),(size(V1,3)/2)*ones(1,size(V1,2)),'-','Color',[0 1 0])
hold on
plot3((size(V1,1)/2)*ones(1,size(V1,3)),(size(V1,2)/2)*ones(1,size(V1,3)),1:size(V1,3),'-','Color',[0 0 1])
set(gca,'FontSize',20)
set(gcf,'Color',[1 1 1])
set(gcf,'Position',[427   -75   937   767])
frame = getframe(1);
writeVideo(writerObj,frame);
close(h1)
waitbar(i/size(cubeT,4),f,'Exporting movie of simulated data...');

end
close(f);
close(writerObj);
