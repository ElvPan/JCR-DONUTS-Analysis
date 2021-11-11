% this script reads in one after the other multiple .mat files, which are
% outputs of 'DONUTS_3D_T_SimulationFormat.m'...It takes azimuth (3D+t)
% profiles and extracts diffusivity at each radius...
% define files 
clear all
[filenames,pathname,indx] = uigetfile('*.mat','MultiSelect','on');
%% read one by one and apply FTCS like analysis to get D(R)
pixelsize=0.1;
timesize=1;
for files=1:length(filenames)
clearvars -except filenames pathname files pixelsize timesize

if length(filenames)==1
load([pathname filenames]);
filename=filenames(1:end-4);
else
load([pathname filenames{files}]);
filename=filenames{files}(1:end-4);
end

%% extracting radial intensity profile for each time points/z-plane and performing the numerical integral (sum)...
% in order to calculate the time evolution of total intensity within a
% radius r from spheroid centre. 
% spheroid centre
dat=tab; %sheet name for data
r=(runique(rind).^0.5)*pixelsize; %sheet name for distance
t2=timesize*(0:size(dat,2)-1)';
clear dat2
% re-organise data by taking cumulative sum at each time pt, up to a
% defined radius...
for curt=1:size(dat,2)
    for curr=2:size(dat,1)
       curint=dat(1:curr,curt); %missing line
       id=find(~isnan(curint));
       curint=curint(id);
        Q = trapz(r(1:max(id)),curint); %variable correct to max(id)
       dat2(curr-1,curt)=Q;

    end  
end
%
%% Apply the FTCS like fit to the averaged intensity, at given radius, as function of time, in order to extract the 
% characteristic diffusion time for this radius...this diffusion time is
% further used to calculate the diffusion coefficient as described in the
% main text.
clear D Rsq adjRsq
ft = fittype('A*(1-exp(-tau*t))','independent','t','coefficients',{'A','tau'});
for cur=1:size(dat2,1);
datcur=dat2(cur,:);
datcur=datcur-datcur(1);
A0=mean(datcur(end-2:end));
tau0=log(2)/t2(max(find(datcur<A0)));
try
[fm,gof] = fit(t2,datcur',ft,'Startpoint',[A0,tau0],'Lower',[0 0],'Upper',[Inf inf]);
fitted=feval(fm,t2);
   
    Rsq(cur)=gof.rsquare;
    adjRsq(cur)=gof.adjrsquare;
    D(cur)=r(cur)^2/(log(2)/fm.tau); 
catch
D(cur)=NaN;
end
    
end
% export current file results into excel sheet
T1=table(r(2:end),D',Rsq',adjRsq','VariableNames',{'r','D','Rsq','adjRsq'});
writetable(T1,[pathname 'Diffusion Coefficient_' filename '.xlsx']);

end