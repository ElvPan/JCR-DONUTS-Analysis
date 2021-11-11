%%
% this script reads in a .mat file, which is
% output of 'DONUTS_3D_T_SimulationFormat.m'...
% It takes azimuth (3D+T) profiles and extracts diffusivity at each radius...
% define and read in the file
clear all
[filename,pathname,indx] = uigetfile('*.mat');
load([pathname filename]);
%% proceed with apply FTCS like analysis to get D(R) 
% pixelsize=0.1;
% timesize=1;
dat=tab; %sheet name for data
r=(runique(rind).^0.5)*pixelsize; %sheet name for distance
t2=timesize*(0:size(dat,2)-1)'; 
clear dat2
for curt=1:size(dat,2)
    for curr=2:size(dat,1)
       curint=dat(1:curr,curt);
       id=find(~isnan(curint));
       curint=curint(id);
       Q = trapz(r(1:max(id)),curint);
       dat2(curr-1,curt)=Q;
    end  
end
%
%% Apply the FTCS like fit to the averaged intensity, at given radius, as function of time, in order to extract the 
% characteristic diffusion time for this radius...this diffusion time is
% further used to calculate the diffusion coefficient as described in the
% main text.
clear D
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
%% output data into excel
T1=table(r(2:end),D',Rsq',adjRsq','VariableNames',{'r','D','Rsq','adjRsq'});
writetable(T1,[pathname 'Diffusion Coefficient_' filename '.xlsx']);

