%%
% this script reads in one after the other multiple excel files, which are
% outputs of 3D_T or 2D_T averaging script...It takes the f(R,t) intensity 
% profiles and extracts diffusivity at each radius using FTCS like
% analysis...
% define files
clear all
[files,path,indx] = uigetfile('*.xlsx','MultiSelect','on');
%read one by one and apply FTCS like analysis to get D(R)
for f=1:length(files)
    clearvars -except files path f 
    if length(files)==1
    file=files;
    else
    file=files{f};
    end
num = importdata([path file]);

dat=(num.data.AveCh2resamp) %sheet name for data, users must change "z1ch3" and must be all the same
r=(num.data.rangeR(:,1)); %sheet name for distance, users must change "z1ch1" and must be all the same
t2=(0:size(dat,2)-1)';
clear dat2
%% extracting radial intensity profile for each time points/z-plane and performing the numerical integral (sum)...
% in order to calculate the time evolution of total intensity within a
% radius r from spheroid centre. 
% spheroid centre
for curt=1:size(dat,2) % loop over time
    for curr=2:size(dat,1) %loop over radii
       curint=dat(1:curr,curt); % current intensity up to this radius and at this time point
       Q = trapz(r(1:curr),curint); % numerical integral
       dat2(curr-1,curt)=Q; %storing the outputs
    end  
end

clear D  Rsq adjRsq
for cur=1:size(dat2,1);
datcur=dat2(cur,:);
datcur=datcur-datcur(1);
A0=mean(datcur(end-2:end));
    tau0=log(2)/t2(max(find(datcur<A0)));
    ft = fittype('A*(1-exp(-tau*t))','independent','t','coefficients',{'A','tau'});
    [fm,gof] = fit(t2,datcur',ft,'Startpoint',[A0,tau0],'Lower',[0 0],'Upper',[Inf inf]);
    fitted=feval(fm,t2);
   
    Rsq(cur)=gof.rsquare;
    adjRsq(cur)=gof.adjrsquare;
    D(cur)=r(cur)^2/(log(2)/fm.tau);    
end
% saving data into excel sheet
T1=table(D',Rsq',adjRsq');
writetable(T1,[path 'Diffusion Coefficient_' file(1:end-4) '.xlsx']);

end