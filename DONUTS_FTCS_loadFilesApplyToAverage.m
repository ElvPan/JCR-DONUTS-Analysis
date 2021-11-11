%% this script allows user to load average (azimuth) intensity profiles
% and further combine them (average them) and apply FTCS like analysis to
% their average...useful for noisier profiles
% select files
clear all
[files,path,indx] = uigetfile('*.xlsx','MultiSelect','on');
datAll=[];
% load files and average int. profiles at equivalent radii R
for f=1:length(files)
    clearvars -except files path f 
    file=files{f};
num = importdata([path file]);
dat=(num.data.z1ch3) %sheet name for data, requires user changes
datAll(:,:,f)=dat;
end
% actual averaging
dat=mean(datAll,3);
%% extracting radial intensity profile for each time points/z-plane and performing the numerical integral (sum)...
% in order to calculate the time evolution of total intensity within a
% radius r from spheroid centre. 
% sphjeroid centre
r=(num.data.z1ch1(:,1)); %sheet name for distance, requires user changes
t2=(0:size(dat,2)-1)';
clear dat2
% re-organise data by taking cumulative sum at each time pt, up to a
% defined radius...
for curt=1:size(dat,2)
    for curr=2:size(dat,1)
       curint=dat(1:curr,curt);
       Q = trapz(r(1:curr),curint);
       dat2(curr-1,curt)=Q;
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
%saving data into an excel file
T1=table(D',Rsq',adjRsq');
writetable(T1,[path 'DiffusionCoefficient_Average_' file(1:end-4) '.xlsx']);

