%% salecting file to load...it is the output of azimuthal averaging script, which is normally saved as excell file.
% this can be modified to accomodate other type of input files...
clear all
[file,path,indx] = uigetfile('*.xlsx');
num = importdata([path file]);
%% extracting radial intensity profile for each time points/z-plane and performing the numerical integral (sum)...
% in order to calculate the time evolution of total intensity within a
% radius r from spheroid centre. 
% sphjeroid centre
dat=(num.data.AveCh2resamp) %sheet name for data
r=(num.data.rangeR(:,1)); %sheet name for distance
t2=(0:size(dat,2)-1)';
clear dat2
% re-organise data by taking cumulative sum at each time pt, up to a
% defined radius...
for curt=1:size(dat,2) % loop over time
    for curr=2:size(dat,1) %loop over radii
       curint=dat(1:curr,curt); % current intensity up to this radius and at this time point
       Q = trapz(r(1:curr),curint); % numerical integral
       dat2(curr-1,curt)=Q; %storing the outputs
    end  
end
%% Apply the FTCS like fit to the averaged intensity, at given radius, as function of time, in order to extract the 
% characteristic diffusion time for this radius...this diffusion time is
% further used to calculate the diffusion coefficient as described in the
% main text.
clear D  Rsq adjRsq
% do for loop fover every radius value
for cur=1:size(dat2,1);
datcur=dat2(cur,:);
datcur=datcur-datcur(1);
A0=mean(datcur(end-2:end));
    % FTCS like analysis and fitting
    tau0=log(2)/t2(max(find(datcur<A0)));
    ft = fittype('A*(1-exp(-tau*t))','independent','t','coefficients',{'A','tau'});
    [fm,gof] = fit(t2,datcur',ft,'Startpoint',[A0,tau0],'Lower',[0 0],'Upper',[Inf inf]);
    fitted=feval(fm,t2);
   % store goodness-of-fit outputs
    Rsq(cur)=gof.rsquare;
    adjRsq(cur)=gof.adjrsquare; 
    D(cur)=r(cur)^2/(log(2)/fm.tau);  % store D extracted  
end
T1=table(D',Rsq',adjRsq'); % save all outputs into a table
%% write the table into an excell sheet
writetable(T1,[path 'Diffusion Coefficient_' file(1:end-4) '.xlsx']);
