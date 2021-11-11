clear variables
close all
%% Experimental conditions
prompt = {'What are particles radii (nm)?','Define how many time points data has?','Define the time increment (h)?'};
dlg_title = 'Define particle size and time points in data';
defaultans = {'30',num2str(24),num2str(0.5)};
answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
particleSize = str2double(answer{1});              %Size of particle (in micrometeres)
totTimePts=str2double(answer{2});
timeInc=str2double(answer{3});
timePointsInData=[timeInc:timeInc:timeInc*totTimePts];
%%
% prompt to open a file
[FileNames,PathName,FilterIndex] = uigetfile('*.xlsx','Select one or more .xlsx files?','MultiSelect','on');
%FileNames = convertCharsToStrings(FileNames);                               % [NEW] Need to use strings for length calculations

 % in case a single file is selected...
    if ~iscell(FileNames)
        numfiles=1;
    else iscell(FileNames)
        numfiles=length(FileNames);
    end

for file=1:numfiles
    clearvars -except particleSize FileNames PathName file numfiles timePointsInData
    % in case a single file is selected...
    if ~iscell(FileNames)
        FileName=FileNames;
    else iscell(FileNames)
        FileName=FileNames{file};
    end
    
    %% Load experimental data
    addpath(genpath(pwd))           %Add directory
    dataFile=[PathName  FileName];
    dataFile=dataFile(1:end-5);
    inputData = xlsread(sprintf('%s.xlsx',dataFile),'AveCh2resamp'); %Load data
    coordinates=xlsread(sprintf('%s.xlsx',dataFile),'rangeR');       %Load co-ordinate data
    coordinates=coordinates(1:size(inputData,1));                    %Distance from centre of cell spheroid
    %% Plot raw experimental data
    figure(2); hold on; plot(coordinates,inputData,'g'); xlabel('Distance from spheroid centre'); ylabel('Arbitrary fluorescence');
    
    %% Set-up model parameters
    D = 1.38e-23*310.15/(3*pi*1e-3*particleSize*1e-9)*1e12; %Particle diffusivity (microns squared per second) in free solution
    
    [~,locationMax] = max(inputData);                       %Choose index of maximum distance of data considered (automatic)
    locationMax=15;                                         %Choose index of maximum distance of data considered (manual)
     nDataTimePoints = numel(timePointsInData);              %[NEW] Number of time points in experiment.
    Parameters.maxPoint = coordinates(mode(locationMax));   %Maximum distance of data considered
    Parameters.nDataPoints = mode(locationMax);             %Number of data points considered per time point
    Parameters.tEnd = max(timePointsInData)*60*60;          %Final time of model solution (seconds)
    Parameters.spheroidRadius = 250;                        %Size of spheroid (microns)
    Parameters.domainMultiple = 7;                          %Number of multiples of spheroid in domain
    Parameters.domainWidth = Parameters.spheroidRadius*Parameters.domainMultiple;   %Size of model domain
    Parameters.nGridPoints = 1000;                          %Number of nodes in model solution
    Parameters.nTimePoints = 2400;                          %Number of timesteps in model solution
    Parameters.dt = Parameters.tEnd/Parameters.nTimePoints; %Timestep length in model solution
    Parameters.outputTimes = ceil(3600/Parameters.dt)*timePointsInData; %[NEW] The timesteps in the simulation to output
    dataScaling = 10;                                       %Scaling of arbitrary fluorescence
    Parameters.edgeConcentration = 0.75;                    %NP concentration at the end of the domain (equivalent to 0.75*dataScaling A.U.)
    
    %% Set-up initial distribution of NPs in solution
    Parameters.initialCondition = zeros(Parameters.nGridPoints,1);
    Parameters.initialCondition(ceil(1*Parameters.nGridPoints/Parameters.domainMultiple):end) = Parameters.edgeConcentration;   %3000 A.U. outside spheroid, 0 A.U. inside
    
    %% Remove any fluorescence artefacts by removing fluorescence at initial timepoint (where there should be 0 NPs in the spheroid)
    if timePointsInData(1) == 0 %[NEW] Only if first data point is captured at 0 hours.
        inputData = inputData(1:Parameters.nDataPoints,:)-inputData(1:Parameters.nDataPoints,1);
    else
        inputData = inputData(1:Parameters.nDataPoints,:);
    end
    
    %% Set-up model solution parameters
    Parameters.grid = linspace(1,Parameters.domainWidth,Parameters.nGridPoints);    %Model domain
    Parameters.dGrid = Parameters.grid(2)-Parameters.grid(1);                       %Space between nodes in model domain
    Parameters.boundaryCondition = 'fixed';                                         %Type of boundary condition (constantly 0.75*dataScaling A.U.)
    opt = optimset('MaxIter',1000,'MaxFunEvals',1000);                              %Curve-fitting parameters
    
    %% Fit model solution to the data
    diffusionParameters = lsqnonlin(@(y) NumericalRadialDiffusionSimulationFunction(@(x) D*(double(x<=Parameters.spheroidRadius).*(x/Parameters.spheroidRadius).^y(1)*(1-y(2))+double(x>Parameters.spheroidRadius)+double(x<=Parameters.spheroidRadius)*y(2)),Parameters)-inputData/dataScaling,[0,0.01],[-10,0],[1e10,1],opt);
    
    %% Set-up model domain and diffusion function to visualise solution
    spheroidGrid = linspace(0,1.2*Parameters.spheroidRadius,300);                   %Model domain for solution visualisation
    diffusivity = D*(double(spheroidGrid<=Parameters.spheroidRadius).*(spheroidGrid/Parameters.spheroidRadius).^diffusionParameters(1)* ...
        (1-diffusionParameters(2))+double(spheroidGrid>Parameters.spheroidRadius)+double(spheroidGrid<=Parameters.spheroidRadius)*diffusionParameters(2)); %Diffusion function
    solutionPlot = NumericalRadialDiffusionSimulationFunction(@(x) D*(double(x<=250).*(x/Parameters.spheroidRadius).^diffusionParameters(1)*(1-diffusionParameters(2))+double(x>Parameters.spheroidRadius)+double(x<=Parameters.spheroidRadius)*diffusionParameters(2)),Parameters); %Calculate predicted fluorescence based on fit parameters
    
    colourOrder = cool(24); %Choose colour map for solution
    
    %% Plot predicted fluorescence and data
    close
    for i = 1:nDataTimePoints
        subplot(1,2,1); plot(coordinates(1:Parameters.nDataPoints),inputData(:,i),'col',colourOrder(i,:),'linewidth',2); ylim([0 dataScaling]); hold on;
        subplot(1,2,2); plot(coordinates(1:Parameters.nDataPoints),solutionPlot(:,i)*dataScaling,'col',colourOrder(i,:),'linewidth',2); ylim([0 dataScaling]); hold on;
        
    end
    subplot(1,2,1); ylabel('Arbitrary fluorescence'); xlabel('Distance from spheroid centre');title('Input data')
    set(gca,'FontSize',14)
    subplot(1,2,2); ylabel('Arbitrary fluorescence'); xlabel('Distance from spheroid centre');title('Numerical Solution')
    set(gcf,'Color',[1 1 1])
    set(gca,'FontSize',14)
    set(gcf,'Position',[300 300 800 300])
    saveas(gcf,sprintf('%s_out1.png',dataFile))
    %% Plot diffusivity
    close
    figure; plot(diffusivity,'LineWidth',3); ylabel('Diffusivity'); xlabel('Distance from spheroid centre'); 
    set(gcf,'Color',[1 1 1])
    set(gca,'FontSize',20)
    saveas(gcf,sprintf('%s_out2.png',dataFile))
    %% Save data output
    save(sprintf('%s_Output.mat',dataFile));
    %% Export data to excel
    clear T1
    T1=table(diffusivity);
    writetable(T1,[PathName 'NumFit_' FileName(1:end-4) '.xlsx'],'Sheet','Diffusivity');
    clear T1
    T1=table(spheroidGrid);
    writetable(T1,[PathName 'NumFit_' FileName(1:end-4) '.xlsx'],'Sheet','Radius (um)');
end