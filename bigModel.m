%final script to model big+ experiments
%CMJ 20200702
clear all


set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultAxesLineWidth',2)

blue=[43 172 226]./256;
orange=[248 149 33]./256;

figureCounter=1;

%import growth data to fit u0 (exp constant)

YPDdata=readtable('YPD_data.txt');
nRows=length(YPDdata{:,1});
k=1;

for i=1:nRows
    if length(YPDdata{i,1}{1})>2
        if strcmp(YPDdata{i,1}{1}(1:3),'<d>')
            toParse=YPDdata{i,1};
            toParse=toParse{1};
            toKeep=strsplit(toParse,'>');
            toKeep=toKeep{2};
            output=strsplit(toKeep,'<');
            output=str2num(output{1});
            if ~isempty(output)
                data(k)=output;
                k=k+1;
            end
        end
    end
end
data=data';


%parse growth data into vectors
vTime=data(1:886); %time in mins
vTime=vTime(1:662);

for i=1:9
    timeOffset=886*2;
    growthOffset=662;
    startIdx=timeOffset+growthOffset*(i-1)+1;
    endIdx=startIdx+growthOffset-1;
    vGrowth(:,i)=data(startIdx:endIdx);
end

%vMax from EC analysis to set u0/u1 ratio
vMaxNaive=[0.00283,0.00300,0.00333];
vMaxBIG=[0.00367,0.00333,0.00350,0.00333,0.00383,0.00350];



for i=1:9
    [yPred, tPred, m(i)] = fitGrowth2(vGrowth(:,i),vTime);
end

%fit aging params (d0,d1) from chron aging exp

vTime=xlsread('20171109agingDataToFit.xlsx','A15:A18');
dataMat=xlsread('20171109agingDataToFit.xlsx','B15:J18');

naiveMean=mean(dataMat(:,1:3),2);
naiveStd=std(dataMat(:,1:3),0,2);%./sqrt(3);
BIGmean=mean(dataMat(:,4:9),2);
BIGstd=std(dataMat(:,4:9),0,2);%./sqrt(6);

%convert to hrs
vTime=vTime*24;

%fit exponential decay constants d0(2) and d1(2) (linear fit to X vs logY)
X=[ones(length(vTime),1) vTime];
d0=X\log10(naiveMean);
d1=X\log10(BIGmean);


%import competition data
vTime2=(0:9)*13;
dataMat2=xlsread('20171109agingDataToFit.xlsx','B24:I33');

naiveMean2=mean(dataMat2(:,5:8),2);
naiveStd2=std(dataMat2(:,5:8),0,2)./sqrt(4);
BIGmean2=mean(dataMat2(:,1:4),2);
BIGstd2=std(dataMat2(:,1:4),0,2)./sqrt(4);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,3)
hold on
plot(vTime2,log10(naiveMean2./BIGmean2),'b')
plot(vTime2,log10(BIGmean2./naiveMean2),'r')
ylim([-3 3])
axis square
ylabel('log_{10}(ratio)')
xlabel('generations')
title('Fig. 2F')

for i=1:4
    plot(vTime2,log10(dataMat2(:,i)./naiveMean2),'r','LineWidth',1)
end

for i=5:8
    plot(vTime2,log10(dataMat2(:,i)./BIGmean2),'b','LineWidth',1)
end



%initialize
params.x00=1;
params.x10=1;

params.u0=mean(m(1:3))*60;
params.u1=params.u0*mean(vMaxBIG)/mean(vMaxNaive);  %vMax ratio from experiment
params.d0=-d0(2); %death constants
params.d1=-d1(2);

%t1=12; %hrs
params.t1=500/2380*40; %set time in exp phase
params.t2=40-params.t1;
%assume no difference in lag phase (~8 hrs)

params.N=19; %passages

%example traces
[t,BIGfrac]=oscillate(params);

%full model
subplot(2,3,3)
hold on
plot(t./(params.t1+params.t2)*6.5,log10(1./BIGfrac),'--b')
plot(t./(params.t1+params.t2)*6.5,log10(BIGfrac),'--r')



%first generate distribution of death rate values
M=1000;
for i=1:M
    
    naiveToFit=naiveMean+randn*naiveStd;
    BIGtoFit=BIGmean+randn*BIGstd;
    
    X=[ones(length(vTime),1) vTime];
    d0sim(i,:)=X\log10(naiveToFit);
    d1sim(i,:)=X\log10(BIGtoFit);
end

%plot distributions
subplot(2,3,1)
histfit(-real(d0sim(:,2)),8)
hold on
histfit(-real(d1sim(:,2)),20)
axis square
xlim([-1e-3 3e-3])
ylim([0 1000])
title('Fig. 2D')

%now generate distribution of growth rate values from growth measurements

vMaxNaiveMean=mean(vMaxNaive);
vMaxNaiveStd=std(vMaxNaive);

vMaxBIGMean=mean(vMaxBIG);
vMaxBIGStd=std(vMaxBIG);

for i=1:M
    
    naiveToFit=vMaxNaiveMean+randn*vMaxNaiveStd;
    BIGtoFit=vMaxBIGMean+randn*vMaxBIGStd;
    
    u0sim(i)=mean(m(1:3))*naiveToFit/vMaxNaiveMean*60;
    u1sim(i)=mean(m(1:3))*BIGtoFit/vMaxNaiveMean*60;
end

%plot distributions
subplot(2,3,2)
histfit(u0sim,25)
hold on
histfit(u1sim,20)
axis square
xlim([0.08 0.16])
ylim([0 160])
title('Fig. 2E')

%generate cloud of projected competition results based on these sims

for i=1:M
%initialize
    params.x00=1;
    params.x10=1;
    params.u0=u0sim(i);
    params.u1=u1sim(i);
    params.d0=d0sim(i,2);
    params.d1=d1sim(i,2);


    %example traces
    [t(:,i),BIGfracSim(:,i)]=oscillate(params);
end

BIGfracSim=real(BIGfracSim);

%plot confidence intervals on model prediction
upper=zeros(params.N,1);
lower=upper;
x=0.95;

for i=1:params.N %for each time point
    
    sorted=sort(BIGfracSim(i,:));
    
    upper(i)=sorted(x*length(sorted));
    lower(i)=sorted(ceil((1-x)*length(sorted)));
    
end

%plot error cloud 
subplot(2,3,3)
vGens=(0:(params.N-1))*6.5;
fill([vGens fliplr(vGens)],log10([upper' fliplr(lower')]),'red')
fill([vGens fliplr(vGens)],log10([1./upper' fliplr(1./lower')]),'blue')
alpha 0.3
xlim([0 125])



set(gcf,'PaperPositionMode','auto')
print(['modelFigs' num2str(figureCounter)],'-dsvg','-r0')
print(['modelFigs' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;



%initialize
params.x00=1;
params.x10=1;
params.u0=1;
params.u1=1.01;
params.d0=1;
params.d1=1.01;

params.t1=1;
params.t2=1;

params.N=25;

params0=params;

%example traces
[t,BIGfrac]=oscillate(params);



figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
%plot
plot(t,log10(BIGfrac),'r')
hold on
plot(t,log10(1./BIGfrac),'b')
xlabel('time')
ylabel('log_{10} BIG+/naive')
title('Fig. S2C')
axis square
ylim([-2 2])

params=params0;
params.t1=1;
params.t2=10;

%example traces
[t,BIGfrac]=oscillate(params);
subplot(2,3,2)
%plot
plot(t,log10(BIGfrac),'r')
hold on
plot(t,log10(1./BIGfrac),'b')
xlabel('time')
ylabel('log_{10} BIG+/naive')
title('Fig. S2A')
axis square
ylim([-2 2])

params=params0;
params.t1=10;
params.t2=1;
%example traces
[t,BIGfrac]=oscillate(params);
subplot(2,3,3)
%plot
plot(t,log10(BIGfrac),'r')
hold on
plot(t,log10(1./BIGfrac),'b')
xlabel('time')
ylabel('log_{10} BIG+/naive')
title('Fig. S2B')
axis square
ylim([-2 2])


params=params0;
params.u1=1.01;
params.d1=1.01;

%sweep params
for i=1:20
    for j=1:20
        params.t1=i/2;
        params.t2=j/2;
        [t,BIGfrac]=oscillate(params);
        ratio(21-i,j)=BIGfrac(end);
    end
end


subplot(2,3,4)
%plot
imagesc(log10(ratio))
xlabel('starvation time')
ylabel('growth time')
colorbar
title('Fig. 2C')
axis square

colormap jet


set(gcf,'PaperPositionMode','auto')
print(['modelFigs' num2str(figureCounter)],'-dsvg','-r0')
print(['modelFigs' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;


%parameter sweeps (2D)

%initialize
params.x00=1;
params.x10=1;
params.u0=1;
params.u1=1.01;
params.d0=1;
params.d1=1.01;

params.t1=1;
params.t2=1;

params.N=25;

params0=params;

paramsToSweep={'u0','u1','d0','d1','t1','t2'};


figure('units','normalized','outerposition',[0 0 1 1])
m=1;
for i=1:length(paramsToSweep)
    
    v1=getfield(params0,paramsToSweep{i})*10.^(-1:0.02:1);
    
    for j=(i+1):length(paramsToSweep)
        
        v2=getfield(params0,paramsToSweep{j})*10.^(-1:0.02:1);
        
        params=params0;
        
        for k=1:length(v1)
            
            for l=1:length(v2)
                
                params=setfield(params,paramsToSweep{i},v1(k));
                params=setfield(params,paramsToSweep{j},v2(l));
                
                [t,BIGfrac]=oscillate(params);
                
                tempFrac(k,l)=BIGfrac(end);
                
            end
            
        end
        
        subplot(4,4,m)
        imagesc(log10(flipud(tempFrac)))
        xlabel(paramsToSweep{j})
        ylabel(paramsToSweep{i})
        axis square
        colorbar
        
        xticks(1:20:length(v1))
        xticklabels(round(v1(1:20:length(v1)),1))
        xtickangle(45)
        
        yticks(1:20:length(v2))
        yticklabels(round(fliplr(v1(1:20:length(v2))),1))
        
        m=m+1;
                
    end
        
    
end





set(gcf,'PaperPositionMode','auto')
print(['modelFigs' num2str(figureCounter)],'-dsvg','-r0')
print(['modelFigs' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;



