clc;
clear all; %close all;
customColorPallet = loadCustomColorPallets;

file = 'E:\logs\IM-1017\2017-08-17\IM-1017_2017-08-17_10-04-41____Ali_Conditioning';

[file,path,FilterIndex]  = uigetfile('*.*');
file = [path filesep file];

[pathstr,nameSession,ext] = fileparts(file);
file = [pathstr nameSession];


software=1;  % Software lockin
figureExport=1;
%saveDirectory='C:\Users\alimo\GDrive\Berke Lab\workflow\Fig graveyard\photometry\Pav';
saveDirectory='/Users/ali/Google Drive/Berke Lab/workflow/Fig graveyard/photometry/Pav';

%% Extract Data
phot=readPhotometryData(file);
[box boxts ratConditioning] = readBoxData_Pav(file);

tau = 10;
filterOrder = 5;


detector = phot.Data(3,:)/2^15*10;
exc1 = phot.Data(1,:)/2^15*10;
exc2 = phot.Data(2,:)/2^15*10;
[sig1,ref1]=lockinDetection(detector,exc1,exc2,phot.SamplingRate,'tau',tau,'filterorder',filterOrder,'detrend',false,'Full',true);

detector = phot.Data(7,:)/2^15*10;
exc1 = phot.Data(3,:)/2^15*10;
exc2 = phot.Data(1,:)/2^15*10;
[sig2,ref2]=lockinDetection(detector,exc1,exc2,phot.SamplingRate,'tau',tau,'filterorder',filterOrder,'detrend',false,'Full',true);

detector = phot.Data(6,:)/2^15*10;
exc1 = phot.Data(7,:)/2^15*10;
exc2 = phot.Data(8,:)/2^15*10;
[sig3,ref3]=lockinDetection(detector,exc1,exc2,phot.SamplingRate,'tau',tau,'filterorder',filterOrder,'detrend',false,'Full',true);

detector = phot.Data(5,:)/2^15*10;
exc1 = phot.Data(4,:)/2^15*10;
exc2 = phot.Data(8,:)/2^15*10;
[sig4,ref4]=lockinDetection(detector,exc1,exc2,phot.SamplingRate,'tau',tau,'filterorder',filterOrder,'detrend',false,'Full',true);


%% Downsample and smooth
Fs = 250; 
bw=25;

sig1 = downsample(sig1,phot.SamplingRate/Fs);
ref1 = downsample(ref1,phot.SamplingRate/Fs);
sig1 = filtfilt(ones(1,bw)./bw,1,sig1);
ref1 = filtfilt(ones(1,bw)./bw,1,ref1);

sig2 = downsample(sig2,phot.SamplingRate/Fs);
ref2 = downsample(ref2,phot.SamplingRate/Fs);
sig2 = filtfilt(ones(1,bw)./bw,1,sig2);
ref2 = filtfilt(ones(1,bw)./bw,1,ref2);

sig3 = downsample(sig3,phot.SamplingRate/Fs);
sig3 = filtfilt(ones(1,bw)./bw,1,sig3);
ref3 = downsample(ref3,phot.SamplingRate/Fs);
ref3 = filtfilt(ones(1,bw)./bw,1,ref3);

sig4 = downsample(sig4,phot.SamplingRate/Fs);
sig4 = filtfilt(ones(1,bw)./bw,1,sig4);
ref4 = downsample(ref4,phot.SamplingRate/Fs);
ref4 = filtfilt(ones(1,bw)./bw,1,ref4);



boxts.event = round(Fs*boxts.event / (phot.SamplingRate));
% 
FoodLine = downsample(boxts.FoodLine,phot.SamplingRate/Fs);
% 
Location1 = [phot.Channel(1).Location(1:20)'];
Location2 = [phot.Channel(2).Location(1:20)'];
% 




%%
    
clear dF1 dF2 dF3 dF4 foodEntry

tWin =[-5 45];

% 
[dF_F1,ref_fitted1,slope] = isosbestic_correction(sig1,ref1,'type','linear');
[dF_F2,ref_fitted2,slope] = isosbestic_correction(sig2,ref2,'type','linear');
[dF_F3,ref_fitted3,slope] = isosbestic_correction(sig3,ref3,'type','linear');
[dF_F4,ref_fitted4,slope] = isosbestic_correction(sig4,ref4,'type','linear');

% 

% 

for trial=1:size(boxts.event,2)
    event=1;baselineWin=[-5 -2];
    trialAvg_sig1 = mean(sig1(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
    trialAvg_sig2 = mean(sig2(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
    trialAvg_sig3 = mean(sig3(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
    trialAvg_sig4 = mean(sig4(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
    trialStd_sig1 = std(sig1(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
    trialStd_sig2 = std(sig2(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
    trialStd_sig3 = std(sig3(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
    trialStd_sig4 = std(sig4(boxts.event(event,trial)+baselineWin(1)*Fs:boxts.event(event,trial)+baselineWin(2)*Fs));
  
    for event=1:size(boxts.event,1)
        
        if ~isnan(boxts.event(event,trial)) && (boxts.event(event,trial)+tWin(1)*Fs>0) && (boxts.event(event,trial)+tWin(2)*Fs<length(sig1))
%             dF1(trial,:,event) = (dF_F1(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs));
%             dF2(trial,:,event) = (dF_F2(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs));
%             dF3(trial,:,event) = (dF_F3(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs));
%             dF4(trial,:,event) = (dF_F4(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs));

            dF1(trial,:,event) = (sig1(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs)-trialAvg_sig1)/trialStd_sig1;
            dF2(trial,:,event) = (sig2(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs)-trialAvg_sig2)/trialStd_sig2;
            dF3(trial,:,event) = (sig3(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs)-trialAvg_sig3)/trialStd_sig3;
            dF4(trial,:,event) = (sig4(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs)-trialAvg_sig4)/trialStd_sig4;

%             
            foodEntry(trial,:,event) = FoodLine(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs);
            
        else
            dF1(trial,:,event) = nan(1,1+Fs*(tWin(2)-tWin(1)));
            dF2(trial,:,event) = nan(1,1+Fs*(tWin(2)-tWin(1)));
            dF3(trial,:,event) = nan(1,1+Fs*(tWin(2)-tWin(1)));
            dF4(trial,:,event) = nan(1,1+Fs*(tWin(2)-tWin(1)));
                        foodEntry(trial,:,event) = nan(1,1+Fs*(tWin(2)-tWin(1)));
        end
    end
end



%% Response breakdown by Rew

h=figure('units','normalized','outerposition',[0 0 1 1]); %[left bottom width height]

xLimWin=[-2 33];
data = dF1;
event = 1;
clear hh hhh

ind = -1*ones(size(ratConditioning.Tone));
ind(ratConditioning.Tone==0) = 1;
ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==75)) = 2;
ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==25)) = 3;
ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==0)) = 4;

rew = ratConditioning.Food;

group = 1;
for cond=1:4
    x = foodEntry(:,:,event);
    x_axis = linspace(tWin(1),tWin(2),1+(tWin(2)-tWin(1))*Fs);
    hh(cond)=subplot('Position', [0.05+(cond-1)*.13 0.88-(2*group-2)*.13 0.11 0.09]); hold on %[left bottom width height]
    imagesc(x_axis,1:sum(ind==cond),1-x(ind==cond,:)>0);colormap(gca,gray)
    xlim(xLimWin);set(gca,'XTICK','','YTICK','');%ylim([1 sum(ind==cond)])
    if cond==2 && group==1,title(nameSession);end
    set(gca,'FontSize',16);grid on
    
    hhh(cond)= subplot('Position', [0.05+(cond-1)*.13 0.88-(2*group-1)*.13 0.11 0.13]); hold on %[left bottom width height]
    bar(x_axis,100*nansum(x(ind==cond,:))/sum(ind==cond),1,'EdgeColor','k','FaceColor','k')
    xlim(xLimWin);
    set(gca,'FontSize',16)
    if cond==1, ylabel('P(Food Entry)');end
    set(gca,'XTICK',[-3 0 3])
end
% linkprop(hh,'CLim')
linkaxes(hhh,'y')



conditionName = {'Unpredicted','75%','25%','0%'};
heatMapColor = customColorPallet.default;

for group=2:5    %         rewarded = find(ratConditioning.Food==1);
    
    for cond=1:4
        
        
        
        switch group
            case 2
                x = dF1(:,:,event);
                offset1=0.86; offset2=0.88;scale=.1;
                signalType='Green Indicator';
            case 3
                x = dF2(:,:,event);
                offset1=0.92; offset2=0.94;scale=.1;
                signalType='Red Indicator';
            case 4
                x = dF3(:,:,event);
                offset1=0.89; offset2=0.91;scale=.1;
                signalType='Green Indicator';
            case 5
                x = dF4(:,:,event);
                offset1=0.95; offset2=0.97;scale=.1;
                signalType='Red Indicator';
        end
        
        
        x_axis = linspace(tWin(1),tWin(2),1+(tWin(2)-tWin(1))*Fs);
        
        
        hh(cond+(group-2)*4)=subplot('Position', [0.05+(cond-1)*.13 offset1-(2*group-2)*scale 0.11 0.06]); hold on %[left bottom width height]
        
        xx = [x(ratConditioning.Food==1 & ind==cond,:);x(ratConditioning.Food==0 & ind==cond,:)];
        for iii=1:size(xx,1)
            xx(iii,:)=medfilt1(xx(iii,:),21);
        end
%         if event==1
%             imagesc(x_axis,1:size(xx,1),xx,[.1*min(mean(xx,1)) .9*max(mean(xx,1))]);hold on;%
%         else
            imagesc(x_axis,1:size(xx,1),xx,[-2 4]);hold on;%
%         end

        set(gcf,'colormap', heatMapColor);shading interp;%view(0,70)
        plot(x_axis,sum(ratConditioning.Food==1 & ind==cond)*ones(size(x_axis)),'linestyle','--','color','w','linewidth',2);hold on
        xlim(xLimWin);set(gca,'XTICK','','YTICK','');ylim([1 size(xx,1)])
        set(gca,'FontSize',16)
        
        hhh(cond+(group-2)*4)= subplot('Position', [0.05+(cond-1)*.13 offset2-(2*group-1)*scale 0.11 0.08]); hold on %[left bottom width height]
        h1=plot(x_axis,medfilt1(nanmean(x(ratConditioning.Food==1 & ind==cond,:)),3),'color','r','linewidth',2);hold on
        h2=plot(x_axis,medfilt1(nanmean(x(ratConditioning.Food==0 & ind==cond,:)),3),'color','b','linewidth',2);
        grid on
        xlabel(conditionName{cond})
        xlim(xLimWin)
        set(gca,'FontSize',16)
        set(gca,'XTICK',[-3 0 3])
        if cond==1, ylabel(signalType);end
        if cond==4,legend([h1,h2],{'Reward','Omission'},'location','northeast');    legend boxoff;end
    end
end
linkprop(hh(1:4),'CLim');linkprop(hh(5:8),'CLim');linkprop(hh(9:12),'CLim');linkprop(hh(13:16),'CLim');

linkaxes(hhh(1:4),'y');linkaxes(hhh(5:8),'y');linkaxes(hhh(9:12),'y');linkaxes(hhh(13:16),'y');


t0=Fs*60:Fs*120;

t_axis=linspace(0,length(sig1)/Fs,length(sig1));
subplot('Position', [0.65 0.55 0.2 0.25]); hold on %[left bottom width height]
m=mean(sig1(t0));h1=plot(t_axis/60,100*(sig1-m)/m,'color','b');xlim([1 t_axis(end)/60])
m=mean(ref1(t0));h2=plot(t_axis/60,100*(ref1-m)/m,'color',[130,0, 200]/255);xlim([1 t_axis(end)/60])
m=mean(sig2(t0));h3=plot(t_axis/60,100*(sig2-m)/m,'color',[255,0, 0]/255);xlim([1 t_axis(end)/60])
ylabel('Percent change in F','FontSize',16)
xlabel('Time(min)','FontSize',16)
legend([h1 h2 h3],{'470nm','405nm','565nm'});legend boxoff
set(gca,'FontSize',16)
title(Location1)

t_axis=linspace(0,length(sig1)/Fs,length(sig1));
subplot('Position', [0.65 0.1 0.2 0.25]); hold on %[left bottom width height]
m=mean(sig3(t0));h1=plot(t_axis/60,100*(sig3-m)/m,'color','b');xlim([1 t_axis(end)/60])
m=mean(ref3(t0));h2=plot(t_axis/60,100*(ref3-m)/m,'color',[130,0, 200]/255);xlim([1 t_axis(end)/60])
m=mean(sig4(t0));h3=plot(t_axis/60,100*(sig4-m)/m,'color',[255,0, 0]/255);xlim([1 t_axis(end)/60])
ylabel('Percent change in F','FontSize',16)
xlabel('Time(min)','FontSize',16)
legend([h1 h2 h3],{'470nm','405nm','565nm'});legend boxoff
set(gca,'FontSize',16)
title(Location2)



set(h,'color','w')

if figureExport
set(h,'Units','inches');
screenposition = get(h,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print('-dpdf','-painters',[saveDirectory filesep nameSession '_'  datestr(now, 'yyyy-mm-dd') '_test'])

end

