%function [] = SingleVarAnalysis_BasedOnCoVarData(MinimumSignal2NoiseStrength,alpha,CutValue,MinimumNumberEstimatableKDs,Experiments,GenerateOutput)
% 
% close all
% clear all

%       Library type	Sample Type	Population	[Gag] concentration
% 1     WT          RNA	Beads                   2000
% 2     WT          RNA	Beads                   200
% 3     WT          RNA	Beads                   20
% 4     Low         RNA	Beads                   2000
% 5     Low         RNA	Beads                   200
% 6     Low         RNA	Beads                   20
% 7     High        RNA	Beads                   2000
% 8     High        RNA	Beads                   200
% 9     High        RNA	Beads                   20
% 10	WT          RNA	Supernatant             2000
% 11	WT          RNA	Supernatant             200
% 12	WT          RNA	Supernatant             20
% 13	Low         RNA	Supernatant             2000
% 14	Low         RNA	Supernatant             200
% 15	Low         RNA	Supernatant             20
% 16	High        RNA	Supernatant             2000
% 17	High        RNA	Supernatant             200
% 18	High        RNA	Supernatant             20
% 19	WT          DNA	Plasmid	
% 20	Low         DNA	Plasmid	
% 21	High        DNA	Plasmid	
% 

close all
clear all
%with Watson-Crick and Wobble base- vs. non Watson-Crick base pairing
%           AA   AC   AG   AT   CA   CC   CG   CT   GA   GC   GG   GT   TA
%           TC   TG   TT

folder = '';

M = dlmread('../RefSeq.txt');
Positions = 199:386;%M(:,1)';

%Experiment = 1;
for counter = 1:38
    Experiment = counter;
    disp(strcat('---Experiment: ',num2str(Experiment),' ---'))
    
    for j = Positions
        disp(strcat('Position: ',num2str(j),' ..'))
        filename1 = strcat('./',folder,num2str(Experiment),'/',num2str(Experiment),'_',num2str(j),'.txt');
        M = dlmread(filename1);
        delete(filename1);
        
        %clear lines where pos 2 <= pos 1
        I = M(:,3)>M(:,1);
        Mtmp = M(I,:);
        dlmwrite(filename1,Mtmp);
    end
end
