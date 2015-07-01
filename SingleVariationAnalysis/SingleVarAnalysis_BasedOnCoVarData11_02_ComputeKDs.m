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
% MinimumSignal2NoiseStrength = floor(MinimumSignal2NoiseStrength); %affects only supernatant counts
% CutValue = floor(CutValue);
% MinimumNumberEstimatableKDs = floor(MinimumNumberEstimatableKDs);
% GenerateOutput = floor(GenerateOutput);
% Experiments = floor(Experiments);
% close all;

%----------------------
close all
clear all
%MinimumSignal2NoiseStrength = 3; %affects only supernatant counts
%CutValue = 30;
CutValueFwd = 30;
CutValueBwd = 45;
ExperimentalPairs = [[4:9]',[4:9]'+9]; 
%WeightThreshold = 0.4;%of maximum read coverage [fraction of maximum] to be included in analysis
%--------------------------

NrPositionsTotal = 535;
bases = 1:4;
%For re-ordering entries in the case where the first-and second base were
%swapt
%-----------------------------

%Read Reference sequence
disp('Read RefSeq..')
Tmp = dlmread('./Data_CoVariation/RefSeq.txt');
RefSeq = Tmp(:,2);

%----------------------------

%load 'Results.mat';

NrComparisons = length(ExperimentalPairs);

%load Error estimates
% gives error matrices type of transition x position
load './ErrorEstimates.mat';
% MedianExpKappa_w2m1_Total
% MeanExpKappa_w2m1_Total 

%individual mutation probabilities, e.g. A -> C
MedianExpKappa_w2m1_Total(MedianExpKappa_w2m1_Total==0) = nan;
MeanExpKappa_w2m1_Total(MeanExpKappa_w2m1_Total==0) = nan;
%all mutations away from wt, e.g. A -> {C,G,T}
MedianExpKappa_w2M1_Total(MedianExpKappa_w2M1_Total==0) = nan;
MeanExpKappa_w2M1_Total(MeanExpKappa_w2M1_Total==0) = nan;

TransitionRule =  zeros(1,20);
TransitionRule(1:4) = [3:4,1:2];
for i = 1:4
    for j = 1:4
        TransitionRule((i-1)*4+j+4) = (j-1)*4+i+4;
    end
end

TotalRel_KD_m_w = nan(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4);
TotalRel_KD_w_m = nan(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4);

% LowerLimitsKD_w_m = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-2*CutValue),4); 
% LowerLimitsKD_m_w = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-2*CutValue),4); 

Signal2Noise_m_w_SuperN = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4);
Signal2Noise_w_m_SuperN = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4);
Signal2Noise_m_w_Beads = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4);
Signal2Noise_w_m_Beads = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4);

TotalRel_KD_M_w = nan(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)));
TotalRel_KD_w_M = nan(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)));

% LowerLimitsKD_w_M = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-2*CutValue)); 
% LowerLimitsKD_M_w = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-2*CutValue)); 

Signal2Noise_M_w_SuperN = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)));
Signal2Noise_w_M_SuperN = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)));
Signal2Noise_M_w_Beads = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)));
Signal2Noise_w_M_Beads = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)));

%
TotalNrSeqBeads = zeros(NrPositionsTotal,NrComparisons);
TotalNrSeqSuperN = zeros(NrPositionsTotal,NrComparisons);
TotalNrSeq = zeros(NrPositionsTotal,NrComparisons);
%
PositionWeightsBeads = zeros(NrPositionsTotal,NrPositionsTotal,NrComparisons); 
PositionWeightsSuperN = zeros(NrPositionsTotal,NrPositionsTotal,NrComparisons); 
PositionWeightsTotal = zeros(NrPositionsTotal,NrPositionsTotal,NrComparisons); 

PositionReadsBeads = zeros(NrPositionsTotal,NrPositionsTotal,NrComparisons); 
PositionReadsSuperN = zeros(NrPositionsTotal,NrPositionsTotal,NrComparisons); 
PositionReadsTotal = zeros(NrPositionsTotal,NrPositionsTotal,NrComparisons); 

% go through positions r_1
for outercounter = CutValueFwd+1:NrPositionsTotal-CutValueBwd %250:250%<----------------------------------------------------------
    
    Position = outercounter
    %disp(num2str(Position));
    
    TotalM1 = zeros(NrPositionsTotal,20,NrComparisons); % take together all experiments [Gag Concentrations]
    TotalM2 = TotalM1; % take together all experiments [Gag Concentrations]
    
    for counter = 1:NrComparisons
        ExperimentalPair = ExperimentalPairs(counter,:);

        filename1 = strcat('./Data_CoVariation/',num2str(ExperimentalPair(1)),'/',num2str(ExperimentalPair(1)),'_',num2str(Position),'.txt');
        M1 = dlmread(filename1);

        filename2 = strcat('./Data_CoVariation/',num2str(ExperimentalPair(2)),'/',num2str(ExperimentalPair(2)),'_',num2str(Position),'.txt');
        M2 = dlmread(filename2);   

        %Matrix M: 
        % pos1 Nt@pos1 pos2 Nt@pos2 AA AC AG AT CA .... TT

        %Normalization
        [~,Phi,Theta] = Normalizationfactors();

        M1(:,5:end) = M1(:,5:end)/Phi(ExperimentalPair(1));
        M2(:,5:end) = M2(:,5:end)/Theta(ExperimentalPair(1));

        [NrEntries,~] = size(M1);
        
        %---------
        %Sort the Matrices
        %re-order if necessary, i.e if not starting with 'Position'
        for i = 1:NrEntries
            if M1(i,3) == Position
                M1(i,:) = M1(i,TransitionRule);
            end
        end
        [NrEntries,NrCollums] = size(M2);
        for i = 1:NrEntries
            if M2(i,3) == Position
                M2(i,:) = M2(i,TransitionRule);
            end
        end
        
        %Sort in ascending Order
        [~,IDX1] = sort(M1(:,3));
        M1_tmp = M1(IDX1,:);
        [~,IDX2] = sort(M2(:,3));
        M2_tmp = M2(IDX2,:);
        
        %---------
        
        M1 = nan(NrPositionsTotal,20);

        M1(:,1) = Position;
        M1(:,2) = RefSeq(Position);
        M1(:,3) = 1:NrPositionsTotal;
        M1(1:length(RefSeq),4) = RefSeq(:);
        if length(RefSeq) < NrPositionsTotal
            M1(length(RefSeq)+1:NrPositionsTotal,4) = 1;
        end

        M2 = M1;

        [~,~,ib] = intersect(M1_tmp(:,3),1:NrPositionsTotal);
        M1(ib,5:end) = M1_tmp(:,5:end);

        [~,~,ib] = intersect(M2_tmp(:,3),1:NrPositionsTotal);
        M2(ib,5:end) = M2_tmp(:,5:end);
        
        TotalM1(:,1:4,counter) = M1(:,1:4);
        TotalM1(:,5:end,counter) = M1(:,5:end);
        TotalM2(:,1:4,counter) = M2(:,1:4);
        TotalM2(:,5:end,counter) = M2(:,5:end);
        
    end % end over experiments
     
    % go through positions r_2
    RunningIndex = 0;
    for i = CutValueFwd+1:NrPositionsTotal-CutValueBwd%NrEntries   
        
        %wt positions
        pos1 = TotalM1(i,2,1);
        pos2 = TotalM1(i,4,1);
        
        idx1 = 1:4;% bases(1:4~=pos1);
        idx2 = 1:4;%bases(1:4~=pos2);
        
        pos_w_w = 4*(pos1-1)+pos2+4;
        pos_m_w = 4*(idx1-1)+pos2+4;
        pos_w_m = 4*(pos1-1)+idx2+4;
        for counter = 1:NrComparisons 
            
            RunningIndex = RunningIndex+1;
            %------------------Censoring-----------------------------------
            NrSeqExperimentPosBeads = nansum(TotalM1(i,5:end,counter));
            NrSeqExperimentPosSuperN = nansum(TotalM2(i,5:end,counter));
            
            TotalNrSeqBeads(i,counter) = NrSeqExperimentPosBeads;
            TotalNrSeqSuperN(i,counter) = NrSeqExperimentPosSuperN;
            TotalNrSeq(i,counter) = NrSeqExperimentPosBeads + NrSeqExperimentPosSuperN;
            %------------------Censoring-----------------------------------
            
            %wt_wt numbers
            Beads_w_w = TotalM1(i,pos_w_w,counter);
            SuperN_w_w = TotalM2(i,pos_w_w,counter);
           
            %-----------------  m_wt numbers  -----------------------
            Beads_m_w = TotalM1(i,pos_m_w,counter);
            SuperN_m_w = TotalM2(i,pos_m_w,counter);
            
            Nominator = SuperN_m_w/SuperN_w_w - MedianExpKappa_w2m1_Total(4*(pos1-1)+idx1,outercounter)';%superN
            Denominator = Beads_m_w/Beads_w_w - MedianExpKappa_w2m1_Total(4*(pos1-1)+idx1,outercounter)';%beads

            Rel_KD_m_w = zeros(1,length(idx1));
            %Signal-2-Noise
            for j = idx1
                Noise = MedianExpKappa_w2m1_Total(4*(pos1-1)+j,outercounter);
                
                Signal2Noise_m_w_SuperN(outercounter,RunningIndex,j) = SuperN_m_w(j)/(SuperN_w_w*Noise);
                Signal2Noise_m_w_Beads(outercounter,RunningIndex,j) = Beads_m_w(j)/(Beads_w_w*Noise);
                
                if SuperN_m_w(j)/(SuperN_w_w*Noise) == inf%% not enough signal at supernatent => signal-to-noise-ratio too bad
                    Rel_KD_m_w(j) = nan; 
                elseif Denominator(j) <= 0 %Denominator <= 0 && Nominator > 0 %not enough signal at beads
                    Rel_KD_m_w(j)  = -1;%Nominator/1e-8; <-------------------------------------------arbitrary
                    %LowerLimitsKD_m_w(outercounter,RunningIndex,j) = 1; 
                elseif Nominator(j) <= 0
                    Rel_KD_m_w(j)  = -10;%Nominator/1e-8; <-------------------------------------------arbitrary
                else
                    Rel_KD_m_w(j) = Nominator(j)/Denominator(j); 
                end
            end

            TotalRel_KD_m_w(outercounter,RunningIndex,:) = Rel_KD_m_w;

            %--------------  M -> wt numbers  -------------------%
            pos_M_w = 4*(bases(1:4~=pos1)-1)+pos2+4;
            Beads_M_w = nansum(TotalM1(i,pos_M_w,counter),2);
            SuperN_M_w = nansum(TotalM2(i,pos_M_w,counter),2);
             
            Nominator = SuperN_M_w/SuperN_w_w - MedianExpKappa_w2M1_Total(outercounter);%superN
            Denominator = Beads_M_w/Beads_w_w - MedianExpKappa_w2M1_Total(outercounter);%beads
            
            Noise = MedianExpKappa_w2M1_Total(outercounter);
            
            Signal2Noise_M_w_SuperN(outercounter,RunningIndex) = SuperN_M_w/(SuperN_w_w*Noise);
            Signal2Noise_M_w_Beads(outercounter,RunningIndex) = Beads_M_w/(Beads_w_w*Noise);
            
            if SuperN_M_w/(SuperN_w_w*Noise) == inf%% not enough signal at supernatent => signal-to-noise-ratio too bad
                Rel_KD_M_w = nan; 
            elseif Denominator <= 0  %Denominator <= 0 && Nominator > 0 %not enough signal at beads
                Rel_KD_M_w  = -1;%Nominator/1e-8; <-------------------------------------------arbitrary
                %LowerLimitsKD_M_w(outercounter,RunningIndex) = 1; 
            elseif Nominator <= 0
                Rel_KD_M_w  = -10;%Nominator/1e-8; <-------------------------------------------arbitrary
            else
                Rel_KD_M_w = Nominator/Denominator; 
            end
            TotalRel_KD_M_w(outercounter,RunningIndex) = Rel_KD_M_w;
            %-----------------------------------------------%
            
            %----------------  wt_m numbers  -------------------%
            Beads_w_m = TotalM1(i,pos_w_m,counter);
            SuperN_w_m = TotalM2(i,pos_w_m,counter); 
           
            Nominator = SuperN_w_m/SuperN_w_w - MedianExpKappa_w2m1_Total(4*(pos2-1)+idx2,i)';%superN
            Denominator = Beads_w_m/Beads_w_w - MedianExpKappa_w2m1_Total(4*(pos2-1)+idx2,i)';%beads
        
            Rel_KD_w_m = zeros(1,length(idx2));
            %Signal-2-Noise
            for j = idx2
                Noise = MedianExpKappa_w2m1_Total(4*(pos2-1)+j,i);
                
                Signal2Noise_w_m_SuperN(outercounter,RunningIndex,j) = SuperN_w_m(j)/(SuperN_w_w*Noise);
                Signal2Noise_w_m_Beads(outercounter,RunningIndex,j) = Beads_w_m(j)/(Beads_w_w*Noise);
                
                if SuperN_w_m(j)/(SuperN_w_w*Noise) == inf%% not enough signal at supernatent => signal-to-noise-ratio too bad
                    Rel_KD_w_m(j) = nan; 
                elseif Denominator(j) <= 0
                    Rel_KD_w_m(j)  = -1;%Nominator/1e-8; <-------------------------------------------arbitrary
                    %LowerLimitsKD_w_m(outercounter,RunningIndex,j) = 1; 
                elseif Nominator(j) <= 0
                    Rel_KD_w_m(j)  = -10;%Nominator/1e-8; <-------------------------------------------arbitrary
                else
                    Rel_KD_w_m(j) = Nominator(j)/Denominator(j); 
                end
            end
            TotalRel_KD_w_m(outercounter,RunningIndex,:) = Rel_KD_w_m;
            
            
            %--------------  wt-M numbers  --------------------%
            pos_w_M = 4*(pos1-1)+bases(1:4~=pos2)+4;  
            Beads_w_M = nansum(TotalM1(i,pos_w_M,counter),2);
            SuperN_w_M = nansum(TotalM2(i,pos_w_M,counter),2);
             
            Nominator = SuperN_w_M/SuperN_w_w - MedianExpKappa_w2M1_Total(i);%superN
            Denominator = Beads_w_M/Beads_w_w - MedianExpKappa_w2M1_Total(i);%beads
            
            Noise = MedianExpKappa_w2M1_Total(i);
            
            Signal2Noise_w_M_SuperN(outercounter,RunningIndex) = SuperN_w_M/(SuperN_w_w*Noise);
            Signal2Noise_w_M_Beads(outercounter,RunningIndex) = Beads_w_M/(Beads_w_w*Noise);
            
            if SuperN_w_M/(SuperN_w_w*Noise) == inf%% not enough signal at supernatent => signal-to-noise-ratio too bad
                Rel_KD_w_M = nan; 
            elseif Denominator <= 0  %Denominator <= 0 && Nominator > 0 %not enough signal at beads
                Rel_KD_w_M  = -1;%Nominator/1e-8; <-------------------------------------------arbitrary
                %LowerLimitsKD_w_M(outercounter,RunningIndex) = 1; 
            elseif Nominator <= 0
                    Rel_KD_w_M  = -10;%Nominator/1e-8; <-------------------------------------------arbitrary
            else
                Rel_KD_w_M = Nominator/Denominator; 
            end
            TotalRel_KD_w_M(outercounter,RunningIndex) = Rel_KD_w_M;
            %------------
            
            
        end%counter
       
    end%end i
    
    for counter = 1:NrComparisons 
    % censoring by Number of Reads
        PositionWeightsBeads(outercounter,:,counter) = (TotalNrSeqBeads(:,counter)'./max(TotalNrSeqBeads(:,counter)));
        PositionWeightsSuperN(outercounter,:,counter) = (TotalNrSeqSuperN(:,counter)'./max(TotalNrSeqSuperN(:,counter)));
        PositionWeightsTotal(outercounter,:,counter) = (TotalNrSeq(:,counter)'./max(TotalNrSeq(:,counter)));
        
        PositionReadsBeads(outercounter,:,counter) = TotalNrSeqBeads(:,counter)';
        PositionReadsSuperN(outercounter,:,counter) = TotalNrSeqSuperN(:,counter)'; 
        PositionReadsTotal(outercounter,:,counter) = TotalNrSeq(:,counter)'; 
    end
    %
  
end%end outercounter

%filename = strcat('./ResultsSN',num2str(MinimumSignal2NoiseStrength),'W',num2str(10*WeightThreshold),'.mat');
filename = './Results.mat';

try delete(filename)
    pause(1);
catch
end
save(filename, 'CutValueFwd','CutValueBwd','NrPositionsTotal','NrComparisons','TotalRel_KD_m_w','TotalRel_KD_M_w',...
    'PositionWeightsBeads','PositionWeightsSuperN','PositionWeightsTotal','Signal2Noise_m_w_SuperN','Signal2Noise_M_w_SuperN',...
    'Signal2Noise_m_w_Beads','Signal2Noise_M_w_Beads','PositionReadsBeads','PositionReadsSuperN','PositionReadsTotal');
