
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

% M = dlmread('InputVariables2D.txt');
% MinimumSignal2NoiseStrength_m_m = M(1);
% MinimumSignal2NoiseStrength_m = M(2);
% alpha =  M(3);
% MinimumNumberEstimatableKDs =  M(4);
% WeightThreshold =  M(5);
% WriteOutput =  M(6);
% ArbitraryValueLowerEstimate =  M(7);
% ArbitraryValueUpperEstimate =  M(8);
% MinimumStemSize =  M(9);
% penalty =  M(10);
% smooth =  M(11);
% Fix253_265 = M(12);
% %----------------------------------------
close all
clear all
MinimumSignal2NoiseStrength_m_m = 6;
MinimumSignal2NoiseStrength_m = 0;
alpha = 0.99
MinimumNumberEstimatableKDs = 20;
WeightThreshold = 0.2;
%ReadThreshold = 5000;
WriteOutput = 1;
ArbitraryValueLowerEstimate = 0;%100;%50;
ArbitraryValueUpperEstimate = 0;%1/100;%1/50;
MinimumStemSize = 3;
penalty = -1; %very sensitive parameter. set to -inf if only connected stems (where there are consecutive Epistasis estimates) should be identified
smooth = 0;
%Fix253_265 = 0;

%----------------------------------------
DisStart = 243;
DisEnd = 277;
DisRegion = DisStart:DisEnd;
FirstPos = 200;
LastPos = 349;
NrPositionsTotal = LastPos-FirstPos;
NrComparisons = 6;
%----------------------------------------


Infolder = './ResultsTriVarWCWB_new/';

WCBasePairs = [4 7 10 13];% %AU, CG, GC, UA
WCSingleBases = [1 4; 2 3; 3 2; 4 1];
WCWBBasePairs = [4 7 10 12 13 15];% %AU, CG, GC, GU, UA, UG
WCWBSingleBases = [1 4; 2 3; 3 2; 3 4; 4 1; 4 3];

Bases = 1:4;

M = dlmread('./Figure4_Data.csv',',',1,0);
BasePairedResidues = M((M(:,3)== 1 | M(:,3)== 2 | M(:,3)== 3),1);
ImportantResidues = dlmread(strcat(Infolder,'ImportantResidues.csv'));
disp('Read RefSeq..')
Tmp = dlmread('./Data_CoVariation/RefSeq.txt');
RefSeq = Tmp(:,2);


NrPairs = dlmread(strcat(Infolder,'NumberPairs.csv'));
Pairs = zeros(NrPairs,2);

%Epistasis only for WC pairs
EpistasisWC = zeros(NrPairs,(NrPositionsTotal-1)*NrComparisons,3);
EpistasisPValueWC = ones(NrPairs,4);% WC - >WC : 3; WB -> WC: 2 
NrEstimatesWC = zeros(NrPairs,4);

%Epistasis only for WC or WB pairs
EpistasisWCWB = zeros(NrPairs,(NrPositionsTotal-1)*NrComparisons,5);
EpistasisPValueWCWB = ones(NrPairs,5);
NrEstimatesWCWB = zeros(NrPairs,5);

%Epistasis for all pairs
Epistasis = zeros(NrPairs,(NrPositionsTotal-1)*NrComparisons,9);
EpistasisPValue = ones(NrPairs,9);
NrEstimates = zeros(NrPairs,9);
AllCombinations = [1 1;1 2;1 3;1 4;2 1;2 2;2 3;2 4;3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];
AllIdx = (AllCombinations(:,1)-1).*4+AllCombinations(:,2);

RunningIndex = 0;
NativeBasePair = zeros(NrPositionsTotal,NrPositionsTotal);

for i = 1:length(ImportantResidues)
    Position1 = ImportantResidues(i);
    %Position1 = 245
    wtPos1 = RefSeq(Position1);
    filename2 = strcat(num2str(Position1),'_Partners.csv');
    try, %catch if Position1 has no partner
        Positions2 = dlmread(strcat(Infolder,filename2));
        %Positions2 = 275
        for j = 1:length(Positions2)
            RunningIndex = RunningIndex +1;  
            disp(strcat('Pair:#',num2str(RunningIndex),'/',num2str(NrPairs)));
            Position2 = Positions2(j);
            Pairs(RunningIndex,1) = Position1;
            Pairs(RunningIndex,2) = Position2;
            
            wtPos2 = RefSeq(Position2);
            filename3 = strcat(Infolder,num2str(Position1),'_',num2str(Position2),'.mat');
            load(filename3);

            %% post-processing based on quality criteria
            % a) Minimum Signal-2-Noise & b) Lower estimates correction
            %m_m_w
            %[n,m] = size(TotalRel_KD_m_m_w);
            LogicIdx = Signal2Noise_m_m_w_SuperN < MinimumSignal2NoiseStrength_m_m & Signal2Noise_m_m_w_Beads < MinimumSignal2NoiseStrength_m_m;
            TotalRel_KD_m_m_w(LogicIdx) = nan;

            LowerLimitsKD_m_m_w = zeros((NrPositionsTotal-1)*NrComparisons,16); 
            LogicIdx = Signal2Noise_m_m_w_SuperN >= MinimumSignal2NoiseStrength_m_m & Signal2Noise_m_m_w_Beads < MinimumSignal2NoiseStrength_m_m;
            LowerLimitsKD_m_m_w(LogicIdx) = 1;

            UpperLimitsKD_m_m_w = zeros((NrPositionsTotal-1)*NrComparisons,16); 
            LogicIdx = Signal2Noise_m_m_w_SuperN < MinimumSignal2NoiseStrength_m_m & Signal2Noise_m_m_w_Beads >= MinimumSignal2NoiseStrength_m_m;
            UpperLimitsKD_m_m_w(LogicIdx) = 1;
            
            %m_w_w
            LogicIdx = Signal2Noise_m_w_w_SuperN < MinimumSignal2NoiseStrength_m & Signal2Noise_m_w_w_Beads < MinimumSignal2NoiseStrength_m;
            TotalRel_KD_m_w_w(LogicIdx) = nan;

            LowerLimitsKD_m_w_w = zeros((NrPositionsTotal-1)*NrComparisons,4); 
            LogicIdx = Signal2Noise_m_w_w_SuperN >= MinimumSignal2NoiseStrength_m & Signal2Noise_m_w_w_Beads < MinimumSignal2NoiseStrength_m;
            LowerLimitsKD_m_w_w(LogicIdx) = 1;

            UpperLimitsKD_m_w_w = zeros((NrPositionsTotal-1)*NrComparisons,4); 
            LogicIdx = Signal2Noise_m_w_w_SuperN < MinimumSignal2NoiseStrength_m & Signal2Noise_m_w_w_Beads >= MinimumSignal2NoiseStrength_m;
            UpperLimitsKD_m_w_w(LogicIdx) = 1;
            
            %w_m_w
            LogicIdx = Signal2Noise_w_m_w_SuperN < MinimumSignal2NoiseStrength_m & Signal2Noise_w_m_w_Beads < MinimumSignal2NoiseStrength_m;
            TotalRel_KD_w_m_w(LogicIdx) = nan;

            LowerLimitsKD_w_m_w = zeros((NrPositionsTotal-1)*NrComparisons,4); 
            LogicIdx = Signal2Noise_w_m_w_SuperN >= MinimumSignal2NoiseStrength_m & Signal2Noise_w_m_w_Beads < MinimumSignal2NoiseStrength_m;
            LowerLimitsKD_w_m_w(LogicIdx) = 1;

            UpperLimitsKD_w_m_w = zeros((NrPositionsTotal-1)*NrComparisons,4); 
            LogicIdx = Signal2Noise_w_m_w_SuperN < MinimumSignal2NoiseStrength_m & Signal2Noise_w_m_w_Beads >= MinimumSignal2NoiseStrength_m;
            UpperLimitsKD_w_m_w(LogicIdx) = 1;

            % c) censoring by minimum coverage
            CensoredPositionsBeadsTmp = PositionWeightsBeads<WeightThreshold;
            CensoredPositionsSuperNTmp = PositionWeightsSuperN<WeightThreshold;
            CensoredPositionsTotalTmp = PositionWeightsTotal<WeightThreshold;
%             
%             CensoredPositionsBeadsTmp = PositionReadsBeads<ReadThreshold;
%             CensoredPositionsSuperNTmp = PositionReadsSuperN<ReadThreshold;
%             CensoredPositionsTotalTmp = PositionReadsTotal<ReadThreshold;
            
            CensoredPositionsBeads = false(1,(NrPositionsTotal-1)*NrComparisons);
            CensoredPositionsSuperN = false(1,(NrPositionsTotal-1)*NrComparisons);
            CensoredPositionsTotal = false(1,(NrPositionsTotal-1)*NrComparisons);


            for counter = 1:NrComparisons 
                idx = counter:NrComparisons:(NrPositionsTotal-1)*NrComparisons;
                CensoredPositionsBeads(idx) = CensoredPositionsBeadsTmp(:,counter)';
                CensoredPositionsSuperN(idx) = CensoredPositionsSuperNTmp(:,counter)';
                CensoredPositionsTotal(idx) = CensoredPositionsTotalTmp(:,counter)';
            end
            %
            CensoredPositionsFinal = max(max(CensoredPositionsBeads,CensoredPositionsSuperN),CensoredPositionsTotal);
            TotalRel_KD_m_m_w(CensoredPositionsFinal,:) = nan;
            LowerLimitsKD_m_m_w(CensoredPositionsFinal,:) = 0;
            UpperLimitsKD_m_m_w(CensoredPositionsFinal,:) = 0;
            
            TotalRel_KD_m_w_w(CensoredPositionsFinal,:) = nan;
            LowerLimitsKD_m_w_w(CensoredPositionsFinal,:) = 0;
            UpperLimitsKD_m_w_w(CensoredPositionsFinal,:) = 0;
            
            TotalRel_KD_w_m_w(CensoredPositionsFinal,:) = nan;
            LowerLimitsKD_w_m_w(CensoredPositionsFinal,:) = 0;
            UpperLimitsKD_w_m_w(CensoredPositionsFinal,:) = 0;
            

            % d) replace arbitrary set lower bounds for KD with realistic values.

            for k = 1:16
                Positions_m_m_w = TotalRel_KD_m_m_w(:,k) == -1 | TotalRel_KD_m_m_w(:,k) == -10;%lower and upper estimates

                minmax_m_m_w = prctile(TotalRel_KD_m_m_w(~Positions_m_m_w,k),[5 50 95]);
                Positions_m_m_w_Low = TotalRel_KD_m_m_w(:,k) == -1;
                Positions_m_m_w_Up = TotalRel_KD_m_m_w(:,k) == -10;
                TotalRel_KD_m_m_w(Positions_m_m_w_Low,k) = max(minmax_m_m_w(3),ArbitraryValueLowerEstimate); %replace with median
                TotalRel_KD_m_m_w(Positions_m_m_w_Up,k) = max(minmax_m_m_w(1),ArbitraryValueUpperEstimate); %replace with median
            end
            
           for k = 1:4
                Positions_m_w_w = TotalRel_KD_m_w_w(:,k) == -1 | TotalRel_KD_m_w_w(:,k) == -10;%lower and upper estimates
                minmax_m_w_w = prctile(TotalRel_KD_m_w_w(~Positions_m_w_w,k),[5 50 95]);
                Positions_m_w_w_Low = TotalRel_KD_m_w_w(:,k) == -1;
                Positions_m_w_w_Up = TotalRel_KD_m_w_w(:,k) == -10;
                TotalRel_KD_m_w_w(Positions_m_w_w_Low,k) = max(minmax_m_w_w(3),ArbitraryValueLowerEstimate); %replace with median
                TotalRel_KD_m_w_w(Positions_m_w_w_Up,k) = max(minmax_m_w_w(1),ArbitraryValueUpperEstimate); %replace with median
                
                Positions_w_m_w = TotalRel_KD_w_m_w(:,k) == -1 | TotalRel_KD_w_m_w(:,k) == -10;%lower and upper estimates
                minmax_w_m_w = prctile(TotalRel_KD_w_m_w(~Positions_w_m_w,k),[5 50 95]);
                Positions_w_m_w_Low = TotalRel_KD_w_m_w(:,k) == -1;
                Positions_w_m_w_Up = TotalRel_KD_w_m_w(:,k) == -10;
                TotalRel_KD_w_m_w(Positions_w_m_w_Low,k) = max(minmax_w_m_w(3),ArbitraryValueLowerEstimate); %replace with median
                TotalRel_KD_w_m_w(Positions_w_m_w_Up,k) = max(minmax_w_m_w(1),ArbitraryValueUpperEstimate); %replace with median 
            end
            %% Epistasis calculation

            wtPos1Pos2 = 4*(wtPos1-1)+wtPos2;
            NativeBasePair(Position1,Position2) = wtPos1Pos2;
            
            % 1) Only for watson-crick pairs
            IDX = find(WCBasePairs~=wtPos1Pos2);
            %WCDoubleMutPairs = WCBasePairs(WCBasePairs~=wtPos1Pos2);
            EpiIndex = 0;
            for k = IDX;%1:length(WCDoubleMutPairs)
                EpiIndex =  EpiIndex+1;
                m_m_index = WCBasePairs(k);
                m_w_index = WCSingleBases(k,1);
                w_m_index = WCSingleBases(k,2); 
                EpistasisWC(RunningIndex,:,EpiIndex) = log(1./TotalRel_KD_m_m_w(:,m_m_index)) - log(1./TotalRel_KD_m_w_w(:,m_w_index)) - log(1./TotalRel_KD_w_m_w(:,w_m_index));
                IgnoreIDX = isnan(EpistasisWC(RunningIndex,:,EpiIndex)) | isinf(EpistasisWC(RunningIndex,:,EpiIndex));
                EpistasisWC(RunningIndex,IgnoreIDX,EpiIndex) = nan;
                NrEstimatesWC(RunningIndex,EpiIndex) = sum(~isnan(EpistasisWC(RunningIndex,:,EpiIndex)));
                if NrEstimatesWC(RunningIndex,EpiIndex) >= MinimumNumberEstimatableKDs %&& isreal(Epi)
                    EpistasisPValueWC1 = (sum(EpistasisWC(RunningIndex,:,EpiIndex)<=0)./NrEstimatesWC(RunningIndex,EpiIndex));
                    EpistasisPValueWC2 = (sum(EpistasisWC(RunningIndex,:,EpiIndex)>=0)./NrEstimatesWC(RunningIndex,EpiIndex));
                    EpistasisPValueWC(RunningIndex,EpiIndex) = min(EpistasisPValueWC1,EpistasisPValueWC2);
                end
            end
            
            % 1) Only for watson-crick or wobble pairs
            IDX = find(WCWBBasePairs~=wtPos1Pos2);
            %WCDoubleMutPairs = WCBasePairs(WCBasePairs~=wtPos1Pos2);
            EpiIndex = 0;
            for k = IDX;%1:length(WCDoubleMutPairs)
                EpiIndex =  EpiIndex+1;
                m_m_index = WCWBBasePairs(k);
                m_w_index = WCWBSingleBases(k,1);
                w_m_index = WCWBSingleBases(k,2); 
                EpistasisWCWB(RunningIndex,:,EpiIndex) = log(1./TotalRel_KD_m_m_w(:,m_m_index)) - log(1./TotalRel_KD_m_w_w(:,m_w_index)) - log(1./TotalRel_KD_w_m_w(:,w_m_index));
                IgnoreIDX = isnan(EpistasisWCWB(RunningIndex,:,EpiIndex)) | isinf(EpistasisWCWB(RunningIndex,:,EpiIndex));
                EpistasisWCWB(RunningIndex,IgnoreIDX,EpiIndex) = nan;
                NrEstimatesWCWB(RunningIndex,EpiIndex) = sum(~isnan(EpistasisWCWB(RunningIndex,:,EpiIndex)));
                if NrEstimatesWCWB(RunningIndex,EpiIndex) >= MinimumNumberEstimatableKDs %&& isreal(Epi)
                    EpistasisPValueWCWB1 = (sum(EpistasisWCWB(RunningIndex,:,EpiIndex)<=0)./NrEstimatesWCWB(RunningIndex,EpiIndex));
                    EpistasisPValueWCWB2 = (sum(EpistasisWCWB(RunningIndex,:,EpiIndex)>=0)./NrEstimatesWCWB(RunningIndex,EpiIndex));
                    EpistasisPValueWCWB(RunningIndex,EpiIndex) = min(EpistasisPValueWCWB1,EpistasisPValueWCWB2);
                end
            end
            
            % 2) For all double mutants
            IDX = find(all([AllCombinations(:,1)~=wtPos1,AllCombinations(:,2)~=wtPos2]'));
            EpiIndex = 0;
            for k = IDX
                EpiIndex =  EpiIndex+1;
                m_m_index = AllIdx(k);
                m_w_index = AllCombinations(k,1);
                w_m_index = AllCombinations(k,2); 
                Epistasis(RunningIndex,:,EpiIndex) = log(1./TotalRel_KD_m_m_w(:,m_m_index)) - log(1./TotalRel_KD_m_w_w(:,m_w_index)) - log(1./TotalRel_KD_w_m_w(:,w_m_index));
                IgnoreIDX = isnan(Epistasis(RunningIndex,:,EpiIndex)) | isinf(Epistasis(RunningIndex,:,EpiIndex));
                Epistasis(RunningIndex,IgnoreIDX,EpiIndex) = nan;
                NrEstimates(RunningIndex,k) = sum(~isnan(Epistasis(RunningIndex,:,EpiIndex)));
                if NrEstimates(RunningIndex,EpiIndex) >= MinimumNumberEstimatableKDs %&& isreal(Epi)
                    %test for sign. increase
                    EpistasisPValue1 = (sum(Epistasis(RunningIndex,:,EpiIndex)<=0)./NrEstimates(RunningIndex,EpiIndex));
                    %test for sign. decrease
                    EpistasisPValue2 = (sum(Epistasis(RunningIndex,:,EpiIndex)>=0)./NrEstimates(RunningIndex,EpiIndex));
                    EpistasisPValue(RunningIndex,EpiIndex)= min(EpistasisPValue1,EpistasisPValue2);
                end
            end
            
            %return
            clear TotalRel_KD_m_w_w TotalRel_KD_w_m_w TotalRel_KD_m_m_w Signal2Noise_m_w_w_SuperN Signal2Noise_m_w_w_Beads
            clear Signal2Noise_w_m_w_SuperN Signal2Noise_w_m_w_Beads Signal2Noise_m_m_w_SuperN Signal2Noise_m_m_w_Beads PositionWeightsBeads
            clear PositionWeightsSuperN PositionWeightsTotal PositionReadsBeads PositionReadsSuperN PositionReadsTotal
        end
    catch
    end
end
%%



SignificantInteractionsWC = any(EpistasisPValueWC'<alpha);


%1) Watson-Crick Pairs
PotentialPairsWC = Pairs(SignificantInteractionsWC,:);
NrPotentialPairsWC = length(PotentialPairsWC);
MedianEpistasisWC = zeros(NrPotentialPairsWC,3);
PrcTileEpistasisWC = zeros(NrPotentialPairsWC,2,3);

EpistasisPValueOut = zeros(NrPotentialPairsWC,3);
NrEstimatesOut = zeros(NrPotentialPairsWC,3);
for k = 1:3
    EpistasisPValueOut(:,k) = EpistasisPValueWC(SignificantInteractionsWC,k);
    NrEstimatesOut(:,k) = NrEstimatesWC(SignificantInteractionsWC,k);
end


for k = 1:3
    MedianEpistasisWC(:,k) = nanmedian(EpistasisWC(SignificantInteractionsWC,:,k),2);
    PrcTileEpistasisWC(:,1:2,k) = (prctile(EpistasisWC(SignificantInteractionsWC,:,k)',[5 95]))';
end

MaxMedianEpistasisWC = max(MedianEpistasisWC,[],2);
%[PotentialPairsWC,MedianEpistasisWC];

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%2) Watson-Crick and Wobble Pairs
SignificantInteractionsWCWB = any(EpistasisPValueWCWB'<alpha);
PotentialPairsWCWB = Pairs(SignificantInteractionsWCWB,:);
NrPotentialPairsWCWB = length(PotentialPairsWCWB);
MedianEpistasisWCWB = zeros(NrPotentialPairsWCWB,3);
PrcTileEpistasisWCWB = zeros(NrPotentialPairsWCWB,2,3);

for k = 1:5
    MedianEpistasisWCWB(:,k) = nanmedian(EpistasisWCWB(SignificantInteractionsWCWB,:,k),2);
    PrcTileEpistasisWCWB(:,1:2,k) = (prctile(EpistasisWCWB(SignificantInteractionsWCWB,:,k)',[5 95]))';
end

MaxMedianEpistasisWCWB = max(MedianEpistasisWCWB,[],2);
%[PotentialPairsWCWB,MedianEpistasisWCWB];
%[Pairs,EpistasisPValueWCWB]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%3) All Double Mutation Pairs
SignificantInteractions = any(EpistasisPValue'<alpha);
%SignificantInteractionsWC = true(1,length(EpistasisPValueWC));
%SignificantInteractionsWC = all(EpistasisPValueWC'<alpha);
PotentialPairs = Pairs(SignificantInteractions,:);
NrPotentialPairs = length(PotentialPairs);
MedianEpistasis = zeros(NrPotentialPairs,3);
PrcTileEpistasis = zeros(NrPotentialPairs,2,3);

for k = 1:9
    MedianEpistasis(:,k) = nanmedian(Epistasis(SignificantInteractions,:,k),2);
    PrcTileEpistasis(:,1:2,k) = (prctile(Epistasis(SignificantInteractions,:,k)',[5 95]))';
end
MaxMedianEpistasis = max(MedianEpistasis,[],2);
%[PotentialPairs,MaxMedianEpistasis];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Write Output (WC Pairs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Output = zeros(NrPotentialPairsWC,6*3+3);
if WriteOutput
    Output(:,1:2) = PotentialPairsWC;
    for i = 1:NrPotentialPairsWC
        idx = find(WCSingleBases(:,1)' == RefSeq(PotentialPairsWC(i,1)));
        WCPair = WCBasePairs(idx);
        Output(i,3) = WCBasePairs(idx);%wt pair
        WCDoubleMutPairs = WCBasePairs(WCBasePairs~=WCPair);
        for k = 1:length(WCDoubleMutPairs)
            m_m_index = WCDoubleMutPairs(k);
            Output(i,(k-1)*6+3+1) = m_m_index;
            Output(i,(k-1)*6+3+2) = MedianEpistasisWC(i,k);
            Output(i,(k-1)*6+3+3) = EpistasisPValueOut(i,k);
            Output(i,(k-1)*6+3+4) = NrEstimatesOut(i,k);
            Output(i,(k-1)*6+3+[5:6]) = PrcTileEpistasisWC(i,:,k);
        end
    end
    filename = strcat('./EpistasisEstimates.csv');
    try,
        delete(filename);
    catch
    end
    %generate Header
    Header = {'First Pos.';'Second Pos.';'wt Pair';...
        'WC mut. Pair #1'; 'median'; 'p-value';'# estimates';'5th prctile';'95th prctile';...
        'WC mut. Pair #2'; 'median'; 'p-value';'# estimates';'5th prctile';'95th prctile';...
        'WC mut. Pair #3'; 'median'; 'p-value';'# estimates';'5th prctile';'95th prctile'};
    
    fid = fopen(filename,'w');
    for i = 1:length(Header)-1
        fprintf(fid,'%s,',char(Header(i)));
    end
    fprintf(fid,'%s\n',char(Header(end)));
    fclose(fid);
    %write Output values
    dlmwrite(filename,Output,'-append')
end

%% Grafics --------------------------------------------------------------------------
% WC Pairs only
ColorMapEpistasis = zeros(NrPositionsTotal+1,NrPositionsTotal+1);
for i = 1:NrPotentialPairsWC
    ColorMapEpistasis(PotentialPairsWC(i,1)-FirstPos+1,PotentialPairsWC(i,2)-FirstPos+1) = MaxMedianEpistasisWC(i);
end
%----------------

IDX = sum((PotentialPairsWC <= DisRegion(end) & PotentialPairsWC >= DisRegion(1)),2) >= 2;
PotentialPairsDis = PotentialPairsWC(IDX,:);
MedianEpistasisDis = MedianEpistasisWC(IDX,:);

ColorMapEpistasisDis = zeros(length(DisRegion)*3,length(DisRegion)*3);
for i = 1:length(PotentialPairsDis)
    Pos1 = PotentialPairsDis(i,1);
    Pos2 = PotentialPairsDis(i,2);
    MutPos1 = Bases(Bases~=RefSeq(Pos1));
    MutPos2 = Bases(Bases~=RefSeq(Pos2));
    counter = 1;   
    for k1 = 1:3
        for k2 = 1:3
            idx_i = (PotentialPairsDis(i,1)-DisRegion(1))*3+k1;
            idx_j = (PotentialPairsDis(i,2)-DisRegion(1))*3+k2;
            if any((MutPos1(k1)-1)*4+MutPos2(k2) == WCBasePairs)
                ColorMapEpistasisDis(idx_i,idx_j) = MedianEpistasisDis(i,counter);
                counter = counter +1;
            end
        end
    end
end
%-------------------G:U Pair--------------------
% if Fix253_265
%     ColorMapEpistasis(253-FirstPos+1,265-FirstPos+1) = 10;
%     idx_i = (253-DisRegion(1))*3+[1 3 3];
%     idx_j = (265-DisRegion(1))*3+[3 2 3];
%     for i = 1:length(idx_i)
%         ColorMapEpistasisDis(idx_i(i),idx_j(i)) = 10;
%     end
% end
%------------------------------------------------

figure(7)
image([FirstPos LastPos],[FirstPos LastPos],ColorMapEpistasis.*30)
colormap(gray)
hold on
plot([FirstPos LastPos],[FirstPos LastPos],'-w','LineWidth',3)
xlabel('first residue')
ylabel('second residue')
title('Epistasis significantly > 0 (only wc)')
%colorbar
fett(7)

ColorMapEpistasisDisSimple =  ColorMapEpistasis([DisStart:DisEnd]-FirstPos+1,[DisStart:DisEnd]-FirstPos+1);
ColorMapEpistasisDisSimple = ColorMapEpistasisDisSimple'+ColorMapEpistasisDisSimple;
figure(5)
image(DisRegion,fliplr(DisRegion),ColorMapEpistasisDisSimple.*5)
%colormap gray
hold on
plot(DisRegion,fliplr(DisRegion),'-w','LineWidth',3)
xlabel('first residue')
ylabel('second residue')
title('WC Epistasis significantly > 0 (WC pairs); Zoom on Dis-region')
colorbar
set(gca,'XTick',(DisStart:3:DisEnd)-0.5,'XTickLabel',num2str((DisStart:3:DisEnd)'))
set(gca,'YTick',(DisStart:3:DisEnd)-0.5,'YTickLabel',num2str((DisEnd:-3:DisStart)'))
fett(5)

ColorMapEpistasisDis = ColorMapEpistasisDis'+ColorMapEpistasisDis;

figure(6)
first = DisStart;
last = DisEnd+((DisEnd-DisStart + 1)*2);
axisvector = [first:last];
image(axisvector,fliplr(axisvector),ColorMapEpistasisDis.*20)
%colormap gray
hold on
plot(axisvector,fliplr(axisvector),'-w','LineWidth',3)
BasesChar = ['A';'C';'G';'U'];
bases = 1:4;
counter1 = 1;
XLabel = repmat(' ',(DisEnd-DisStart)*3,1); 
for i = DisStart:DisEnd
    wtbase = RefSeq(i);
    XLabel(counter1:counter1+2) = BasesChar(bases~=wtbase);
    counter1 = counter1 +3;
end
XLabelRev = flipud(XLabel);
set(gca,'XTick',(first:1:last),'XTickLabel',XLabel,'TickDir','in','FontSize',6,'FontWeight','bold')
set(gca,'YTick',(first:1:last),'YTickLabel',XLabelRev,'FontSize',6,'FontWeight','bold')
%add grid
for i = 1:3:length(axisvector)
    line([axisvector(i) axisvector(i)]-0.5,[first last],'LineStyle',':','Color','w')
    line([first last],[axisvector(i) axisvector(i)]-0.5,'LineStyle',':','Color','w')
    text(axisvector(i),last+4,num2str(DisRegion(ceil(i/3))),'FontWeight','bold')
    text(first-3.5,axisvector(i)+1,num2str(DisRegion(DisEnd-DisStart-ceil(i/3)+2)),'FontWeight','bold')
end
title('WC Epistasis significantly > 0 (WC pairs); Zoom on Dis-region', 'FontSize',14,'FontWeight','bold')
%% Grafics --------------------------------------------------------------------------
% WC & WB Pairs only
NrPositionsTotal = LastPos-FirstPos+1;

ColorMapEpistasis = zeros(NrPositionsTotal,NrPositionsTotal);
for i = 1:NrPotentialPairsWCWB
    ColorMapEpistasis(PotentialPairsWCWB(i,1)-FirstPos+1,PotentialPairsWCWB(i,2)-FirstPos+1) = MaxMedianEpistasisWCWB(i);
end
%----------------

IDX = sum((PotentialPairsWCWB <= DisRegion(end) & PotentialPairsWCWB >= DisRegion(1)),2) >= 2;
PotentialPairsDis = PotentialPairsWCWB(IDX,:);
MedianEpistasisDis = MedianEpistasisWCWB(IDX,:);

ColorMapEpistasisDis = zeros(length(DisRegion)*3,length(DisRegion)*3);
for i = 1:length(PotentialPairsDis)
    Pos1 = PotentialPairsDis(i,1);
    Pos2 = PotentialPairsDis(i,2);
    MutPos1 = Bases(Bases~=RefSeq(Pos1));
    MutPos2 = Bases(Bases~=RefSeq(Pos2));
    counter = 1;   
    for k1 = 1:3
        for k2 = 1:3
            idx_i = (PotentialPairsDis(i,1)-DisRegion(1))*3+k1;
            idx_j = (PotentialPairsDis(i,2)-DisRegion(1))*3+k2;
            if any((MutPos1(k1)-1)*4+MutPos2(k2) == WCWBBasePairs)
                ColorMapEpistasisDis(idx_i,idx_j) = MedianEpistasisDis(i,counter);
                counter = counter +1;
            end
        end
    end
end
%-------------------G:U Pair--------------------
% if Fix253_265
%     ColorMapEpistasis(253-FirstPos+1,265-FirstPos+1) = 10;
%     idx_i = (253-DisRegion(1))*3+[1 3 3];
%     idx_j = (265-DisRegion(1))*3+[3 2 3];
%     for i = 1:length(idx_i)
%         ColorMapEpistasisDis(idx_i(i),idx_j(i)) = 10;
%     end
% end
%------------------------------------------------

ColorMapEpistasisPlot = ColorMapEpistasis'+ColorMapEpistasis;
dlmwrite('Figure4bData.txt',[[nan,(FirstPos:LastPos)]; [[(FirstPos:LastPos)]',ColorMapEpistasisPlot]]);
%
xvalues = [FirstPos LastPos];
figure(2)
image(xvalues,fliplr(xvalues),ColorMapEpistasisPlot.*30)
colormap(gray)
hold on
plot(xvalues,fliplr(xvalues),'-w','LineWidth',2)
% for i = 1:length(BasePairedResidues)
%     line([BasePairedResidues(i) BasePairedResidues(i)]-0.5,xvalues,'Color','r','LineStyle',':')
%     line(xvalues,[BasePairedResidues(i) BasePairedResidues(i)]+0.5,'Color','r','LineStyle',':')
% end
title('Epistasis significantly > 0 (only wc & wb)','FontSize',24,'FontWeight','bold')
xlabel('first residue','FontSize',24,'FontWeight','bold')
ylabel('second residue','FontSize',24,'FontWeight','bold')
set(gca,'XTick',(FirstPos:20:LastPos)+0.5,'XTickLabel',num2str((FirstPos:20:LastPos)'),'FontSize',18,'FontWeight','bold')
set(gca,'YTick',(FirstPos:20:LastPos)-0.5,'YTickLabel',num2str((LastPos:-20:FirstPos)'),'FontSize',18,'FontWeight','bold')


h = colorbar;
Ticks = get(h,'YTick');
Ticks = [1 Ticks(end)];
%maxim = max(max(ColorMapEpistasisPlot));
%vec = [1 maxim/2 maxim];
colorbar('YTick',Ticks,'YTickLabel',num2str([0; round((Ticks(end)./30.*10)./10)]),'FontSize',24,'FontWeight','bold');
print(2,'-depsc2','EpistasisOverview.eps')

%set(h,'YTick',vec,'YTickLabel',num2str(round(vec'./30).*10))
%fett(2)

ColorMapEpistasisDisSimple =  ColorMapEpistasis([DisStart:DisEnd]-FirstPos+1,[DisStart:DisEnd]-FirstPos+1);
ColorMapEpistasisDisSimple = ColorMapEpistasisDisSimple'+ColorMapEpistasisDisSimple;
figure(3)
image(DisRegion,fliplr(DisRegion),ColorMapEpistasisDisSimple.*5)
%colormap gray
hold on
plot(DisRegion,fliplr(DisRegion),'-w','LineWidth',3)
xlabel('first residue')
ylabel('second residue')
title('WC Epistasis significantly > 0 (WC & WB pairs); Zoom on Dis-region')
colorbar
set(gca,'XTick',(DisStart:3:DisEnd)-0.5,'XTickLabel',num2str((DisStart:3:DisEnd)'))
set(gca,'YTick',(DisStart:3:DisEnd)-0.5,'YTickLabel',num2str((DisEnd:-3:DisStart)'))
fett(3)

ColorMapEpistasisDis = ColorMapEpistasisDis'+ColorMapEpistasisDis;

figure(4)
first = DisStart;
last = DisEnd+((DisEnd-DisStart + 1)*2);
axisvector = [first:last];
image(axisvector,fliplr(axisvector),ColorMapEpistasisDis.*20)
%colormap gray
hold on
plot(axisvector,fliplr(axisvector),'-w','LineWidth',3)
BasesChar = ['A';'C';'G';'U'];
bases = 1:4;
counter1 = 1;
XLabel = repmat(' ',(DisEnd-DisStart)*3,1); 
for i = DisStart:DisEnd
    wtbase = RefSeq(i);
    XLabel(counter1:counter1+2) = BasesChar(bases~=wtbase);
    counter1 = counter1 +3;
end
XLabelRev = flipud(XLabel);
set(gca,'XTick',(first:1:last),'XTickLabel',XLabel,'TickDir','in','FontSize',6,'FontWeight','bold')
set(gca,'YTick',(first:1:last),'YTickLabel',XLabelRev,'FontSize',6,'FontWeight','bold')
%add grid
for i = 1:3:length(axisvector)
    line([axisvector(i) axisvector(i)]-0.5,[first last],'LineStyle',':','Color','w')
    line([first last],[axisvector(i) axisvector(i)]-0.5,'LineStyle',':','Color','w')
    text(axisvector(i),last+4,num2str(DisRegion(ceil(i/3))),'FontWeight','bold')
    text(first-3.5,axisvector(i)+1,num2str(DisRegion(DisEnd-DisStart-ceil(i/3)+2)),'FontWeight','bold')
end
title('WC Epistasis significantly > 0 (WC & WB pairs); Zoom on Dis-region', 'FontSize',14,'FontWeight','bold')

%% All Mutations
IDX = sum((PotentialPairs <= DisRegion(end) & PotentialPairs >= DisRegion(1)),2) >= 2;
PotentialPairsDis = PotentialPairs(IDX,:);
MedianEpistasisDis = MedianEpistasis(IDX,:);

ColorMapEpistasisDis = zeros(length(DisRegion)*3,length(DisRegion)*3);
for i = 1:length(PotentialPairsDis)
    for k1 = 1:3
        for k2 = 1:3
            idx_i = (PotentialPairsDis(i,1)-DisRegion(1))*3+k1;
            idx_j = (PotentialPairsDis(i,2)-DisRegion(1))*3+k2;
            ColorMapEpistasisDis(idx_i,idx_j) = MedianEpistasisDis(i,(k1-1)*3+k2);
        end
    end
end
% %-------------------G:U Pair--------------------
% if Fix253_265
%     ColorMapEpistasis(253-FirstPos+1,265-FirstPos+1) = 10;
%     idx_i = (253-DisRegion(1))*3+[1 3 3];
%     idx_j = (265-DisRegion(1))*3+[3 2 3];
%      for i = 1:length(idx_i)
%         ColorMapEpistasisDis(idx_i(i),idx_j(i)) = 10;
%     end
% end
%------------------------------------------------

ColorMapEpistasisDis = ColorMapEpistasisDis'+ColorMapEpistasisDis;

figure(1)
first = 243;
last = 277+((277-243 + 1)*2);
axisvector = [first:last];
image(axisvector,fliplr(axisvector),ColorMapEpistasisDis.*20)
%colormap gray
hold on
plot(axisvector,fliplr(axisvector),'-w','LineWidth',3)
set(gca,'XTick',(first:1:last),'XTickLabel',XLabel,'TickDir','in','FontSize',6,'FontWeight','bold')
set(gca,'YTick',(first:1:last),'YTickLabel',XLabelRev,'FontSize',6,'FontWeight','bold')
for i = 1:3:length(axisvector)
    line([axisvector(i) axisvector(i)]-0.5,[first last],'LineStyle',':','Color','w')
    line([first last],[axisvector(i) axisvector(i)]-0.5,'LineStyle',':','Color','w')
    text(axisvector(i),last+4,num2str(DisRegion(ceil(i/3))),'FontWeight','bold')
    text(first-4,axisvector(i)+1,num2str(DisRegion(DisEnd-DisStart-ceil(i/3)+2)),'FontWeight','bold')
end
title('Epistasis significantly > 0 (all double mutants); Zoom on Dis-region','FontSize',14,'FontWeight','bold')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Stems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%filter for maximum likelihood interactions
LastIdx = find([PotentialPairsWCWB(2:end,1)-PotentialPairsWCWB(1:end-1,1);NrPotentialPairsWCWB]);
FirstIdx = [1;LastIdx(1:end-1)+1];


%Do Regularization (a penalty for each added position  to the stem, where
%there is actually no epistasis value)
K = ones(NrPositionsTotal);
%Option I. Leave out if all extentions should be penalized
ColorMapEpistasis(isnan(ColorMapEpistasis)) = 0;
K(logical(ColorMapEpistasis)) = 0;
EpiNew2 = ColorMapEpistasis + K*penalty;
%Smoothing ------------------------------------------------------------
if smooth == 1
    disp('smoothing...')
    EpiNew2 = smooth2d(EpiNew2,1);
end
%Smoothing ------------------------------------------------------------
EpiNew2 = fliplr(EpiNew2);

AllStarts_i = zeros(1,NrPositionsTotal);
AllStarts_j = zeros(1,NrPositionsTotal);
AllEnds_i = zeros(1,NrPositionsTotal);
AllEnds_j = zeros(1,NrPositionsTotal);
AllScores = zeros(1,NrPositionsTotal);

i = NrPositionsTotal; %index in first dimension (row) 
j = 1; %index in 2nd dimension (collumn)
counter = 1;

%go through all diagonals (possible stems) and assign start & end positions
for i = -NrPositionsTotal+2:1:NrPositionsTotal-2
    %all diagonals 
    DiagonalElements = diag(EpiNew2,i); %corresponding elements
    NrDiagonalElements = length(DiagonalElements);
   
    if any(DiagonalElements ~= penalty) % so that one doesn't need to scan through all parameters
        % For each diagonal = possible stem, go through all subsets of the diagonal K:L
        L = (DiagonalElements>penalty)';%logical indices of entries
        if sum(L) > 0 % there are entries
            %start and end positions of consecutive entries (poss. stems)
            K = [0 L(1:end)]-[L(1:end) 0];
            Starts = find(K == -1);
            Ends = find(K == 1)-1;
            for k = 1:length(Starts)
                if Ends(k)-Starts(k) > MinimumStemSize-1 %if minimum length 2; store start and end  
                    if i < 0
                        idx_i = -i+1:NrPositionsTotal;
                        idx_j = 1:length(idx_i);
                    elseif i == 0
                        idx_i = 1:NrPositionsTotal;
                        idx_j = 1:NrPositionsTotal;
                    elseif i > 0
                        idx_j = i+1:NrPositionsTotal;
                        idx_i = 1:length(idx_j);
                    end
                    AllStarts_i(counter) = idx_i(Starts(k));
                    AllEnds_i(counter) = idx_i(Ends(k));
                    AllStarts_j(counter) = idx_j(Starts(k));
                    AllEnds_j(counter) = idx_j(Ends(k));
                    AllScores(counter) = sum(DiagonalElements(Starts(k):Ends(k)));
                    counter = counter+1;
                end
            end
        end
    end
end
%trim Outputs
AllScores = AllScores(logical(AllScores));
AllStarts_j = AllStarts_j(logical(AllStarts_j));
AllStarts_i = AllStarts_i(logical(AllStarts_i));
AllEnds_j = AllEnds_j(logical(AllEnds_j));
AllEnds_i = AllEnds_i(logical(AllEnds_i));

IDX = AllScores>0;
AllScores = AllScores(IDX);
AllStarts_j = AllStarts_j(IDX);
AllStarts_i = AllStarts_i(IDX);
AllEnds_j = AllEnds_j(IDX);
AllEnds_i = AllEnds_i(IDX);
%the left-right flipping has to be reversed
AllStarts_j = NrPositionsTotal -AllStarts_j + 1;
AllEnds_j = NrPositionsTotal -AllEnds_j + 1;


% % Post-Check if all identified pairs are Watson-Crick
% FinalColorMap = zeros(NrPositionsTotal); % for illustration
% NrIdentified = length(AllStarts_i);
% IsWatsonCrick = ones(1,NrIdentified);
% for i = 1: NrIdentified
%     LengthOfStem = AllEnds_i(i) - AllStarts_i(i) + 1;
%     idx_i = AllStarts_i(i):1:AllEnds_i(i);
%     idx_j = AllStarts_j(i):-1:AllEnds_j(i);
%     for j = 1:LengthOfStem 
%         if ~(any(WCBasePairs == NativeBasePair(idx_i(j),idx_j(j))))%not a WC base Pair
%             IsWatsonCrick(i) = 0;
%         end
%         FinalColorMap(idx_i(j),idx_j(j)) = AllScores(i);
%     end
%     if ~IsWatsonCrick(i) % those stems where non-WC base pairings have been introduced are removed
%         for j = 1:LengthOfStem
%             FinalColorMap(idx_i(j),idx_j(j)) = 0;
%         end
%     end
% end
% if ~all(IsWatsonCrick)
%     disp('Increase penalty value and ReRun') % warning
% end



%-----------------All stems with length > theta have been selected--------

%
FinalColorMap2 = zeros(NrPositionsTotal); % for illustration
NrStems = length(AllScores);
% % 6) Check for Conflics
ResiduesI = {}; % for resolving Conflicts
ResiduesJ = {}; % for resolving Conflicts
for i = 1: NrStems
    LengthOfStem = AllEnds_i(i) - AllStarts_i(i) + 1;
    idx_i = AllStarts_i(i):1:AllEnds_i(i);
    idx_j = AllStarts_j(i):-1:AllEnds_j(i);
    ResiduesI.(strcat('S',num2str(i))) = idx_i;
    ResiduesJ.(strcat('S',num2str(i))) = idx_j;
    %if IsWatsonCrick(i)
        for j = 1:LengthOfStem 
            FinalColorMap2(idx_i(j),idx_j(j)) = AllScores(i);
        end
    %end
end
FinalColorMap2Plot = FinalColorMap2' + FinalColorMap2;

figure(8)
image(xvalues,fliplr(xvalues),FinalColorMap2Plot.*4)
colormap gray
hold on
plot([FirstPos LastPos],fliplr([FirstPos LastPos]),'-w','LineWidth',2)
for i = 1:length(BasePairedResidues)
    line([BasePairedResidues(i) BasePairedResidues(i)]+0.4,[xvalues(1) xvalues(1)+4],'Color','r','LineStyle','-','LineWidth',4)
    line([BasePairedResidues(i) BasePairedResidues(i)]+0.4,[xvalues(2) xvalues(2)-4],'Color','r','LineStyle','-','LineWidth',4)
    line([xvalues(1) xvalues(1)+4],LastPos-[BasePairedResidues(i) BasePairedResidues(i)]+FirstPos,'Color','w','LineStyle','-','LineWidth',4)
    line([xvalues(2) xvalues(2)-4],LastPos-[BasePairedResidues(i) BasePairedResidues(i)]+FirstPos,'Color','w','LineStyle','-','LineWidth',4)
end

for i = 1: NrStems
    idx_j =ResiduesI.(strcat('S',num2str(i)))+FirstPos-1; 
    idx_i = ResiduesJ.(strcat('S',num2str(i)))+FirstPos-1;
    for j = 1:length(idx_j)
        if any(BasePairedResidues == idx_j(j)) || any(BasePairedResidues == idx_i(j))
            plot(idx_j(j)+0.3,LastPos-idx_i(j)+FirstPos-0.5,'wo','MarkerSize',10);%,'MarkerFaceColor','w')
            plot(idx_i(j)+0.3,LastPos-idx_j(j)+FirstPos-0.5,'wo','MarkerSize',10);%,,'MarkerFaceColor','w')
%         elseif any(BasePairedResidues == idx_i(j))
%             plot(idx_j(j),LastPos-idx_i(j)+FirstPos,'dr','MarkerFaceColor
%             ','r')
        end
    end
end

title('Stem score','FontSize',24,'FontWeight','bold')
xlabel('first residue','FontSize',24,'FontWeight','bold')
ylabel('second residue','FontSize',24,'FontWeight','bold')
% set(gca,'XTick',(FirstPos:2:LastPos)+0.5,'XTickLabel',num2str((FirstPos:2:LastPos)'),'FontSize',18,'FontWeight','bold')
% set(gca,'YTick',(FirstPos:2:LastPos)-0.5,'YTickLabel',num2str((LastPos:-2:FirstPos)'),'FontSize',18,'FontWeight','bold')
set(gca,'XTick',(FirstPos:20:LastPos)+0.5,'XTickLabel',num2str((FirstPos:20:LastPos)'),'FontSize',18,'FontWeight','bold')
set(gca,'YTick',(FirstPos:20:LastPos)-0.5,'YTickLabel',num2str((LastPos:-20:FirstPos)'),'FontSize',18,'FontWeight','bold')
h = colorbar;
Ticks = get(h,'YTick');
Ticks = [1 Ticks(end)];
%maxim = max(max(FinalColorMap2Plot));
%vec = [1 maxim/2 maxim];
colorbar('YTick',Ticks,'YTickLabel',num2str([0; round((Ticks(end)./4.*10)./10)]),'FontSize',24,'FontWeight','bold');
print(8,'-depsc2','Stems1.eps')


borderLOW = 247;
borderUP = 256;
ColorR = BasePairedResidues(BasePairedResidues <= borderUP);

figure(10)
image(xvalues,fliplr(xvalues),FinalColorMap2Plot.*4)
colormap gray
hold on
plot([FirstPos LastPos],fliplr([FirstPos LastPos]),'-w','LineWidth',2)
for i = 1:length(BasePairedResidues)
    line([BasePairedResidues(i) BasePairedResidues(i)]+0.4,[-254 -251]+LastPos+FirstPos+1,'Color','r','LineStyle','-','LineWidth',4)
    line([BasePairedResidues(i) BasePairedResidues(i)]+0.4,[-320 -319]+LastPos+FirstPos+1,'Color','r','LineStyle','-','LineWidth',4)
    line([borderLOW borderLOW+1],LastPos-[BasePairedResidues(i) BasePairedResidues(i)]-0.5+FirstPos,'Color','w','LineStyle','-','LineWidth',4)
    line([borderUP borderUP-1],LastPos-[BasePairedResidues(i) BasePairedResidues(i)]-0.5+FirstPos,'Color','w','LineStyle','-','LineWidth',4)
end
for i = 1: NrStems
    idx_j =ResiduesI.(strcat('S',num2str(i)))+FirstPos-1; 
    idx_i = ResiduesJ.(strcat('S',num2str(i)))+FirstPos-1;
    for j = 1:length(idx_j)
        if any(BasePairedResidues == idx_j(j)) || any(BasePairedResidues == idx_i(j))
           IDXMarker = [any(BasePairedResidues == idx_j(j)),any(BasePairedResidues == idx_i(j))];
           if IDXMarker(1) == 1
               IDX = idx_j(j);
           else
               IDX = idx_i(j);
           end
           if sum(IDX <= ColorR) > 0
                plot(idx_j(j)+0.3,LastPos-idx_i(j)+FirstPos-0.5,'ro','MarkerSize',14);%,'MarkerFaceColor','w')
                plot(idx_i(j)+0.3,LastPos-idx_j(j)+FirstPos-0.5,'ro','MarkerSize',14);%,'MarkerFaceColor','w')
           else
                plot(idx_j(j)+0.3,LastPos-idx_i(j)+FirstPos-0.5,'wo','MarkerSize',14);%,'MarkerFaceColor','w')
                plot(idx_i(j)+0.3,LastPos-idx_j(j)+FirstPos-0.5,'wo','MarkerSize',14);%,'MarkerFaceColor','w')
           end
        end
    end
end

title('Stem score','FontSize',24,'FontWeight','bold')
xlabel('first residue','FontSize',24,'FontWeight','bold')
ylabel('second residue','FontSize',24,'FontWeight','bold')
% set(gca,'XTick',(FirstPos:2:LastPos)+0.5,'XTickLabel',num2str((FirstPos:2:LastPos)'),'FontSize',18,'FontWeight','bold')
% set(gca,'YTick',(FirstPos:2:LastPos)-0.5,'YTickLabel',num2str((LastPos:-2:FirstPos)'),'FontSize',18,'FontWeight','bold')
set(gca,'XTick',(FirstPos:3:LastPos)+0.5,'XTickLabel',num2str((FirstPos:3:LastPos)'),'FontSize',18,'FontWeight','bold')
set(gca,'YTick',(FirstPos:4:LastPos)-0.5,'YTickLabel',num2str((LastPos:-4:FirstPos)'),'FontSize',18,'FontWeight','bold')
h = colorbar;
Ticks = get(h,'YTick');
Ticks = [1 Ticks(end)];
%maxim = max(max(FinalColorMap2Plot));
%vec = [1 maxim/2 maxim];
colorbar('YTick',Ticks,'YTickLabel',num2str([0; round((Ticks(end)./4.*10)./10)]),'FontSize',24,'FontWeight','bold');
xlim([borderLOW borderUP])
ylim([(LastPos-320)+FirstPos+1 LastPos-253+FirstPos+1])
print(10,'-depsc2','Stems1a.eps')
%axis([min(AllEnds_j)-1 max(AllStarts_j)+1 min(AllStarts_i)-1 max(AllEnds_i)+1]+FirstPos-1)
%axis([230 290 230 290])
%fett(8)
%

% Write Stems out -------------------------------------
Outfolder = './Stems/';
try,
    cd(Outfolder)
    cd('../')
catch
    mkdir(Outfolder);
end
Outfile = nan(max(AllEnds_i-AllStarts_i+1),length(AllStarts_i)*4);
for i = 1:length(AllStarts_i)
    LengthOfStem = AllEnds_i(i)-AllStarts_i(i)+1;
    %Positions 1
    Outfile(1:LengthOfStem,(i-1)*4+1) = (AllStarts_i(i):AllEnds_i(i))+FirstPos-1;
    %Positions 2
    Outfile(1:LengthOfStem,(i-1)*4+2) = (AllStarts_j(i):-1:AllEnds_j(i))+FirstPos-1;
    %Score
    Outfile(1:LengthOfStem,(i-1)*4+3) = AllScores(i); 
    %Individual Epistasis Value
    idx_i = (AllStarts_i(i):1:AllEnds_i(i))+FirstPos-1;
    idx_j = (AllStarts_j(i):-1:AllEnds_j(i))+FirstPos-1;
    for j = 1:LengthOfStem
        Idx = find(PotentialPairsWCWB(:,1) == idx_i(j) & PotentialPairsWCWB(:,2) == idx_j(j));
        if ~isempty(Idx)
            Epi = MaxMedianEpistasisWCWB(Idx);
            Outfile(j,(i-1)*4+4) = Epi;
        end
    end
end
Header = {'First Pos.';'Second Pos.';'Stem Score';...
        'median Epistasis'};

filename = strcat(Outfolder,'StemsBeforeResolution.csv');
try,
    delete(filename);
catch
end
       
fid = fopen(filename,'w');
for j = 1:length(AllStarts_i)
    for i = 1:length(Header)
        fprintf(fid,'%s,',char(Header(i)));
    end
end
fprintf(fid,'%s\n','');
fclose(fid);
%write Output values
dlmwrite(filename,Outfile,'-append')
%----------------------------------------------------


%Resolve Conflicts

Isfeasible = true(1,NrStems);
Conflicts = zeros(NrStems,NrStems);

for i = 1:NrStems-1
    res_i = ResiduesI.(strcat('S',num2str(i)));
    res_j = ResiduesJ.(strcat('S',num2str(i)));
    for j = i+1:NrStems
        res_i2 = ResiduesI.(strcat('S',num2str(j)));
        res_j2 = ResiduesJ.(strcat('S',num2str(j)));
        c1 = intersect(res_i,res_i2);
        c2 = intersect(res_j,res_j2);
        c12 = intersect(res_i,res_j2);
        c21 = intersect(res_j,res_i2);
        if ~isempty(c1) || ~isempty(c2) || ~isempty(c12) || ~isempty(c21) 
            Conflicts(i,j) = 1;
%             if AllScores(i) > AllScores(j)
%                Isfeasible(j) = false;
%             elseif AllScores(i) <= AllScores(j)
%                Isfeasible(i) = false;
%             end
        end
    end
end
Conflicts = Conflicts + Conflicts';
figure(20)
image(Conflicts*30)
xlabel('stem #')
ylabel('stem #')
title('conflicting stems')

[AllScores_sorted,Idx] = sort(AllScores,'descend');

for i = 1:NrStems
    stem = Idx(i);
    if AllScores(stem) > 0
        Conflicting = find(Conflicts(stem,:));
        if ~isempty(Conflicting)
            Isfeasible(Conflicting) = false;
            AllScores(Conflicting) = 0;
        end
    end
end

AllScores = AllScores(Isfeasible);
AllStarts_j = AllStarts_j(Isfeasible);
AllStarts_i = AllStarts_i(Isfeasible);
AllEnds_j = AllEnds_j(Isfeasible);
AllEnds_i = AllEnds_i(Isfeasible);
%IsWatsonCrick = IsWatsonCrick(Isfeasible);

NrStems = length(AllScores);
FinalColorMap3 = zeros(NrPositionsTotal); % for illustration

for i = 1: NrStems
    LengthOfStem = AllEnds_i(i) - AllStarts_i(i) + 1;
    idx_i = AllStarts_i(i):1:AllEnds_i(i);
    idx_j = AllStarts_j(i):-1:AllEnds_j(i);
    %if IsWatsonCrick(i)
        for j = 1:LengthOfStem 
            FinalColorMap3(idx_i(j),idx_j(j)) = AllScores(i);
        end
    %end
end

FinalColorMap3Plot = FinalColorMap3'+FinalColorMap3;
figure(9)
image(xvalues,fliplr(xvalues),FinalColorMap3Plot.*4)
colormap gray
hold on
plot([FirstPos LastPos]-0.5,fliplr([FirstPos LastPos])-0.5,'-w','LineWidth',2)
title('Stem score','FontSize',24,'FontWeight','bold')
xlabel('first residue','FontSize',24,'FontWeight','bold')
ylabel('second residue','FontSize',24,'FontWeight','bold')
set(gca,'XTick',(FirstPos:20:LastPos)+0.5,'XTickLabel',num2str((FirstPos:20:LastPos)'),'FontSize',18,'FontWeight','bold')
set(gca,'YTick',(FirstPos:20:LastPos)-0.5,'YTickLabel',num2str((LastPos:-20:FirstPos)'),'FontSize',18,'FontWeight','bold')
h = colorbar;
Ticks = get(h,'YTick');
Ticks = [1 Ticks(end)];
%maxim = max(max(FinalColorMap3Plot));
%vec = [1 maxim/2 maxim];
colorbar('YTick',Ticks,'YTickLabel',num2str([0; round((Ticks(end)./4.*10)./10)]),'FontSize',24,'FontWeight','bold');
%axis([min(AllEnds_j)-1 max(AllStarts_j)+1 min(AllStarts_i)-1 max(AllEnds_i)+1]+FirstPos-1)
%axis([31 505 31 505])
%fett(9)
print(9,'-depsc2','Stems2.eps')
%-----------------All conflicts have been resolved-------------

%%
% Write Stems out
% Outfolder = './Stems/';
% try,
%     cd(Outfolder)
%     cd('../')
% catch
%     mkdir(Outfolder);
% end
Outfile = nan(max(AllEnds_i-AllStarts_i+1),length(AllStarts_i)*4);
for i = 1:length(AllStarts_i)
    LengthOfStem = AllEnds_i(i)-AllStarts_i(i)+1;
    %Positions 1
    Outfile(1:LengthOfStem,(i-1)*4+1) = (AllStarts_i(i):AllEnds_i(i))+FirstPos-1;
    %Positions 2
    Outfile(1:LengthOfStem,(i-1)*4+2) = (AllStarts_j(i):-1:AllEnds_j(i))+FirstPos-1;
    %Score
    Outfile(1:LengthOfStem,(i-1)*4+3) = AllScores(i); 
    %Individual Epistasis Value
    idx_i = (AllStarts_i(i):1:AllEnds_i(i))+FirstPos-1;
    idx_j = (AllStarts_j(i):-1:AllEnds_j(i))+FirstPos-1;
    for j = 1:LengthOfStem
        Idx = find(PotentialPairsWCWB(:,1) == idx_i(j) & PotentialPairsWCWB(:,2) == idx_j(j));
        if ~isempty(Idx)
            Epi = MaxMedianEpistasisWCWB(Idx);
            Outfile(j,(i-1)*4+4) = Epi;
        end
    end
end
% Header = {'First Pos.';'Second Pos.';'Stem Score';...
%         'median Epistasis'};

filename = strcat(Outfolder,'Stems.csv');
try,
    delete(filename);
catch
end
       
fid = fopen(filename,'w');
for j = 1:length(AllStarts_i)
    for i = 1:length(Header)
        fprintf(fid,'%s,',char(Header(i)));
    end
end
fprintf(fid,'%s\n','');
fclose(fid);
%write Output values
dlmwrite(filename,Outfile,'-append')
