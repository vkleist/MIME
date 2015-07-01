%function [] = SingleVarAnalysis_BasedOnCoVarData08_07_PlotKDs()

close all
clear all

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
M = dlmread('InputVariables.txt');


MinimumSignal2NoiseStrength = M(1);
alpha = M(2);
MinimumNumberEstimatableKDs = M(3);
GenerateOutput = M(4);
WeightThreshold = M(5);
FullPlot = M(6);
PlotStartRegion = M(7);
PlotEndRegion = M(8);
PutMedian = M(9);%for inputation
SummaryPlot = M(10);% 0 = plot mean effect over all mutations away from wt; 1 = plot effect of mutation that maximally affects Kd
clear M

FirstXTick = 0;
XTickIncrement= 100;

%--------------------------

%NrPositionsTotal = 535;
bases = 1:4;
%For re-ordering entries in the case where the first-and second base were
%swapt
%-----------------------------

%Read Reference sequence
disp('Read RefSeq..')
Tmp = dlmread('./Data_CoVariation/RefSeq.txt');
RefSeq = Tmp(:,2);

%----------------------------
%filename = strcat('./ResultsSN',num2str(MinimumSignal2NoiseStrength),'W',num2str(10*WeightThreshold),'.mat');

filename = './Results.mat';
load(filename);
%load 'ResultsSN1W5.mat';
%SingleVarAnalysis_BasedOnCoVarData07_07_ComputeKDs;

%% post-processing based on quality criteria

% a) Minimum Signal-2-Noise & b) Lower estimates correction
%m_w
LogicIdx = Signal2Noise_m_w_SuperN < MinimumSignal2NoiseStrength & Signal2Noise_m_w_Beads < MinimumSignal2NoiseStrength;
for outercounter = CutValueFwd+1:NrPositionsTotal-CutValueBwd
    for j = 1:4
        idx = LogicIdx(outercounter,:,j);
        TotalRel_KD_m_w(outercounter,idx,j) = nan;
    end
end

LowerLimitsKD_m_w = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4); 
LogicIdx = Signal2Noise_m_w_SuperN >= MinimumSignal2NoiseStrength & Signal2Noise_m_w_Beads < MinimumSignal2NoiseStrength;
LowerLimitsKD_m_w(LogicIdx) = 1;

UpperLimitsKD_m_w = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)),4); 
LogicIdx = Signal2Noise_m_w_SuperN < MinimumSignal2NoiseStrength & Signal2Noise_m_w_Beads >= MinimumSignal2NoiseStrength;
UpperLimitsKD_m_w(LogicIdx) = 1;

%M_w
LogicIdx = Signal2Noise_M_w_SuperN < MinimumSignal2NoiseStrength & Signal2Noise_M_w_Beads < MinimumSignal2NoiseStrength;
for outercounter = CutValueFwd+1:NrPositionsTotal-CutValueBwd
   idx = LogicIdx(outercounter,:);
   TotalRel_KD_M_w(outercounter,idx) = nan;
end
LowerLimitsKD_M_w = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd)));
LogicIdx = Signal2Noise_M_w_SuperN >= MinimumSignal2NoiseStrength & Signal2Noise_M_w_Beads < MinimumSignal2NoiseStrength;
LowerLimitsKD_M_w(LogicIdx) = 1;

UpperLimitsKD_M_w = zeros(NrPositionsTotal,NrComparisons*(NrPositionsTotal-(CutValueFwd+CutValueBwd))); 
LogicIdx = Signal2Noise_M_w_SuperN < MinimumSignal2NoiseStrength & Signal2Noise_M_w_Beads >= MinimumSignal2NoiseStrength;
UpperLimitsKD_M_w(LogicIdx) = 1;

%
%
% c) censoring by minimum coverage
CensoredPositionsBeadsTmp = PositionWeightsBeads<WeightThreshold;
CensoredPositionsSuperNTmp = PositionWeightsSuperN<WeightThreshold;
CensoredPositionsTotalTmp = PositionWeightsTotal<WeightThreshold;

CensoredPositionsBeads = false(NrPositionsTotal,(NrPositionsTotal-(CutValueFwd+CutValueBwd))*NrComparisons);
CensoredPositionsSuperN = false(NrPositionsTotal,(NrPositionsTotal-(CutValueFwd+CutValueBwd))*NrComparisons);
CensoredPositionsTotal = false(NrPositionsTotal,(NrPositionsTotal-(CutValueFwd+CutValueBwd))*NrComparisons);
    

for counter = 1:NrComparisons 
    idx = counter:NrComparisons:(NrPositionsTotal-(CutValueFwd+CutValueBwd))*NrComparisons;
    CensoredPositionsBeads(:,idx) = CensoredPositionsBeadsTmp(:,CutValueFwd+1:NrPositionsTotal-CutValueBwd,counter);
    CensoredPositionsSuperN(:,idx) = CensoredPositionsSuperNTmp(:,CutValueFwd+1:NrPositionsTotal-CutValueBwd,counter);
    CensoredPositionsTotal(:,idx) = CensoredPositionsTotalTmp(:,CutValueFwd+1:NrPositionsTotal-CutValueBwd,counter);
end
%
% CensoredPositionsSum = CensoredPositionsBeads+CensoredPositionsSuperN+CensoredPositionsTotal;
% CensoredPositionsFinal = CensoredPositionsSum>0;
CensoredPositionsFinal = max(max(CensoredPositionsBeads,CensoredPositionsSuperN),CensoredPositionsTotal);
TotalRel_KD_M_w(CensoredPositionsFinal) = nan;
LowerLimitsKD_M_w(CensoredPositionsFinal) = 0;
UpperLimitsKD_M_w(CensoredPositionsFinal) = 0;

for outercounter = CutValueFwd+1:NrPositionsTotal-CutValueBwd
    for j = 1:4
        idx = CensoredPositionsFinal(outercounter,:);
        TotalRel_KD_m_w(outercounter,idx,j) = nan;
        LowerLimitsKD_m_w(outercounter,idx,j) = 0;
        UpperLimitsKD_m_w(outercounter,idx,j) = 0;
    end
end
%
% d) replace arbitrary set lower bounds for KD with realistic values.

for outercounter = CutValueFwd+1:NrPositionsTotal-CutValueBwd
    for j = 1:4
        Positions_m_w = TotalRel_KD_m_w(outercounter,:,j) == -1 | TotalRel_KD_m_w(outercounter,:,j) == -10;%lower and upper estimates
        minmax_m_w = prctile(TotalRel_KD_m_w(outercounter,~Positions_m_w,j),[5 50 95]);
        Positions_m_w_Low = TotalRel_KD_m_w(outercounter,:,j) == -1;
        Positions_m_w_Up = TotalRel_KD_m_w(outercounter,:,j) == -10;
        
        if PutMedian
            TotalRel_KD_m_w(outercounter,Positions_m_w_Low,j) = max(minmax_m_w(2),0);
            TotalRel_KD_m_w(outercounter,Positions_m_w_Up,j) = max(minmax_m_w(2),0);
        else
            TotalRel_KD_m_w(outercounter,Positions_m_w_Low,j) = max(minmax_m_w(3),0);
            TotalRel_KD_m_w(outercounter,Positions_m_w_Up,j) = max(minmax_m_w(1),0);
        end 
    end
    Positions_M_w = TotalRel_KD_M_w(outercounter,:) == -1 | TotalRel_KD_M_w(outercounter,:) == -10;
    minmax_M_w = prctile(TotalRel_KD_M_w(outercounter,~Positions_M_w),[5 50 95]);
    Positions_M_w_Low = TotalRel_KD_M_w(outercounter,:) == -1;
    Positions_M_w_Up = TotalRel_KD_M_w(outercounter,:) == -10;
    
    if PutMedian
        TotalRel_KD_M_w(outercounter,Positions_M_w_Low) = max(minmax_M_w(2),0);
        TotalRel_KD_M_w(outercounter,Positions_M_w_Up) = max(minmax_M_w(2),0);
    else
        TotalRel_KD_M_w(outercounter,Positions_M_w_Low) = max(minmax_M_w(3),0);
        TotalRel_KD_M_w(outercounter,Positions_M_w_Up) = max(minmax_M_w(1),0);
    end
end

%-------------------------------
%Generate Output

TotalNumberPValuesComputed = 0;
RawPvalues = nan(NrPositionsTotal,4);
OutputMatrix = nan(NrPositionsTotal-(CutValueFwd+CutValueBwd),31);
OutputMatrix(:,1) = CutValueFwd+1:NrPositionsTotal-CutValueBwd;%Position
OutputMatrix(:,2) = RefSeq(CutValueFwd+1:NrPositionsTotal-CutValueBwd);%wt

OutputMatrix2 = nan(NrPositionsTotal-(CutValueFwd+CutValueBwd),10);
OutputMatrix2(:,1:2)=OutputMatrix(:,1:2);

%Go over all mutant
%-> A; -> C; -> G; -> U
for j = 1:4 
    TMP = nanmedian(TotalRel_KD_m_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:,j),2); 
    TMP(isnan(TMP)) = 1;
    OutputMatrix(:,(j-1)*7 + 4) = TMP;
    TMP_5_95 = prctile(TotalRel_KD_m_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:,j),[5 95],2);
    TMP_5_95(isnan(TMP_5_95)) = 1;
    OutputMatrix(:,(j-1)*7 + [9:10]) = TMP_5_95;

    LowerEstimates = LowerLimitsKD_m_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:,j);
    ResWithLowerEstimates = sum(LowerEstimates,2);
    OutputMatrix(:,(j-1)*7+7) = ResWithLowerEstimates;

    UpperEstimates = UpperLimitsKD_m_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:,j);
    ResWithUpperEstimates = sum(UpperEstimates,2);
    OutputMatrix(:,(j-1)*7+8) = ResWithUpperEstimates;
    
    %Positions where the wildtype is 'j'
    Positions = find(RefSeq == j);
    %Corresponding possible mutations
    SubPlotIdx = bases(1:4~=j);
    for q = SubPlotIdx %Go over MutantBases (3 in total)        
        Mut = log2(TotalRel_KD_m_w(Positions,:,q)');
        % test for significant increase of KD (p-values)
        NumberCalls = sum(~isnan(Mut));%number of resamplings
        NumberCallsGreaterthreshold = NumberCalls>=MinimumNumberEstimatableKDs;%eligible resamplings

        Pvalue = min((sum(Mut<=0)./sum(~isnan(Mut))),(sum(Mut>=0)./sum(~isnan(Mut))));%raw p-value
        Pvalue(~NumberCallsGreaterthreshold) = nan;%eligible p-values

        NrComparisons = sum(~isnan(Pvalue));
        TotalNumberPValuesComputed = TotalNumberPValuesComputed+NrComparisons;%counter total nr. of elegible p-value calculations
        RawPvalues(Positions,q) = Pvalue; % raw P-Value

        TMP_PositionsLogic = (Positions > CutValueFwd & Positions < NrPositionsTotal - CutValueBwd);
        TMP_Positions = Positions(TMP_PositionsLogic);
        OutputMatrix(TMP_Positions-CutValueFwd,(q-1)*7+6) = NumberCalls(TMP_PositionsLogic)';% assign nr. of resamplings
    end
end

% P-value correction for multiple testing
PvaluesTMP = RawPvalues(:);
[Pvalues_Sorted,I] = sort(PvaluesTMP);

%Benjamini Hochberg False Discovery Rate
for i = (TotalNumberPValuesComputed:-1:1)
    candidate = Pvalues_Sorted(i)*TotalNumberPValuesComputed/i;
    Pvalues_Sorted(i) = candidate;
end
PvaluesTMP(I) = Pvalues_Sorted;
Pvalues = reshape(PvaluesTMP,NrPositionsTotal,4);
Positions = 1:NrPositionsTotal;
TMP_PositionsLogic = (Positions > CutValueFwd & Positions < NrPositionsTotal - CutValueBwd+1);
TMP_Positions = Positions(TMP_PositionsLogic);
%assign corrected p-value
OutputMatrix(TMP_Positions-CutValueFwd,((1:4)-1)*7+5) = Pvalues(TMP_PositionsLogic,:);

%%% find mutation that has the maximal effect on KD mmax
% maximum effect for significant positions only
indexmmax = ones(length(TMP_Positions),1);
for i = TMP_Positions-CutValueFwd
    if any(OutputMatrix(i,((1:4)-1)*7+5) < alpha)
        idx = find(OutputMatrix(i,((1:4)-1)*7+5) < alpha);
    else
        idx = 1:4;
    end
    [~,k] = max(abs(log2(OutputMatrix(i,((idx)-1)*7 + 4)))); 
    indexmmax(i) = idx(k);
end
% %maximum effect
%[~,indexmmax] = max(abs(log2(OutputMatrix(:,((1:4)-1)*7 + 4))),[],2);   
OutputMatrix(:,3) = indexmmax;
OutputMatrix2(:,3) = OutputMatrix(:,3);

if SummaryPlot == 0
%%% Option 1; Depict & write as output the effect of any mutation away from wt
    %-> Any mutated
    % median, percentiles
    TMP = nanmedian(TotalRel_KD_M_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:),2);%-> Any mutated
    TMP(isnan(TMP)) = 1;
    OutputMatrix2(:,4) =  TMP;
    TMP_5_95 = prctile(TotalRel_KD_M_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:),[5 95],2);
    TMP_5_95(isnan(TMP_5_95)) = 1;
    OutputMatrix2(:,9:10) = TMP_5_95;
    % number of resamplings that are lower or upper estimates
    LowerEstimates = LowerLimitsKD_M_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:);
    ResWithLowerEstimates = sum(LowerEstimates,2);
    OutputMatrix2(:,7) = ResWithLowerEstimates;

    UpperEstimates = UpperLimitsKD_M_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:);
    ResWithUpperEstimates = sum(UpperEstimates,2);
    OutputMatrix2(:,8) = ResWithUpperEstimates;
    % number of resamplings in total
    NumberCalls = sum(~isnan(TotalRel_KD_M_w(CutValueFwd+1:NrPositionsTotal-CutValueBwd,:)'))';
    OutputMatrix2(:,6) = NumberCalls;
    %p-value
    Mut = log2(TotalRel_KD_M_w');
    % test for significant increase of KD
    NumberCallsGreaterthreshold = NumberCalls>=MinimumNumberEstimatableKDs;

    Pvalue_M_w = min((sum(Mut<=0)./sum(~isnan(Mut))),(sum(Mut>=0)./sum(~isnan(Mut))));
    Pvalue_M_w(~NumberCallsGreaterthreshold) = nan;
    NrComparisons = sum(~isnan(Pvalue_M_w));

    % P-value correction for multiple testing
    [Pvalue_M_w_Sorted,I] = sort(Pvalue_M_w);
    %Benjamini Hochberg False Discovery Rate
    for i = (NrComparisons:-1:1)
        candidate = Pvalue_M_w_Sorted(i)*NrComparisons/i;
        Pvalue_M_w_Sorted(i) = candidate;
    end
    Pvalue_M_w(I) = Pvalue_M_w_Sorted;
    TMP_PositionsLogic = (1:NrPositionsTotal > CutValueFwd & 1:NrPositionsTotal < NrPositionsTotal - CutValueBwd+1);
    OutputMatrix2(:,5) =  Pvalue_M_w(TMP_PositionsLogic)';

elseif SummaryPlot == 1
%%% Option 2: Depict & write as output the mutation with the maximum effect on Kd.
    for i = 1:length(indexmmax)
        OutputMatrix2(i,4) = OutputMatrix(i,(indexmmax(i)-1)*7 + 4);% median
        OutputMatrix2(i,9:10) = OutputMatrix(i,(indexmmax(i)-1)*7 + [9:10]);% percentiles
        OutputMatrix2(i,7) = OutputMatrix(i,(indexmmax(i)-1)*7+7);% # lower estimates
        OutputMatrix2(i,8) = OutputMatrix(i,(indexmmax(i)-1)*7+8);% # upper estimates
        OutputMatrix2(i,6) = OutputMatrix(i,(indexmmax(i)-1)*7+6);% total # estimates
        OutputMatrix2(i,5) = OutputMatrix(i,(indexmmax(i)-1)*7+5);% p-value
    end    
else
    disp('Change Option for summary plot: "0": average effect of mutation wt -> any mut.; "1": effect of mutation that maximally impacts Kd')
end


%% Statistic Test based on Re-sampling

TitleLabel = ['A','C','G','U'];

BoxPlotColors = [1 0 0; % red
                 0 1 0; % green
                 0 0 1; % blue
                 1 0.7 0.1];%yellow
Region = PlotStartRegion:PlotEndRegion;
InterestingRefSeq = RefSeq(Region);
[n,m,z] = size(TotalRel_KD_m_w);
%Plot individual mutation Effects

YMarkSignificant = 5.5;
for j = 1:4 %wt Base
    SubPlotIdx = bases(1:4~=j);
    Positions = find(InterestingRefSeq == j)+PlotStartRegion-1;
    counter = 1;
    if ~isempty(Positions)
        
        for q = SubPlotIdx %Go over MutantBases (3 in total)
            Mut2 = log2(TotalRel_KD_m_w(Positions,:,q)');
            LowerEstimates2 = LowerLimitsKD_m_w(Positions,:,q)';
            UpperEstimates2 = UpperLimitsKD_m_w(Positions,:,q)';
            
            if FullPlot
                figure(j)
                hold on
                boxplot(Mut2,'plotstyle','compact','colors',BoxPlotColors(q,:),'positions',[1:length(Positions)]+(counter-2)/5);
                set(gca,'XTickLabel',{' '})
            end
 
            SignificantlyGreaterThanZero2 = nan(1,length(Positions));
            ResWithLowerEstimates2 = nan(1,length(Positions));
            NumberCallsGreaterthreshold2 = sum(~isnan(Mut2))>=MinimumNumberEstimatableKDs;
            Pvalue2 = Pvalues(Positions,q)';
            
            SignificantlyGreaterThanZeroTmp = all([(Pvalue2 < alpha);(nanmedian(Mut2) > 0)]); % 
            SignificantlyGreaterThanZero2(all([SignificantlyGreaterThanZeroTmp;NumberCallsGreaterthreshold2])) = YMarkSignificant;
            ResWithLowerEstimatesTmp = sum(LowerEstimates2);
            ResWithLowerEstimates2(all([SignificantlyGreaterThanZeroTmp;NumberCallsGreaterthreshold2;ResWithLowerEstimatesTmp])) = 1; 
            
            if FullPlot
                figure(j)
                plot([1:length(Positions)]+(counter-2)/5,SignificantlyGreaterThanZero2+(counter-2)/20,'*','Color',BoxPlotColors(q,:),'MarkerSize',12)
                plot([1:length(Positions)]+(counter-2)/5,ResWithLowerEstimates2*YMarkSignificant+(counter-2)/20,'s','Color',BoxPlotColors(q,:),'MarkerSize',12)
            end

            SignificantlySmallerThanZero2 = nan(1,length(Positions));
            ResWithGreaterEstimates2 = nan(1,length(Positions));
            
            SignificantlySmallerThanZeroTmp = all([(Pvalue2 < alpha);(nanmedian(Mut2) < 0)]);
            SignificantlySmallerThanZero2(all([SignificantlySmallerThanZeroTmp;NumberCallsGreaterthreshold2])) = -2;
            ResWithGreaterEstimatesTmp = sum(UpperEstimates2);
            ResWithGreaterEstimates2(all([SignificantlySmallerThanZeroTmp;NumberCallsGreaterthreshold2;ResWithGreaterEstimatesTmp])) = 1;
             
            if FullPlot
                figure(j)
                plot([1:length(Positions)]+(counter-2)/5,SignificantlySmallerThanZero2+(counter-2)/20,'+','Color',BoxPlotColors(q,:),'MarkerSize',12)
                plot([1:length(Positions)]+(counter-2)/5,ResWithGreaterEstimates2*-2+(counter-2)/20,'o','Color',BoxPlotColors(q,:),'MarkerSize',12)

            end
            
            clear Mut LowerEstimates  NumberCallsGreaterthreshold Pvalue SignificantlyGreaterThanZero ResWithLowerEstimates SignificantlySmallerThanZero
            counter = counter +1;
        end
        if FullPlot
            figure(j)
            hold on
            line([1 length(Positions)], [0 0],'Color','k','LineStyle','--')
            axis([0.5 length(Positions)+0.5 -3 7])
            set(gca,'XTick',1:length(Positions),'XTickLabel',num2str(Positions(1:end)))
            xlabel('position')
            title(strcat('Wt Base: ', TitleLabel(j), '; Relative binding affinity ->A (red), ->C (green), ->G (blue) and ->U (yellow)'))
            ylabel('log_2(KD_{m}/KD_{w})')
            set(get(gca,'xlabel'),'Fontsize',12);
            set(get(gca,'ylabel'),'Fontsize',12);
            set(get(gca,'title'),'Fontsize',12);
        end
    end

end

%% Tome-like plot
RefSeqAlphabetic = strrep(strrep(strrep(strrep(strrep(num2str(RefSeq'),'1','A'),'2','C'),'3','G'),'4','U'),' ','');
for j = 1:4 %wt Base
    Positions = find(InterestingRefSeq == j)+PlotStartRegion-1;
    SubPlotIdx = bases(1:4~=j);
    counter = 1;
    if ~isempty(Positions)
        for q= SubPlotIdx
            Mut = nan(m,length(Region));
            Mut(:,Positions-PlotStartRegion+1) = log2(TotalRel_KD_m_w(Positions,:,q)');
            ErrorBarData = prctile(Mut,[25 50 75]);
            if length(Region) == 1
                ErrorBarData = ErrorBarData';
            end
            
            figure(10)
            hold on
            if q == 4
                errorbar([Region]+(counter-2)/5, ErrorBarData(2,:),ErrorBarData(2,:)-ErrorBarData(1,:),ErrorBarData(2,:)-ErrorBarData(3,:),'LineStyle','none','Marker','o','Color',BoxPlotColors(q,:),'MarkerFaceColor',BoxPlotColors(q,:),'MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2);
            else
                errorbar([Region]+(counter-2)/5, ErrorBarData(2,:),ErrorBarData(2,:)-ErrorBarData(1,:),ErrorBarData(2,:)-ErrorBarData(3,:),'LineStyle','none','Marker','o','Color',BoxPlotColors(q,:),'MarkerFaceColor',BoxPlotColors(q,:),'MarkerSize',10,'LineWidth',2);
            end
            set(gca,'XTickLabel',{' '})
            
            figure(11)
            hold on
            boxplot(Mut,'plotstyle','compact','colors',BoxPlotColors(q,:),'positions',Region+(counter-2)/4,'symbol','');
            set(gca,'XTickLabel',{' '})
            
            figure(12)
            hold on
            IDX = any(~isnan(Mut));
            if sum(IDX)>0
                distributionPlot(Mut,'showMM',0,'xValues',Region+(counter-2)/4,'color',BoxPlotColors(q,:),'distWidth',0.4,'globalNorm',0)%compares shapes; if 3 then all have same area.
                plot(Region+(counter-2)/4,nanmedian(Mut),'ko','MarkerFaceColor','w');%,'LineWidth',0.5)
                plot(Region+(counter-2)/4,nanmedian(Mut),'k.','LineWidth',2)
            end
            
            counter = counter+1;
            clear Mut ErrorBarData   
        end
        %
        for i = 1:length(Positions)
            if j == 4
                figure(10)
                text(Positions(i),-2,RefSeqAlphabetic(Positions(i)),'Color',[1 0.7 0.1],'FontSize',14,'FontWeight','bold','HorizontalAlignment','center')
                figure(11)
                text(Positions(i),-2,RefSeqAlphabetic(Positions(i)),'Color',[1 0.7 0.1],'FontSize',14,'FontWeight','bold','HorizontalAlignment','center')
                figure(12)
                text(Positions(i),-2,RefSeqAlphabetic(Positions(i)),'Color',[1 0.7 0.1],'FontSize',14,'FontWeight','bold','HorizontalAlignment','center')
            else
                figure(10)
                text(Positions(i),-2,RefSeqAlphabetic(Positions(i)),'Color',BoxPlotColors(j,:),'FontSize',14,'FontWeight','bold','HorizontalAlignment','center')
                figure(11)
                text(Positions(i),-2,RefSeqAlphabetic(Positions(i)),'Color',BoxPlotColors(j,:),'FontSize',14,'FontWeight','bold','HorizontalAlignment','center')
                figure(12)
                text(Positions(i),-2,RefSeqAlphabetic(Positions(i)),'Color',BoxPlotColors(j,:),'FontSize',14,'FontWeight','bold','HorizontalAlignment','center')
            end
        end
        counter = counter+1; 
    end
end
    

figure(10)
hold on
line([PlotStartRegion-0.5 PlotEndRegion+0.5], [0 0],'Color','k','LineStyle','--')
if length(Region) == 1
    xlim([PlotStartRegion-3 PlotEndRegion+3]);
else
    xlim([PlotStartRegion-1 PlotEndRegion+1])
end
ylim([-2.5 6])
%axis([PlotStartRegion-1 PlotEndRegion+1 -4 8])
Xtickers = PlotStartRegion:max(round((PlotEndRegion-PlotStartRegion)/10),1):PlotEndRegion;
set(gca,'XTick',Xtickers,'XTickLabel',num2str(Xtickers'))
xlabel('position')
%title(strcat('Wt Base: ', TitleLabel(j), '; Relative binding affinity ->A (blue), ->C (green), ->G (red) and ->U (magenta)'))
ylabel('log_2(KD_{m}/KD_{w})')
set(get(gca,'xlabel'),'Fontsize',12);
set(get(gca,'ylabel'),'Fontsize',12);
set(get(gca,'title'),'Fontsize',12);

figure(11)
hold on
line([PlotStartRegion-0.5 PlotEndRegion+0.5], [0 0],'Color','k','LineStyle','--')
if length(Region) == 1
    xlim([PlotStartRegion-2 PlotEndRegion+2]);
else
    xlim([PlotStartRegion-0.5 PlotEndRegion+0.5])
end
ylim([-3.5 8])
set(gca,'XTick',Xtickers,'XTickLabel',num2str(Xtickers'))
%set(gca,'XTick',FirstXTick:XTickIncrement:NrPositionsTotal,'XTickLabel',num2str([1:XTickIncrement:NrPositionsTotal]'))
xlabel('position')
title(strcat('Relative binding affinity wt -> mutant'))
ylabel('log_2(KD_{m}/KD_{w})')
set(get(gca,'xlabel'),'Fontsize',12);
set(get(gca,'ylabel'),'Fontsize',12);
set(get(gca,'title'),'Fontsize',12);

figure(12)
hold on
line([PlotStartRegion-0.5 PlotEndRegion+0.5], [0 0],'Color','k','LineStyle','--')
if length(Region) == 1
    xlim([PlotStartRegion-3 PlotEndRegion+3]);
else
    xlim([PlotStartRegion-1 PlotEndRegion+1])
end
ylim([-2.5 6])
set(gca,'XTick',Xtickers,'XTickLabel',num2str(Xtickers'))
xlabel('position')
title(strcat('Relative binding affinity wt -> mutant'))
ylabel('log_2(KD_{m}/KD_{w})')
set(get(gca,'xlabel'),'Fontsize',12);
set(get(gca,'ylabel'),'Fontsize',12);
set(get(gca,'title'),'Fontsize',12);

%%
%Plot Effect of any mutation away from wt.
YMarkSignificant = 3;
Positions = 1:NrPositionsTotal;

if SummaryPlot == 0
    %%% Option 1: "average" effect of any mutation away from wt.
    Mut = log2(TotalRel_KD_M_w(Positions,:)');
    LowerEstimates = LowerLimitsKD_M_w(Positions,:)';
    Pvalue =  Pvalue_M_w;
    figtitle = ' wt -> any mutant';
    ylabelstring = 'M';
elseif SummaryPlot == 1
    %%% Option 2: Effect of mutation that maximally affects Kd.
    [n,m,z] = size(TotalRel_KD_m_w);
    TotalRel_KD_mmax_w = nan(n,m);
    LowerLimitsKD_mmax_w = nan(n,m);
    Pvalue = nan(1,n);
    IndexmMax = [ones(CutValueFwd,1);indexmmax;ones(CutValueBwd,1)];
    for i = Positions
        Pvalue(i) = Pvalues(i,IndexmMax(i));
        TotalRel_KD_mmax_w(i,:) = TotalRel_KD_m_w(i,:,IndexmMax(i));
        LowerLimitsKD_mmax_w(i,:) = LowerLimitsKD_m_w(i,:,IndexmMax(i));
    end
    Mut = log2(TotalRel_KD_mmax_w(Positions,:)');
    LowerEstimates = LowerLimitsKD_mmax_w(Positions,:)';
    figtitle = ' wt -> m_{max}';
    ylabelstring = 'm_{max}';
else
    disp('Change Option for summary plot: "0": average effect of mutation wt -> any mut.; "1": effect of mutation that maximally affects Kd')
end
    
if FullPlot
    figure(5)
    hold on
    boxplot(Mut,'plotstyle','compact','colors','k','positions',Positions);
    set(gca,'XTickLabel',{' '})
end

SignificantlyGreaterThanZeroTmp = Pvalue < alpha & nanmedian(Mut)>0; % 
SignificantlyGreaterThanZero = nan(1,length(Positions));
SignificantlyGreaterThanZero(SignificantlyGreaterThanZeroTmp) = YMarkSignificant;

ResWithLowerEstimatesTmp = sum(LowerEstimates);
ResWithLowerEstimates = nan(1,length(Positions));
ResWithLowerEstimates(all([SignificantlyGreaterThanZeroTmp;ResWithLowerEstimatesTmp])) = 1;


if FullPlot
    plot(Positions,SignificantlyGreaterThanZero,'*','Color','k','MarkerSize',12)
    plot(Positions,ResWithLowerEstimates.*YMarkSignificant,'s','Color','k','MarkerSize',12)
end

SignificantlySmallerThanZeroTmp = Pvalue < alpha & nanmedian(Mut)<0; % 
SignificantlySmallerThanZero = nan(1,length(Positions));
SignificantlySmallerThanZero(SignificantlySmallerThanZeroTmp) = -1.5;

%plot those ones where the KD estimate is a lower estimate
if FullPlot
    plot(Positions,SignificantlySmallerThanZero,'*','Color','k','MarkerSize',12)
end

if FullPlot
    figure(5)
    hold on
    line([1 length(Positions)], [0 0],'Color','k','LineStyle','--')
    axis([1 length(Positions) -3 7])
    set(gca,'XTick',FirstXTick:XTickIncrement:NrPositionsTotal,'XTickLabel',num2str([1:XTickIncrement:NrPositionsTotal]'))
    xlabel('position')
    title(strcat('Relative binding affinity',figtitle, ' (black)'))
    ylabel(strcat('log_2(KD_{',ylabelstring,'}/KD_{w})'))
    set(get(gca,'xlabel'),'Fontsize',12);
    set(get(gca,'ylabel'),'Fontsize',12);
    set(get(gca,'title'),'Fontsize',12);
end

figure(6)
hold on
XPositions = 1:NrPositionsTotal;
TMP_5_50_95 = prctile(Mut',[5 50 95 25 75],2);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(:,3),-3.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(:,5),-3.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(:,4),-3.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(:,1),-3.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated

plot(XPositions,TMP_5_50_95(:,2),'r','LineWidth',2)
plot(XPositions(SignificantlyGreaterThanZero==YMarkSignificant),YMarkSignificant,'kv')
plot(XPositions(SignificantlySmallerThanZero==-1.5),-1.5,'k^')
line([1 NrPositionsTotal],[0 0],'Color','k','LineStyle','--')
set(gca,'XTick',FirstXTick:XTickIncrement:NrPositionsTotal,'XTickLabel',num2str([FirstXTick:XTickIncrement:NrPositionsTotal]'))
set(gca,'ylim',[-2 4])
set(gca,'xlim',[-10 NrPositionsTotal+10])
xlabel('position')
title(strcat('Relative binding affinity',figtitle))
ylabel(strcat('log_2(KD_{',ylabelstring,'}/KD_{w})'))
set(get(gca,'xlabel'),'Fontsize',12);
set(get(gca,'ylabel'),'Fontsize',12);
set(get(gca,'title'),'Fontsize',12);
set(gca,'FontWeight','b');

%Smoothed figure
datafig2 = zeros(NrPositionsTotal,5);
windowSize = 5;
figure(7)
hold on
data = TMP_5_50_95(:,3);
SmoothedData = filter(ones(1,windowSize)/windowSize,1,data);
datafig2(:,3) = SmoothedData;
area(XPositions,SmoothedData,-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
data = TMP_5_50_95(:,5);
SmoothedData = filter(ones(1,windowSize)/windowSize,1,data);
datafig2(:,5) = SmoothedData;
area(XPositions,SmoothedData,-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
data = TMP_5_50_95(:,4);
SmoothedData = filter(ones(1,windowSize)/windowSize,1,data);
datafig2(:,4) = SmoothedData;
area(XPositions,SmoothedData,-1.9,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated

data = TMP_5_50_95(:,1);
SmoothedData = filter(ones(1,windowSize)/windowSize,1,data);
datafig2(:,1) = SmoothedData;
area(XPositions,SmoothedData,-1.9,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
data = TMP_5_50_95(:,2);
SmoothedData = filter(ones(1,windowSize)/windowSize,1,data);
datafig2(:,2) = SmoothedData;
plot(XPositions,SmoothedData,'r','LineWidth',2)
plot(XPositions(SignificantlyGreaterThanZero==YMarkSignificant),YMarkSignificant,'kv')
plot(XPositions(SignificantlySmallerThanZero==-1.5),-1.5,'k^')
line([1 NrPositionsTotal],[0 0],'Color','k','LineStyle','--')
set(gca,'XTick',FirstXTick:XTickIncrement:NrPositionsTotal,'XTickLabel',num2str([FirstXTick:XTickIncrement:NrPositionsTotal]'))
set(gca,'ylim',[-2 4])
set(gca,'xlim',[-10 NrPositionsTotal+10])
xlabel('position')
title(strcat('Relative binding affinity',figtitle,' (smoothed; window size =',num2str(windowSize),')'))
ylabel(strcat('log_2(KD_{',ylabelstring,'}/KD_{w})'))
set(get(gca,'xlabel'),'Fontsize',12);
set(get(gca,'ylabel'),'Fontsize',12);
set(get(gca,'title'),'Fontsize',12);
set(gca,'FontWeight','b');

dlmwrite('Figure2data_smoothed.txt',datafig2);
dlmwrite('Figure2data_raw.txt',[(1:NrPositionsTotal)',2.^TMP_5_50_95]);
%return
%%
figure(8)
hold on

idx = Pvalue>0.05 & nanmedian(Mut)> 0;
plot(XPositions(idx),Pvalue(idx),'k.')
idx = Pvalue>0.01 & Pvalue<0.05 & nanmedian(Mut)>0;
plot(XPositions(idx),Pvalue(idx),'ko','MarkerSize',5)
idx = Pvalue>0.001 & Pvalue<0.01 & nanmedian(Mut)>0;
plot(XPositions(idx),Pvalue(idx),'ko','MarkerSize',5,'MarkerFaceColor','k')
if any(Pvalue<0.001)
    idx = Pvalue<=0.001 & Pvalue>1e-4 & nanmedian(Mut)>0;
    plot(XPositions(idx),Pvalue(idx),'ko','MarkerSize',8,'MarkerFaceColor','k')
    if any(Pvalue<0.0001)
        idx = Pvalue <= 1e-4 & nanmedian(Mut)>0;
        plot(XPositions(idx),1e-4,'ko','MarkerSize',8,'MarkerFaceColor','k')
    end
end

line([XPositions(1) XPositions(end)], [0.05 0.05],'Color','k','LineStyle','--')
line([XPositions(1) XPositions(end)], [0.01 0.01],'Color','k','LineStyle','--')
line([XPositions(1) XPositions(end)], [0.001 0.001],'Color','k','LineStyle','--')

%plot(XPositions,Pvalue,'k','LineWidth',2)
set(gca,'YScale','log','ylim',[0.00008 1.1])
set(gca,'XTick',FirstXTick:XTickIncrement:NrPositionsTotal,'XTickLabel',num2str([FirstXTick:XTickIncrement:NrPositionsTotal]'),'YTick',[1e-4 1e-3 1e-2 1e-1 1], 'YTickLabel', {'< 1e-4'; '0.001'; '0.01'; '0.1'; '1'})
xlabel('position')
title(strcat('Significant increase on relative binding affinity',figtitle))
ylabel('BHFDR corrected p-value')
set(get(gca,'xlabel'),'Fontsize',12);
set(get(gca,'ylabel'),'Fontsize',12);
set(get(gca,'title'),'Fontsize',12);
set(gca,'FontWeight','b');
xlim([0 NrPositionsTotal+1]);

K = nanmedian(Mut);
idx = Pvalue<0.05 &K>0;
dlmwrite('IdentifiedPositionsFullSet_increasingKd.csv',[XPositions(idx)',K(idx)']);

K = nanmedian(Mut);
idx = Pvalue<0.05 &K<0;
dlmwrite('IdentifiedPositionsFullSet_decreasingKd.csv',[XPositions(idx)',K(idx)']);

%-----Zoom-------
ZoomRegionStart = PlotStartRegion;
ZoomRegionEnd = PlotEndRegion;
XPositionsZoom = ZoomRegionStart:ZoomRegionEnd;
PvaluesZoom = Pvalue(ZoomRegionStart:ZoomRegionEnd);
MutZoom = Mut(:,XPositionsZoom);
figure(9)
hold on
idx = PvaluesZoom>0.05 & nanmedian(MutZoom)>0;
plot(XPositionsZoom(idx),PvaluesZoom(idx),'k.')
idx = PvaluesZoom>0.01 & PvaluesZoom<0.05 & nanmedian(MutZoom)>0;
plot(XPositionsZoom(idx),PvaluesZoom(idx),'ko','MarkerSize',5)
Pos = find(idx)+ZoomRegionStart-1;
for i = 1:sum(idx)
    if rem(i,2)
        text(Pos(i),PvaluesZoom(Pos(i)-ZoomRegionStart+1),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color','r')
    else
        text(Pos(i),PvaluesZoom(Pos(i)-ZoomRegionStart+1),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'Color','r')
    end
end
idx = PvaluesZoom>0.001 & PvaluesZoom<0.01 & nanmedian(MutZoom)>0;
plot(XPositionsZoom(idx),PvaluesZoom(idx),'ko','MarkerSize',5,'MarkerFaceColor','k')
Pos = find(idx)+ZoomRegionStart-1;
for i = 1:sum(idx)
    %if rem(i,2)
        text(Pos(i),PvaluesZoom(Pos(i)-ZoomRegionStart+1),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color','r')
   % else
    %    text(Pos(i),PvaluesZoom(Pos(i)-ZoomRegionStart+1),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'Color','r')
    %end
end
idx = PvaluesZoom<=0.001 & PvaluesZoom>1e-4 & nanmedian(MutZoom)>0;
plot(XPositionsZoom(idx),PvaluesZoom(idx),'ko','MarkerSize',8,'MarkerFaceColor','k')
Pos = find(idx)+ZoomRegionStart-1;
for i = 1:sum(idx)
    if rem(i,2)
        text(Pos(i),PvaluesZoom(Pos(i)-ZoomRegionStart+1),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color','r')
    else
        text(Pos(i),PvaluesZoom(Pos(i)-ZoomRegionStart+1),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'Color','r')
    end
end
idx = PvaluesZoom <= 1e-4 & nanmedian(MutZoom)>0;
plot(XPositionsZoom(idx),ones(1,sum(idx)).*1e-4,'ko','MarkerSize',8,'MarkerFaceColor','k')
Pos = find(idx)+ZoomRegionStart-1;
for i = 1:sum(idx)
    if rem(i,2)
        text(Pos(i),1e-4*(1+(rand-0.5)/5),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color','r')
    else
        text(Pos(i),1e-4*(1+(rand-0.5)/5),num2str(Pos(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'Color','r')
    end
end

line([XPositionsZoom(1)-0.5 XPositionsZoom(end)+0.5], [0.05 0.05],'Color','k','LineStyle','--')
line([XPositionsZoom(1)-0.5 XPositionsZoom(end)+0.5], [0.01 0.01],'Color','k','LineStyle','--')
line([XPositionsZoom(1)-0.5 XPositionsZoom(end)+0.5], [0.001 0.001],'Color','k','LineStyle','--')

set(gca,'YScale','log','ylim',[0.00005 1.1])
Xtickers = max(round((ZoomRegionEnd-ZoomRegionStart)/11),1);
set(gca,'XTick',ZoomRegionStart:Xtickers:ZoomRegionEnd,'XTickLabel',num2str([ZoomRegionStart:Xtickers:ZoomRegionEnd]'),'YTick',[1e-4 1e-3 1e-2 1e-1 1], 'YTickLabel', {'< 1e-4'; '0.001'; '0.01'; '0.1'; '1'})
xlabel('position')
title(strcat('Significant increase on relative binding affinity',figtitle))
ylabel('BHFDR corrected p-value')
xlim([ZoomRegionStart-0.5 ZoomRegionEnd+0.5])
set(get(gca,'xlabel'),'Fontsize',12);
set(get(gca,'ylabel'),'Fontsize',12);
set(get(gca,'title'),'Fontsize',12);
set(gca,'FontWeight','b');
%%

if GenerateOutput
    %filename = strcat('./PositionWiseKDEstimatesSN',num2str(MinimumSignal2NoiseStrength),'W',num2str(10*WeightThreshold),'.csv');
    filename = strcat('./PositionWiseKDEstimates.csv');
    try,
        delete(filename);
    catch
    end
    %generate Header
    Header = {'Pos.';'wt Base';'mut. with max effect';' -> A (median)'; 'p-value';'# estimates';'# lower estim.';'# upper estim.';' -> A (5th prctile)';' -> A (95th prctile)';...
        ' -> C (median)'; 'p-value';'# estimates';'# lower estim.';'# upper estim.';' -> C (5th prctile)';' -> C (95th prctile)';...
        ' -> G (median)'; 'p-value';'# estimates';'# lower estim.';'# upper estim.';' -> G (5th prctile)';' -> G (95th prctile)';...
        ' -> U (median)'; 'p-value';'# estimates';'# lower estim.';'# upper estim.';' -> U (5th prctile)';' -> U (95th prctile)'};
    
    fid = fopen(filename,'w');
    for i = 1:length(Header)-1
        fprintf(fid,'%s,',char(Header(i)));
    end
    fprintf(fid,'%s\n',char(Header(end)));
    fclose(fid);
    %write Output values
    dlmwrite(filename,OutputMatrix,'-append')
    %%%%%-------------------
    
    filename = strcat('./PositionWiseKDEstimatesSummary.csv');
    try,
        delete(filename);
    catch
    end
    %generate Header
    Header = {'Pos.';'wt Base';'mut. with max effect';strcat(figtitle,' (median)');...
        'p-value';'# estimates';'# lower estim.';'# upper estim.';strcat(figtitle,' (5th prctile)');strcat(figtitle,' (95th prctile)')};
    fid = fopen(filename,'w');
    for i = 1:length(Header)-1
        fprintf(fid,'%s,',char(Header(i)));
    end
    fprintf(fid,'%s\n',char(Header(end)));
    fclose(fid);
    %write Output values
    dlmwrite(filename,OutputMatrix2,'-append')
end
TestPositionWiseDifferences130215;

%save('RelKds.mat','TotalRel_KD_m_w')