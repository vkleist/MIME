%function [] = CoSelectionAnalysis2(Position,MinimumNrReads,MinimumNrEpistasisValues)

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

%close all
%clear all
%with Watson-Crick and Wobble base- vs. non Watson-Crick base pairing
%           AA   AC   AG   AT   CA   CC   CG   CT   GA   GC   GG   GT   TA   TC   TG   TT
markers = {'ob';'ok';'oy';'pk';'or';'og';'pk';'ok';'ok';'pr';'oc';'sk';'pb';'om';'sb';'ok'};
FaceColors = {'b';'r';'k';'b';'r';  'g'; 'r'; 'm';  'k'; 'r';'c';'b';'b';'m';'b';'y'};
%markers = {'ob';'or';'ok';'og';'sb';'sr';'sk';'sg';'^b';'^r';'^k';'^g';'vb';'vr';'vk';'vg'};

%----------------------
ExperimentalPairs = [[1:3]',[1:3]'+9]; %wtExperiments
[NrComparisons,~] = size(ExperimentalPairs);
%Position = 266

MinimumNrReads = 5;
NrPositionsTotal = 535;
bases = 1:4;
CutValueFwd = 30;
CutValueBwd = 45;
PlotOutput = 0;
WeightThreshold = 0.6;%of minimum read coverage (for censoring)
%--------------------------

Exp_Random_VariableBeads = zeros(NrPositionsTotal,16,NrPositionsTotal); % Pos. x type of transition x Pos. 
Exp_RateOfFalseDetectionBeads = zeros(NrPositionsTotal,16,NrPositionsTotal); 
Exp_Random_VariableSuperN = zeros(NrPositionsTotal,16,NrPositionsTotal); %
Exp_RateOfFalseDetectionSuperN = zeros(NrPositionsTotal,16,NrPositionsTotal); %
%
Exp_Random_VariableTotal = zeros(NrPositionsTotal,16,NrPositionsTotal); %
Exp_RateOfFalseDetectionTotal = zeros(NrPositionsTotal,16,NrPositionsTotal); %


NativeBasePair = zeros(NrPositionsTotal,NrPositionsTotal);
Pos_m1_w = zeros(NrPositionsTotal,NrPositionsTotal,3); %Positions where only pos1 is mutated and pos 2 is in wt 
Pos_w_m1 = zeros(NrPositionsTotal,NrPositionsTotal,3); %Positions where only pos2 is mutated and pos 1 is in wt 
TotalNrSeqBeads = zeros(NrPositionsTotal,NrPositionsTotal);
TotalNrSeqSuperN = zeros(NrPositionsTotal,NrPositionsTotal);
%
TotalNrSeq = zeros(NrPositionsTotal,NrPositionsTotal);

%For re-ordering entries in the case where the first-and second base were
%swapt
TransitionRule =  zeros(1,20);
TransitionRule(1:4) = [3:4,1:2];
for i = 1:4
    for j = 1:4
        TransitionRule((i-1)*4+j+4) = (j-1)*4+i+4;
    end
end
%-----------------------------

%Read Reference sequence
disp('Read RefSeq..')
Tmp = dlmread('./Data_CoVariation/RefSeq.txt');
RefSeq = Tmp(:,2);

%----------------------------
%Binding Partners
PartnerNt = {};
PartnerNt{1} = 4;%A:U
PartnerNt{2} = 3;%C:G
PartnerNt{3} = [2 4];%G:C,G:U
PartnerNt{4} = [1 3];%U:A,U:G

for outercounter = CutValueFwd+1:NrPositionsTotal-CutValueBwd %250:250%<----------------------------------------------------------
    
    Position = outercounter

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

        [NrEntries,~] = size(M1);
        bases = 1:4;
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
        
        % --------------------------------Careful------------------------------ 
        M1(M1 == 0) = nan;
        M2(M2 == 0) = nan;
        
    
       
        for i = CutValueFwd+1:NrPositionsTotal-CutValueBwd%NrEntries   
            %wt positions
            pos1 = M1(i,2);
            pos2 = M1(i,4);
            
            idx1 = bases(1:4~=pos1);
            %wt_wt numbers
            Beads_w_w = M1(i,4*(pos1-1)+pos2+4);
            SuperN_w_w = M2(i,4*(pos1-1)+pos2+4);
            
            NativeBasePair(outercounter,i) = 4*(pos1-1)+pos2;
            Pos_m1_w(outercounter,i,:) = 4*(idx1-1)+pos2;
            
            NrSeqExperimentPosBeads = nansum(M1(i,5:end));
            NrSeqExperimentPosSuperN = nansum(M2(i,5:end));
            
            TotalNrSeqBeads(outercounter,i) = TotalNrSeqBeads(outercounter,i)+ NrSeqExperimentPosBeads;
            TotalNrSeqSuperN(outercounter,i) = TotalNrSeqSuperN(outercounter,i)+ NrSeqExperimentPosSuperN;
            TotalNrSeq(outercounter,i) = TotalNrSeq(outercounter,i) + NrSeqExperimentPosBeads + NrSeqExperimentPosSuperN;
            
            %Sum of probability of false detection sum(kappa) (weighted)
            Exp_RateOfFalseDetectionBeads(outercounter,:,i) = Exp_RateOfFalseDetectionBeads(outercounter,:,i) + M1(i,5:end)./Beads_w_w;
            Exp_RateOfFalseDetectionSuperN(outercounter,:,i) = Exp_RateOfFalseDetectionSuperN(outercounter,:,i) + M2(i,5:end)./SuperN_w_w;
            Exp_RateOfFalseDetectionTotal(outercounter,:,i) = Exp_RateOfFalseDetectionTotal(outercounter,:,i) + (M1(i,5:end)+ M2(i,5:end))/(Beads_w_w+SuperN_w_w);
            
            
            %Sum of probability of random numbers  sum(X)
            Exp_Random_VariableBeads(outercounter,:,i) = Exp_Random_VariableBeads(outercounter,:,i) + M1(i,5:end)./Phi(ExperimentalPair(1));
            Exp_Random_VariableSuperN(outercounter,:,i) =  Exp_Random_VariableSuperN(outercounter,:,i) + M2(i,5:end)./Theta(ExperimentalPair(1));
            Exp_Random_VariableTotal(outercounter,:,i) =  Exp_Random_VariableTotal(outercounter,:,i) +  M1(i,5:end)./Phi(ExperimentalPair(1))+M2(i,5:end)./Theta(ExperimentalPair(1));
        end
    end % end over experiments
end
%%
%average over experiments (including the weighing)
Exp_RateOfFalseDetectionBeads = Exp_RateOfFalseDetectionBeads./NrComparisons;%expected value over different experiments
Exp_RateOfFalseDetectionSuperN = Exp_RateOfFalseDetectionSuperN./NrComparisons;%expected value over different experiments
Exp_RateOfFalseDetectionTotal = Exp_RateOfFalseDetectionTotal./NrComparisons;

Exp_Random_VariableBeads = Exp_Random_VariableBeads./NrComparisons;%expected value over different experiments
Exp_Random_VariableSuperN = Exp_Random_VariableSuperN./NrComparisons;%expected value over different experiments
Exp_Random_VariableTotal = Exp_Random_VariableTotal./(NrComparisons*2);

PositionWeightsBeads = TotalNrSeqBeads./(max(TotalNrSeqBeads,[],2)*ones(1,NrPositionsTotal));
PositionWeightsSuperN = TotalNrSeqSuperN./(max(TotalNrSeqSuperN,[],2)*ones(1,NrPositionsTotal));
PositionWeightsTotal = TotalNrSeq./(max(TotalNrSeq,[],2)*ones(1,NrPositionsTotal));

% %%
% 
% %Evaluation of E(kappa_w->m1,w)....expected value over covarying positions
% 
% %Positions of potential binding partners (i.e. A:U)
% %This is because CoVar data only available for potential binding pairs
% Partners = {};
% 
% k = PartnerNt{1};
% PartnersX = [];
% for i = 1:length(k)
%     PartnersX = cat(1,PartnersX,find(RefSeq(CutValue+1:NrPositionsTotal-CutValue)==k(i))+CutValue);
% end
% Partners{1} = sort(PartnersX);
% 
% k = PartnerNt{2};
% PartnersX = [];
% for i = 1:length(k)
%     PartnersX = cat(1,PartnersX,find(RefSeq(CutValue+1:NrPositionsTotal-CutValue)==k(i))+CutValue);
% end
% Partners{2} = sort(PartnersX);
% 
% k = PartnerNt{3};
% PartnersX = [];
% for i = 1:length(k)
%     PartnersX = cat(1,PartnersX,find(RefSeq(CutValue+1:NrPositionsTotal-CutValue)==k(i))+CutValue);
% end
% Partners{3} = sort(PartnersX);
% 
% k = PartnerNt{4};
% PartnersX = [];
% for i = 1:length(k)
%     PartnersX = cat(1,PartnersX,find(RefSeq(CutValue+1:NrPositionsTotal-CutValue)==k(i))+CutValue);
% end
% Partners{4} = sort(PartnersX);
% 
% NrPartners= [length(Partners{1}),length(Partners{2}),length(Partners{3}),length(Partners{4})];
%%

%kappa: Beads
ExpKappa_w2m1_Beads_A2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'A'
ExpKappa_w2m1_Beads_C2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'C'
ExpKappa_w2m1_Beads_G2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'G'
ExpKappa_w2m1_Beads_T2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'T'

ExpKappa_w2m1_Beads_A2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'A'
ExpKappa_w2m1_Beads_C2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'C'
ExpKappa_w2m1_Beads_G2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'G'
ExpKappa_w2m1_Beads_T2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'T'

%random variable X: Beads
ExpX_w2m1_Beads_A2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'A'
ExpX_w2m1_Beads_C2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'C'
ExpX_w2m1_Beads_G2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'G'
ExpX_w2m1_Beads_T2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'T'

%kappa: SuperN
ExpKappa_w2m1_SuperN_A2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'A'
ExpKappa_w2m1_SuperN_C2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'C'
ExpKappa_w2m1_SuperN_G2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'G'
ExpKappa_w2m1_SuperN_T2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'T'

ExpKappa_w2m1_SuperN_A2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'A'
ExpKappa_w2m1_SuperN_C2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'C'
ExpKappa_w2m1_SuperN_G2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'G'
ExpKappa_w2m1_SuperN_T2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'T'

%random variable X: SuperN
ExpX_w2m1_SuperN_A2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'A'
ExpX_w2m1_SuperN_C2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'C'
ExpX_w2m1_SuperN_G2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'G'
ExpX_w2m1_SuperN_T2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'T'

%kappa: Pooled (beads + superN)
ExpKappa_w2m1_Total_A2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'A'
ExpKappa_w2m1_Total_C2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'C'
ExpKappa_w2m1_Total_G2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'G'
ExpKappa_w2m1_Total_T2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'T'

ExpKappa_w2m1_Total_A2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'A'
ExpKappa_w2m1_Total_C2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'C'
ExpKappa_w2m1_Total_G2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'G'
ExpKappa_w2m1_Total_T2Any = nan(NrPositionsTotal,NrPositionsTotal);%w->M1 where 'w' = 'T'

%random variable X:  Pooled (beads + superN)
ExpX_w2m1_Total_A2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'A'
ExpX_w2m1_Total_C2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'C'
ExpX_w2m1_Total_G2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'G'
ExpX_w2m1_Total_T2X = nan(NrPositionsTotal,NrPositionsTotal,3);%w->m1 where 'w' = 'T'


for outercounter = CutValueFwd+1:NrPositionsTotal-CutValueBwd %250:250%<----------------------------------------------------------
    wtbase = RefSeq(outercounter);
    %k = Partners{wtbase};
    for i =CutValueFwd+1:NrPositionsTotal-CutValueBwd%NrEntries  
        idx = reshape(Pos_m1_w(outercounter,i,:),1,3,1);
        
        Q1 = Exp_RateOfFalseDetectionBeads(outercounter,idx,i);
        Q2 = Exp_RateOfFalseDetectionSuperN(outercounter,idx,i);
        
        Q3 = Exp_Random_VariableBeads(outercounter,idx,i);
        Q4 = Exp_Random_VariableSuperN(outercounter,idx,i);
        
        Q5 = Exp_RateOfFalseDetectionTotal(outercounter,idx,i);
        Q6 = Exp_Random_VariableTotal(outercounter,idx,i);
        
        %<-----------censoring---------->
        if PositionWeightsBeads(outercounter,i)<WeightThreshold 
            Q1 = nan(1,3,1);
            Q3 = Q1;
        end
        if PositionWeightsSuperN(outercounter,i)<WeightThreshold 
            Q2 = nan(1,3,1);
            Q4 = Q1;   
        end
        
        if PositionWeightsTotal(outercounter,i)<WeightThreshold 
            Q5 = nan(1,3,1);
            Q6 = Q5;
        end
        %<---------------------------------->
        
        if wtbase == 1
            ExpKappa_w2m1_Beads_A2X(outercounter,i,:) = Q1;
            ExpKappa_w2m1_SuperN_A2X(outercounter,i,:) = Q2;
            ExpKappa_w2m1_Total_A2X(outercounter,i,:) = Q5;
            
            ExpKappa_w2m1_Beads_A2Any(outercounter,i) = nansum(Q1); 
            ExpKappa_w2m1_SuperN_A2Any(outercounter,i) = nansum(Q2); 
            ExpKappa_w2m1_Total_A2Any(outercounter,i) = nansum(Q5); 
            
            ExpX_w2m1_Beads_A2X(outercounter,i,:) = Q3;
            ExpX_w2m1_SuperN_A2X(outercounter,i,:) = Q4;
            ExpX_w2m1_Total_A2X(outercounter,i,:) = Q6;
        elseif wtbase == 2
            ExpKappa_w2m1_Beads_C2X(outercounter,i,:) = Q1;
            ExpKappa_w2m1_SuperN_C2X(outercounter,i,:) = Q2;
            ExpKappa_w2m1_Total_C2X(outercounter,i,:) = Q5;
            
            ExpKappa_w2m1_Beads_C2Any(outercounter,i) =  nansum(Q1); 
            ExpKappa_w2m1_SuperN_C2Any(outercounter,i) =  nansum(Q2); 
            ExpKappa_w2m1_Total_C2Any(outercounter,i) =  nansum(Q5); 
            
            ExpX_w2m1_Beads_C2X(outercounter,i,:) = Q3;
            ExpX_w2m1_SuperN_C2X(outercounter,i,:) = Q4;
            ExpX_w2m1_Total_C2X(outercounter,i,:) = Q6;
        elseif wtbase == 3
            ExpKappa_w2m1_Beads_G2X(outercounter,i,:) = Q1;
            ExpKappa_w2m1_SuperN_G2X(outercounter,i,:) = Q2;
            ExpKappa_w2m1_Total_G2X(outercounter,i,:) = Q5;
            
            ExpKappa_w2m1_Beads_G2Any(outercounter,i) = nansum(Q1); 
            ExpKappa_w2m1_SuperN_G2Any(outercounter,i) = nansum(Q2); 
            ExpKappa_w2m1_Total_G2Any(outercounter,i) = nansum(Q5); 
            
            ExpX_w2m1_Beads_G2X(outercounter,i,:) = Q3;
            ExpX_w2m1_SuperN_G2X(outercounter,i,:) = Q4;
            ExpX_w2m1_Total_G2X(outercounter,i,:) = Q6;
        elseif wtbase == 4
            ExpKappa_w2m1_Beads_T2X(outercounter,i,:) = Q1;
            ExpKappa_w2m1_SuperN_T2X(outercounter,i,:) = Q2;
            ExpKappa_w2m1_Total_T2X(outercounter,i,:) = Q5;
            
            ExpKappa_w2m1_Beads_T2Any(outercounter,i) =  nansum(Q1); 
            ExpKappa_w2m1_SuperN_T2Any(outercounter,i) =  nansum(Q2); 
            ExpKappa_w2m1_Total_T2Any(outercounter,i) =  nansum(Q5); 
            
            ExpX_w2m1_Beads_T2X(outercounter,i,:) = Q3;
            ExpX_w2m1_SuperN_T2X(outercounter,i,:) = Q4;
            ExpX_w2m1_Total_T2X(outercounter,i,:) = Q6;
        end
    end
    
end

ExpKappa_w2m1_Beads_A2Any(ExpKappa_w2m1_Beads_A2Any==0) = nan; 
ExpKappa_w2m1_SuperN_A2Any(ExpKappa_w2m1_SuperN_A2Any== 0) = nan; 
ExpKappa_w2m1_Total_A2Any(ExpKappa_w2m1_Total_A2Any==0) = nan; 

ExpKappa_w2m1_Beads_C2Any(ExpKappa_w2m1_Beads_C2Any==0) = nan; 
ExpKappa_w2m1_SuperN_C2Any(ExpKappa_w2m1_SuperN_C2Any== 0) = nan; 
ExpKappa_w2m1_Total_C2Any(ExpKappa_w2m1_Total_C2Any==0) = nan; 

ExpKappa_w2m1_Beads_G2Any(ExpKappa_w2m1_Beads_G2Any==0) = nan; 
ExpKappa_w2m1_SuperN_G2Any(ExpKappa_w2m1_SuperN_G2Any== 0) = nan; 
ExpKappa_w2m1_Total_G2Any(ExpKappa_w2m1_Total_G2Any==0) = nan; 

ExpKappa_w2m1_Beads_T2Any(ExpKappa_w2m1_Beads_T2Any==0) = nan; 
ExpKappa_w2m1_SuperN_T2Any(ExpKappa_w2m1_SuperN_T2Any== 0) = nan; 
ExpKappa_w2m1_Total_T2Any(ExpKappa_w2m1_Total_T2Any==0) = nan; 

%%
TitleLabel = ['A','C','G','T'];

% %check
pos = 273;
k = RefSeq(pos);
idx = bases(1:4~=k);
switch k
    case 1
        K1 = ExpKappa_w2m1_Beads_A2X(pos,:,:);
        K2 = ExpX_w2m1_Beads_A2X(pos,:,:);
        
        K3 = ExpKappa_w2m1_SuperN_A2X(pos,:,:);
        K4 = ExpX_w2m1_SuperN_A2X(pos,:,:);
        
        K5 = ExpKappa_w2m1_Total_A2X(pos,:,:);
        K6 = ExpX_w2m1_Total_A2X(pos,:,:);
    case 2
        K1 = ExpKappa_w2m1_Beads_C2X(pos,:,:);
        K2 = ExpX_w2m1_Beads_C2X(pos,:,:);
        
        K3 = ExpKappa_w2m1_SuperN_C2X(pos,:,:);
        K4 = ExpX_w2m1_SuperN_C2X(pos,:,:);
        
        K5 = ExpKappa_w2m1_Total_C2X(pos,:,:);
        K6 = ExpX_w2m1_Total_C2X(pos,:,:);
    case 3
        K1 = ExpKappa_w2m1_Beads_G2X(pos,:,:);
        K2 = ExpX_w2m1_Beads_G2X(pos,:,:);
        
        K3 = ExpKappa_w2m1_SuperN_G2X(pos,:,:);
        K4 = ExpX_w2m1_SuperN_G2X(pos,:,:);
        
        K5 = ExpKappa_w2m1_Total_G2X(pos,:,:);
        K6 = ExpX_w2m1_Total_G2X(pos,:,:);
    case 4
        K1 = ExpKappa_w2m1_Beads_T2X(pos,:,:);
        K2 = ExpX_w2m1_Beads_T2X(pos,:,:);
        
        K3 = ExpKappa_w2m1_SuperN_T2X(pos,:,:);
        K4 = ExpX_w2m1_SuperN_T2X(pos,:,:);
        
        K5 = ExpKappa_w2m1_Total_T2X(pos,:,:);
        K6 = ExpX_w2m1_Total_T2X(pos,:,:);
end
%%
XValues = CutValueFwd+1:NrPositionsTotal-CutValueBwd;
figure(20)
plot(XValues,K1(1,XValues,1),'o-')
hold on; 
plot(XValues,K1(1,XValues,2),'o-r')
plot(XValues,K1(1,XValues,3),'o-g')

plot(XValues,K3(1,XValues,1),'s-')
plot(XValues,K3(1,XValues,2),'s-r')
plot(XValues,K3(1,XValues,3),'s-g')

plot(XValues,K5(1,XValues,1),'d-')
plot(XValues,K5(1,XValues,2),'d-r')
plot(XValues,K5(1,XValues,3),'d-g')

legend(strcat(TitleLabel(k),'->',TitleLabel(idx(1)),'(b)'),strcat(TitleLabel(k),...
    '->',TitleLabel(idx(2)),'(b)'),strcat(TitleLabel(k),'->',TitleLabel(idx(3)),...
    '(b)'),strcat(TitleLabel(k),'->',TitleLabel(idx(1)),'(S)'),strcat(TitleLabel(k),...
    '->',TitleLabel(idx(2)),'(S)'),strcat(TitleLabel(k),'->',TitleLabel(idx(3)),'(S)'),...
    strcat(TitleLabel(k),'->',TitleLabel(idx(1)),'(P)'),strcat(TitleLabel(k),...
    '->',TitleLabel(idx(2)),'(P)'),strcat(TitleLabel(k),'->',TitleLabel(idx(3)),'(P)'))


title(strcat('Prob. of false detect. position ',num2str(pos), ' ' ,TitleLabel(k),'->X'))
ylabel('E(\kappa_{w \rightarrow m1,w})')
xlabel('Position')
set(gca,'YScale','log')

figure(21)
plot(XValues,K2(1,XValues,1),'o-')
hold on; 
plot(XValues,K2(1,XValues,2),'o-r')
plot(XValues,K2(1,XValues,3),'o-g')

plot(XValues,K4(1,XValues,1),'s-')
plot(XValues,K4(1,XValues,2),'s-r')
plot(XValues,K4(1,XValues,3),'s-g')

plot(XValues,K6(1,XValues,1),'d-')
plot(XValues,K6(1,XValues,2),'d-r')
plot(XValues,K6(1,XValues,3),'d-g')

legend(strcat(TitleLabel(k),'->',TitleLabel(idx(1)),'(b)'),strcat(TitleLabel(k),'->',...
    TitleLabel(idx(2)),'(b)'),strcat(TitleLabel(j),'->',TitleLabel(idx(3)),'(b)'),...
    strcat(TitleLabel(k),'->',TitleLabel(idx(1)),'(S)'),strcat(TitleLabel(k),'->',...
    TitleLabel(idx(2)),'(S)'),strcat(TitleLabel(j),'->',TitleLabel(idx(3)),'(S)'),...
    strcat(TitleLabel(k),'->',TitleLabel(idx(1)),'(P)'),strcat(TitleLabel(k),'->',...
    TitleLabel(idx(2)),'(P)'),strcat(TitleLabel(j),'->',TitleLabel(idx(3)),'(P)'))

title(strcat('Random Variable. position ',num2str(pos), ' ' ,TitleLabel(k),'->X'))
ylabel('E(X_{w \rightarrow m1,w})')
xlabel('Position')
set(gca,'YScale','log')

%%


BoxPlotColors = ['b','g','r','m'];
for j = 1:4 %wt Base
    SubPlotIdx = bases(1:4~=j);

    Positions = find(RefSeq == j);
    counter = 1;
    for q = 1:3%SubPlotIdx %Go over MutantBases (3 in total)
        c = SubPlotIdx(q);
        switch j
            case 1
                Mut1 = ExpKappa_w2m1_Beads_A2X(Positions,:,q);
                Mut2 = ExpKappa_w2m1_SuperN_A2X(Positions,:,q);
                Mut3 = ExpKappa_w2m1_Total_A2X(Positions,:,q);
            case 2
                Mut1 = ExpKappa_w2m1_Beads_C2X(Positions,:,q);
                Mut2 = ExpKappa_w2m1_SuperN_C2X(Positions,:,q);
                Mut3 = ExpKappa_w2m1_Total_C2X(Positions,:,q);
            case 3
                Mut1 = ExpKappa_w2m1_Beads_G2X(Positions,:,q);
                Mut2 = ExpKappa_w2m1_SuperN_G2X(Positions,:,q);
                Mut3 = ExpKappa_w2m1_Total_G2X(Positions,:,q);
            case 4
                Mut1 = ExpKappa_w2m1_Beads_T2X(Positions,:,q);
                Mut2 = ExpKappa_w2m1_SuperN_T2X(Positions,:,q);
                Mut3 = ExpKappa_w2m1_Total_T2X(Positions,:,q);
        end
        figure(j)
        hold on
%         boxplot(log10(Mut1)','plotstyle','compact','colors',BoxPlotColors(c),'positions',[1:length(Positions)]+(counter-2)/5);
%         set(gca,'XTickLabel',{' '})
%         boxplot(log10(Mut2)','plotstyle','compact','colors',BoxPlotColors(c),'positions',[1:length(Positions)]+(counter-2)/5);
%         set(gca,'XTickLabel',{' '})
        boxplot(log10(Mut3)','plotstyle','compact','colors',BoxPlotColors(c),'positions',[1:length(Positions)]+(counter-2)/5);
        set(gca,'XTickLabel',{' '})
        plot([1:length(Positions)]+(counter-2)/5,log10(nanmean(Mut3,2)),'Color',BoxPlotColors(c))
        set(gca,'XTickLabel',{' '})
        
        %Mut2 = reshape(RandomVariable(Positions,q,:),length(Positions),k,1);
        figure(j+4)
        hold on
%         plot([1:length(Positions)]+(counter-2)/5,nanstd(Mut1,1,2)./nanmean(Mut1,2).*100,'Color',BoxPlotColors(c))
%         set(gca,'XTickLabel',{' '})
%         plot([1:length(Positions)]+(counter-2)/5,nanstd(Mut2,1,2)./nanmean(Mut2,2).*100,'Color',BoxPlotColors(c),'LineStyle','--')
%         set(gca,'XTickLabel',{' '})
        plot([1:length(Positions)]+(counter-2)/5,nanstd(Mut3,1,2)./nanmean(Mut3,2).*100,'Color',BoxPlotColors(c),'LineStyle','-')
        set(gca,'XTickLabel',{' '})
        
        
        figure(j+8)
        hold on
        plot([1:length(Positions)]+(counter-2)/5,(nanmean(Mut1,2)-nanmean(Mut3,2))./nanmean(Mut3,2).*100,'Color',BoxPlotColors(c))
        set(gca,'XTickLabel',{' '})
        plot([1:length(Positions)]+(counter-2)/5,(nanmean(Mut2,2)-nanmean(Mut3,2))./nanmean(Mut3,2).*100,'Color',BoxPlotColors(c),'LineStyle','--')
        set(gca,'XTickLabel',{' '})
        
        counter = counter +1;
    end
    figure(j)
    axis([1 length(Positions) -4.7 -2])
    set(gca,'XTick',1:10:length(Positions),'XTickLabel',num2str(Positions(1:10:end)))
    %set(gca,'XTick',1:length(Positions),'XTickLabel',num2str(Positions(1:end)))
    xlabel('position')
    title(strcat('Wt Base: ', TitleLabel(j), '; False detection rate ->A (blue), ->C (green), ->G (red) and ->T (magenta)'))
    ylabel('log_{10} \kappa_{w \rightarrow m}')
    set(get(gca,'xlabel'),'Fontsize',12);
    set(get(gca,'ylabel'),'Fontsize',12);
    set(get(gca,'title'),'Fontsize',12);
    
    figure(j+4)
    axis([1 length(Positions) 0 200])
    set(gca,'XTick',1:10:length(Positions),'XTickLabel',num2str(Positions(1:10:end)))
    %set(gca,'XTick',1:length(Positions),'XTickLabel',num2str(Positions(1:end)))
    xlabel('position')
    title(strcat('Wt Base: ', TitleLabel(j), '; False detection rate ->A (blue), ->C (green), ->G (red) and ->T (magenta)'))
    ylabel('Coefficient of Variation \kappa_{w \rightarrow m} (% of mean)')
    set(get(gca,'xlabel'),'Fontsize',12);
    set(get(gca,'ylabel'),'Fontsize',12);
    set(get(gca,'title'),'Fontsize',12);
    
    figure(j+8)
    axis([1 length(Positions) -200 200])
    set(gca,'XTick',1:10:length(Positions),'XTickLabel',num2str(Positions(1:10:end)))
    %set(gca,'XTick',1:length(Positions),'XTickLabel',num2str(Positions(1:end)))
    xlabel('position')
    title(strcat('Wt Base: ', TitleLabel(j), '; False detection rate ->A (blue), ->C (green), ->G (red) and ->T (magenta)'))
    ylabel('difference in means between beads and SuperN \kappa_{w \rightarrow m} (%) with pooled data')
    set(get(gca,'xlabel'),'Fontsize',12);
    set(get(gca,'ylabel'),'Fontsize',12);
    set(get(gca,'title'),'Fontsize',12);
end
%%
MedianExpKappa_w2m1_Total_A2X = (reshape(nanmedian(ExpKappa_w2m1_Total_A2X,2),NrPositionsTotal,3,1))';%type of transition x position
MedianExpKappa_w2m1_Total_C2X = (reshape(nanmedian(ExpKappa_w2m1_Total_C2X,2),NrPositionsTotal,3,1))';%type of transition x position
MedianExpKappa_w2m1_Total_G2X = (reshape(nanmedian(ExpKappa_w2m1_Total_G2X,2),NrPositionsTotal,3,1))';%type of transition x position
MedianExpKappa_w2m1_Total_T2X = (reshape(nanmedian(ExpKappa_w2m1_Total_T2X,2),NrPositionsTotal,3,1))';%type of transition x position

MedianExpKappa_w2m1_Total = zeros(16,NrPositionsTotal);
MedianExpKappa_w2m1_Total(2:4,:) = MedianExpKappa_w2m1_Total_A2X;
MedianExpKappa_w2m1_Total([5,7:8],:) = MedianExpKappa_w2m1_Total_C2X;
MedianExpKappa_w2m1_Total([9:10,12],:) = MedianExpKappa_w2m1_Total_G2X;
MedianExpKappa_w2m1_Total([13:15],:) = MedianExpKappa_w2m1_Total_T2X;

%
MedianExpKappa_w2m1_Total_A2Any = nanmedian(ExpKappa_w2m1_Total_A2Any,2);%all types of transition x position
MedianExpKappa_w2m1_Total_C2Any = nanmedian(ExpKappa_w2m1_Total_C2Any,2);%all types of transition x position
MedianExpKappa_w2m1_Total_G2Any = nanmedian(ExpKappa_w2m1_Total_G2Any,2);%all types of transition x position
MedianExpKappa_w2m1_Total_T2Any = nanmedian(ExpKappa_w2m1_Total_T2Any,2);%all types of transition x position

MedianExpKappa_w2M1_Total= zeros(NrPositionsTotal,1);
MedianExpKappa_w2M1_Total(RefSeq == 1) = MedianExpKappa_w2m1_Total_A2Any(RefSeq == 1);
MedianExpKappa_w2M1_Total(RefSeq == 2) = MedianExpKappa_w2m1_Total_C2Any(RefSeq == 2);
MedianExpKappa_w2M1_Total(RefSeq == 3) = MedianExpKappa_w2m1_Total_G2Any(RefSeq == 3);
MedianExpKappa_w2M1_Total(RefSeq == 4) = MedianExpKappa_w2m1_Total_T2Any(RefSeq == 4);
%MedianExpKappa_w2M1_Total = [MedianExpKappa_w2m1_Total_A2Any,MedianExpKappa_w2m1_Total_C2Any,MedianExpKappa_w2m1_Total_G2Any,MedianExpKappa_w2m1_Total_T2Any];
%

MeanExpKappa_w2m1_Total_A2X = (reshape(nanmean(ExpKappa_w2m1_Total_A2X,2),NrPositionsTotal,3,1))';%type of transition x position
MeanExpKappa_w2m1_Total_C2X = (reshape(nanmean(ExpKappa_w2m1_Total_C2X,2),NrPositionsTotal,3,1))';%type of transition x position
MeanExpKappa_w2m1_Total_G2X = (reshape(nanmean(ExpKappa_w2m1_Total_G2X,2),NrPositionsTotal,3,1))';%type of transition x position
MeanExpKappa_w2m1_Total_T2X = (reshape(nanmean(ExpKappa_w2m1_Total_T2X,2),NrPositionsTotal,3,1))';%type of transition x position

MeanExpKappa_w2m1_Total = zeros(16,NrPositionsTotal);
MeanExpKappa_w2m1_Total(2:4,:) = MeanExpKappa_w2m1_Total_A2X;
MeanExpKappa_w2m1_Total([5,7:8],:) = MeanExpKappa_w2m1_Total_C2X;
MeanExpKappa_w2m1_Total([9:10,12],:) = MeanExpKappa_w2m1_Total_G2X;
MeanExpKappa_w2m1_Total([13:15],:) = MeanExpKappa_w2m1_Total_T2X;

%
MeanExpKappa_w2m1_Total_A2Any = nanmean(ExpKappa_w2m1_Total_A2Any,2);%all types of transition x position
MeanExpKappa_w2m1_Total_C2Any = nanmean(ExpKappa_w2m1_Total_C2Any,2);%all types of transition x position
MeanExpKappa_w2m1_Total_G2Any = nanmean(ExpKappa_w2m1_Total_G2Any,2);%all types of transition x position
MeanExpKappa_w2m1_Total_T2Any = nanmean(ExpKappa_w2m1_Total_T2Any,2);%all types of transition x position

MeanExpKappa_w2M1_Total= zeros(NrPositionsTotal,1);
MeanExpKappa_w2M1_Total(RefSeq == 1) = MeanExpKappa_w2m1_Total_A2Any(RefSeq == 1);
MeanExpKappa_w2M1_Total(RefSeq == 2) = MeanExpKappa_w2m1_Total_C2Any(RefSeq == 2);
MeanExpKappa_w2M1_Total(RefSeq == 3) = MeanExpKappa_w2m1_Total_G2Any(RefSeq == 3);
MeanExpKappa_w2M1_Total(RefSeq == 4) = MeanExpKappa_w2m1_Total_T2Any(RefSeq == 4);
%MeanExpKappa_w2M1_Total = [MeanExpKappa_w2m1_Total_A2Any,MeanExpKappa_w2m1_Total_C2Any,MeanExpKappa_w2m1_Total_G2Any,MeanExpKappa_w2m1_Total_T2Any];
%

filename = './ErrorEstimates.mat';

try delete(filename)
    pause(1);
catch
end
save('ErrorEstimates.mat', 'MedianExpKappa_w2m1_Total','MeanExpKappa_w2m1_Total','MedianExpKappa_w2M1_Total','MeanExpKappa_w2M1_Total');
