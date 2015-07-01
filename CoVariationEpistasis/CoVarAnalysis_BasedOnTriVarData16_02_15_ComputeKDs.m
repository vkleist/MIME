function [] = CoVarAnalysis_BasedOnTriVarData16_02_15_ComputeKDs(FirstPosition,SecondPositions,FirstPos3,LastPos3,ExperimentalPairs,Outfolder)


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

%----------------------
NrPositionsTotal3 = LastPos3-FirstPos3;
%ExperimentalPairs = [[4]',[4]'+9]; 
%ExperimentalPairs = [[4:9]',[4:9]'+9]; 
bases = 1:4;
%-----------------------------

%Read Reference sequence
% disp('Read RefSeq..')
% Tmp = dlmread('./Data_CoVariation/RefSeq.txt');
% RefSeq = Tmp(:,2);
% %Binding Partners
% PartnerNt = {};
% PartnerNt{1} = 4;%A:U
% PartnerNt{2} = 3;%C:G
% PartnerNt{3} = [2 4];%G:C,G:U
% PartnerNt{4} = [1 3];%U:A,U:G
%----------------------------
[NrComparisons,~] = size(ExperimentalPairs);

%load Error estimates: gives error matrices type of transition x position
load './ErrorEstimates.mat';

%individual mutation probabilities, e.g. A -> C
MedianExpKappa_w2m1_Total(MedianExpKappa_w2m1_Total==0) = nan;
%MeanExpKappa_w2m1_Total(MeanExpKappa_w2m1_Total==0) = nan;

InFolder = './DataTriVariation_200/';
%tic
% go through positions r_1

Position1 = FirstPosition;
%disp(num2str(Position1));
for outercounter2 = SecondPositions %250:250%<-----
    fprintf('%s', strcat('pos2:',num2str(outercounter2)));
    Position2 = outercounter2;
    TotalM1 = zeros(NrPositionsTotal3-1,70,NrComparisons); % take together all experiments [Gag Concentrations]
    TotalM2 = TotalM1; % take together all experiments [Gag Concentrations]

    for counter = 1:NrComparisons %start over experiments
        ExperimentalPair = ExperimentalPairs(counter,:);

        InFolderExp = strcat(InFolder,num2str(ExperimentalPair(1)),'/');
        filename = strcat(num2str(Position1),'_',num2str(Position2));
        M1 = dlmread(strcat(InFolderExp,filename,'.txt')); %Data @ Beads

        InFolderExp = strcat(InFolder,num2str(ExperimentalPair(2)),'/');
        %filename2 = strcat(num2str(Position1),'_',num2str(Position2),'.txt');
        M2 = dlmread(strcat(InFolderExp,filename,'.txt')); %Data @ S/N

        %Matrix M: 
        % pos1 Nt@pos1 pos2 Nt@pos2 pos3 Nr@pos3 AAA AAC AAG AAT
        % ACA...TTT

        %[NrEntries,~] = size(M1);

        %---------

        TotalM1(:,1:6,counter) = M1(:,1:6);
        TotalM1(:,7:end,counter) = M1(:,7:end);
        TotalM2(:,1:6,counter) = M2(:,1:6);
        TotalM2(:,7:end,counter) = M2(:,7:end);

    end % end over experiments

    TotalRel_KD_m_w_w = nan((NrPositionsTotal3-1)*NrComparisons,4);
    TotalRel_KD_w_m_w = nan((NrPositionsTotal3-1)*NrComparisons,4);
    TotalRel_KD_m_m_w = nan((NrPositionsTotal3-1)*NrComparisons,16);

    %Epistasis = nan(NrPositionsTotal3-1,16);

    Signal2Noise_m_w_w_SuperN = zeros((NrPositionsTotal3-1)*NrComparisons,4);
    Signal2Noise_w_m_w_SuperN = zeros((NrPositionsTotal3-1)*NrComparisons,4);
    Signal2Noise_m_m_w_SuperN = zeros((NrPositionsTotal3-1)*NrComparisons,16);

    Signal2Noise_m_w_w_Beads = Signal2Noise_m_w_w_SuperN ;
    Signal2Noise_w_m_w_Beads = Signal2Noise_w_m_w_SuperN ;
    Signal2Noise_m_m_w_Beads = Signal2Noise_m_m_w_SuperN ;

    NrOccurances_m_m_w_Beads = Signal2Noise_m_m_w_SuperN ;
    
    PositionReadsBeads = zeros(NrPositionsTotal3-1,NrComparisons); 
    PositionReadsSuperN = zeros(NrPositionsTotal3-1,NrComparisons); 
    PositionReadsTotal = zeros(NrPositionsTotal3-1,NrComparisons); 

    % go through positions r_3
    [n,~,~] = size(TotalM1);

    %wt positions
    pos1 = TotalM1(1,2,1);
    pos2 = TotalM1(1,4,1);
    %-----------------------------------------------------------------
    RunningIndex = 0;
    %----------------------------------------------------------------------
    % set these ones to nan
    Positions3 = TotalM1(:,5,1);
    Indices = ~(Positions3 < Position1 | Positions3 > Position2);
    %----------------------------------------------------------------------

    for i = 1:n%   
        %wt positions
        pos3 = TotalM1(i,6,1);

        idx1 = 1:4;% bases(1:4~=pos1);
        idx2 = 1:4;%bases(1:4~=pos2);
        %idx3 = 1:4;%bases(1:4~=pos2);

        pos_w_w_w = (pos1-1)*16 + (pos2-1)*4 + pos3 + 6;
        pos_m_w_w = (idx1-1).*16 + (pos2-1)*4 + pos3 + 6;
        pos_w_m_w = (pos1-1)*16 + (idx2-1).*4 + pos3 + 6;

        for counter = 1:NrComparisons 

            idx1 = 1:4;% bases(1:4~=pos1);
            idx2 = 1:4;%bases(1:4~=pos2);
            idx3 = 1:4;%bases(1:4~=pos2);

            RunningIndex = RunningIndex+1;
            %------------------Censoring-----------------------------------
            NrSeqExperimentPosBeads = nansum(TotalM1(i,7:end,counter));
            NrSeqExperimentPosSuperN = nansum(TotalM2(i,7:end,counter));

            PositionReadsBeads(i,counter) = NrSeqExperimentPosBeads;
            PositionReadsSuperN(i,counter) = NrSeqExperimentPosSuperN;
            PositionReadsTotal(i,counter) = NrSeqExperimentPosBeads + NrSeqExperimentPosSuperN;
            %------------------Censoring-----------------------------------

            %-----------------  wt_wt numbers  -----------------------
            Beads_w_w_w = TotalM1(i,pos_w_w_w,counter);
            SuperN_w_w_w = TotalM2(i,pos_w_w_w,counter);

            %-----------------  m_wt numbers  -----------------------
            Beads_m_w_w = TotalM1(i,pos_m_w_w,counter);
            SuperN_m_w_w = TotalM2(i,pos_m_w_w,counter);

            %-----------------  wt_m numbers  -----------------------
            Beads_w_m_w = TotalM1(i,pos_w_m_w,counter);
            SuperN_w_m_w = TotalM2(i,pos_w_m_w,counter);

            %----------------------------------------------------------
            %   KD(m_wt_wt)
            %----------------------------------------------------------

            Nominator = SuperN_m_w_w/SuperN_w_w_w - MedianExpKappa_w2m1_Total(4*(pos1-1)+idx1,Position1)';%superN
            Denominator = Beads_m_w_w/Beads_w_w_w - MedianExpKappa_w2m1_Total(4*(pos1-1)+idx1,Position1)';%beads

            Rel_KD_m_w_w = zeros(1,length(idx1));
            %Signal-2-Noise
            for j = idx1
                Noise = MedianExpKappa_w2m1_Total(4*(pos1-1)+j,Position1);

                Signal2Noise_m_w_w_SuperN(RunningIndex,j) = SuperN_m_w_w(j)/(SuperN_w_w_w*Noise);
                Signal2Noise_m_w_w_Beads(RunningIndex,j) = Beads_m_w_w(j)/(Beads_w_w_w*Noise);

                if SuperN_m_w_w(j)/(SuperN_w_w_w*Noise) == inf%% not enough signal at supernatent => signal-to-noise-ratio too bad
                    Rel_KD_m_w_w(j) = nan; 
                elseif Denominator(j) <= 0 %Denominator <= 0 && Nominator > 0 %not enough signal at beads
                    Rel_KD_m_w_w(j)  = -1;%Nominator/1e-8; <-------------------------------------------arbitrary
                    %LowerLimitsKD_m_w(Position1,RunningIndex,j) = 1; 
                elseif Nominator(j) <= 0
                    Rel_KD_m_w_w(j)  = -10;%Nominator/1e-8; <-------------------------------------------arbitrary
                else
                    Rel_KD_m_w_w(j) = Nominator(j)/Denominator(j); 
                end
            end

            TotalRel_KD_m_w_w(RunningIndex,:) = Rel_KD_m_w_w;

            %----------------------------------------------------------
            %   KD(wt_m_wt)
            %----------------------------------------------------------

            Nominator = SuperN_w_m_w/SuperN_w_w_w - MedianExpKappa_w2m1_Total(4*(pos2-1)+idx2,Position2)';%superN
            Denominator = Beads_w_m_w/Beads_w_w_w - MedianExpKappa_w2m1_Total(4*(pos2-1)+idx2,Position2)';%beads

            Rel_KD_w_m_w = zeros(1,length(idx2));
            %Signal-2-Noise
            for j = idx2
                Noise = MedianExpKappa_w2m1_Total(4*(pos2-1)+j,Position2);

                Signal2Noise_w_m_w_SuperN(RunningIndex,j) = SuperN_w_m_w(j)/(SuperN_w_w_w*Noise);
                Signal2Noise_w_m_w_Beads(RunningIndex,j) = Beads_w_m_w(j)/(Beads_w_w_w*Noise);

                if SuperN_w_m_w(j)/(SuperN_w_w_w*Noise) == inf%% not enough signal at supernatent => signal-to-noise-ratio too bad
                    Rel_KD_w_m_w(j) = nan; 
                elseif Denominator(j) <= 0
                    Rel_KD_w_m_w(j)  = -1;%Nominator/1e-8; <-------------------------------------------arbitrary
                    %LowerLimitsKD_w_m(Position1,RunningIndex,j) = 1; 
                elseif Nominator(j) <= 0
                    Rel_KD_w_m_w(j)  = -10;%Nominator/1e-8; <-------------------------------------------arbitrary
                else
                    Rel_KD_w_m_w(j) = Nominator(j)/Denominator(j); 
                end
            end
            TotalRel_KD_w_m_w(RunningIndex,:) = Rel_KD_w_m_w;

            %----------------------------------------------------------
            %   KD(m_m_wt)
            %----------------------------------------------------------

            %         if any([WCBasePairs,WobleBasePairs] == NativeBasePair(outercounter,i))
            idx1 = bases(1:4~=pos1);
            idx2 = bases(1:4~=pos2);
            for k = idx1% nt@pos1
                for l = idx2% nt@pos2
                    %----------------  m_m numbers  -------------------
                    pos_m_m_w = (k-1).*16 + (l-1).*4 + pos3 + 6;
                    Beads_m_m_w = TotalM1(i,pos_m_m_w,counter);
                    SuperN_m_m_w = TotalM2(i,pos_m_m_w,counter);

                    %------------------ Noise -------------------------
                    %a) double mutation ~ k_wt_m*k_m_wt
                    Noise_w1w2_m1m2 = MedianExpKappa_w2m1_Total(4*(pos1-1)+k,Position1)*MedianExpKappa_w2m1_Total(4*(pos2-1)+l,Position2);

                    %b) additional single mutation2 starting from single
                    %mutant 1: R_m_w_w/R_w_w_w * k_w_m
                    Noise_w1m2_m1m2_Beads = Beads_m_w_w(k)./Beads_w_w_w*MedianExpKappa_w2m1_Total(4*(pos2-1)+l,Position2);
                    Noise_w1m2_m1m2_SuperN = SuperN_m_w_w(k)./SuperN_w_w_w*MedianExpKappa_w2m1_Total(4*(pos2-1)+l,Position2);
                    %c) additional single mutation1 starting from single
                    %mutant 2: R_w_m_w/R_w_w_w * k_m_w
                    Noise_m1w2_m1m2_Beads = Beads_w_m_w(l)./Beads_w_w_w*MedianExpKappa_w2m1_Total(4*(pos1-1)+k,Position1);
                    Noise_m1w2_m1m2_SuperN = SuperN_w_m_w(l)./SuperN_w_w_w*MedianExpKappa_w2m1_Total(4*(pos1-1)+k,Position1);

                    Noise_Beads = Noise_w1w2_m1m2 + Noise_w1m2_m1m2_Beads + Noise_m1w2_m1m2_Beads;
                    Noise_SuperN = Noise_w1w2_m1m2 + Noise_w1m2_m1m2_SuperN + Noise_m1w2_m1m2_SuperN;
                    
                    Nominator = SuperN_m_m_w/SuperN_w_w_w - Noise_SuperN;%superN
                    Denominator = Beads_m_m_w/Beads_w_w_w - Noise_Beads;%beads
                    
                    NrOccurances_m_m_w_Beads(RunningIndex,4*(k-1)+l) = Beads_m_m_w;
                    
                    Signal2Noise_m_m_w_SuperN(RunningIndex,4*(k-1)+l) = SuperN_m_m_w/(SuperN_w_w_w*Noise_SuperN);
                    Signal2Noise_m_m_w_Beads(RunningIndex,4*(k-1)+l) = Beads_m_m_w/(Beads_w_w_w*Noise_Beads);

                    if SuperN_m_m_w/(SuperN_w_w_w*Noise_SuperN) == inf%% not enough signal at supernatent => signal-to-noise-ratio too bad
                        Rel_KD_m_m_w = nan; 
                        %Rel_KA_m_m_w = nan; 
                    elseif Denominator <= 0
                        Rel_KD_m_m_w  = -1;%Nominator/1e-8; <-------------------------------------------arbitrary
                        %LowerLimitsKD_w_m(Position1,RunningIndex,j) = 1; 
                        %Rel_KA_m_m_w = nan; 
                    elseif Nominator <= 0
                        Rel_KD_m_m_w  = -10;%Nominator/1e-8; <-------------------------------------------arbitrary
                        %Rel_KA_m_m_w = nan;
                    else
                        Rel_KD_m_m_w = Nominator/Denominator; 
                        %Rel_KA_m_m_w = Denominator/Nominator; 
                    end

                    TotalRel_KD_m_m_w(RunningIndex,4*(k-1)+l) = Rel_KD_m_m_w;
                    %Epistasis(RunningIndex,pos_m_m_w) = log(1/Rel_KD_m_m_w) - log(1./Rel_KD_m_w_w(k)*1/Rel_KD_w_m_w(l));
%                      end
                end
            end

        end%counter (experiments)
    end% i (pos 3)

    % censoring by Number of Reads
    PositionWeightsBeads = PositionReadsBeads./(ones(n,1)*max(PositionReadsBeads));
    PositionWeightsSuperN = PositionReadsSuperN./(ones(n,1)*max(PositionReadsSuperN));
    PositionWeightsTotal = PositionReadsTotal./(ones(n,1)*max(PositionReadsTotal));
    
    
    PositionWeightsBeads(Indices) = 0;
    PositionWeightsSuperN(Indices) = 0;
    PositionWeightsWeightsTotal(Indices) = 0;
    
% toc
%          %write *.mat file i
    save(strcat(Outfolder,filename,'.mat'),'FirstPos3','LastPos3','NrPositionsTotal3','NrComparisons','TotalRel_KD_m_w_w','TotalRel_KD_w_m_w',...
        'TotalRel_KD_m_m_w','Signal2Noise_m_w_w_SuperN','Signal2Noise_m_w_w_Beads','Signal2Noise_w_m_w_SuperN','Signal2Noise_w_m_w_Beads',...
    'Signal2Noise_m_m_w_SuperN','Signal2Noise_m_m_w_Beads','PositionWeightsBeads','PositionWeightsSuperN','PositionWeightsTotal',...
    'PositionReadsBeads','PositionReadsSuperN','PositionReadsTotal','NrOccurances_m_m_w_Beads');
    clear TotalRel_KD_m_w_w TotalRel_KD_w_m_w TotalRel_KD_m_m_w Signal2Noise_m_w_w_SuperN Signal2Noise_m_w_w_Beads
    clear Signal2Noise_w_m_w_SuperN Signal2Noise_w_m_w_Beads Signal2Noise_m_m_w_SuperN Signal2Noise_m_m_w_Beads PositionWeightsBeads
    clear PositionWeightsSuperN PositionWeightsTotal PositionReadsBeads PositionReadsSuperN PositionReadsTotal NrOccurances_m_m_w_Beads

end% outercounter2
