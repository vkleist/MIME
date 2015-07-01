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

%----------------------
close all
clear all
alpha = 0.05;%significance criteria for KD values
FirstPos = 200;
LastPos = 349;
%NrPositionsTotal = LastPos-FirstPos;
ExperimentalPairs = [[4:9]',[4:9]'+9];
Outfolder = './ResultsTriVarWCWB_new/';
bases = 1:4;
OnlyWC = 2; %0 = all base pairs, 1 = only potential WC Base pairs, 2 = WC or WB base pairs
MinimumEffect = 1;
%--------------------------

try,
    cd(Outfolder)
    cd('../')
catch
    mkdir(Outfolder);
end

%Read Reference sequence
disp('Read RefSeq..')
Tmp = dlmread('./Data_CoVariation/RefSeq.txt');
RefSeq = Tmp(:,2);

disp('Reading xls....');
M = dlmread('./PositionWiseKDEstimates.csv',',',1,0);
disp('done');

%--- 1. Extract important positions (currently based on: wt -> any mut)
%     Pvalues = M(:,32);
%     MedianEffect = M(:,31);
%     ImportantResiduesLogic = find(all([(Pvalues < alpha)';(MedianEffect>MinimumEffect)']));
Pvalues = M(:,[5 12 19 26]);
MedianEffect = M(:,[4 11 18 25]);
ImportantResiduesLogic = find(sum((Pvalues < alpha)+(MedianEffect>MinimumEffect) >= 2,2));

PositionNumbers = M(:,1);

ImportantResidues = PositionNumbers(ImportantResiduesLogic);
ImportantResidues = ImportantResidues(ImportantResidues>=FirstPos&ImportantResidues<=LastPos)';

dlmwrite(strcat(Outfolder,'ImportantResidues.csv'),ImportantResidues');
BaseOfImportantResidues = RefSeq(ImportantResidues);

if OnlyWC == 1
    %WC binding partners
    PartnerNt = zeros(4,1);
    PartnerNt(1) = 4;%A:U
    PartnerNt(2) = 3;%C:G
    PartnerNt(3) = 2;%G:C,G:U
    PartnerNt(4) = 1;%U:A,U:G
elseif OnlyWC == 2
    PartnerNt = zeros(4,2);
    PartnerNt(1,:) = [4 0];%A:U
    PartnerNt(2,:) = [3 0];%C:G
    PartnerNt(3,:) = [2 4];%G:C,G:U
    PartnerNt(4,:) = [1 3];%U:A,U:G
end

%--- 2. Extract potential interaction partners
%--- non-redundant upper triangular matrix.
NrPairs = 0;
for i = 1:length(ImportantResidues)-1
    FirstPosition = ImportantResidues(i);
    disp(strcat('Pos1: ',num2str(FirstPosition)))
    BaseofResidue1 = BaseOfImportantResidues(i);
    if OnlyWC == 1
        PotentialBindingPartnersLogic = all([(BaseOfImportantResidues == PartnerNt(BaseofResidue1))';ImportantResidues > ImportantResidues(i)]);
    elseif OnlyWC == 0
        PotentialBindingPartnersLogic = ImportantResidues > ImportantResidues(i);
    elseif OnlyWC == 2
        PotentialBindingPartnersLogic = all([any([(BaseOfImportantResidues == PartnerNt(BaseofResidue1,1))';(BaseOfImportantResidues == PartnerNt(BaseofResidue1,2))']);ImportantResidues > ImportantResidues(i)]);
    end
    PotentialBindingPartners = ImportantResidues(PotentialBindingPartnersLogic)
    dlmwrite(strcat(Outfolder,num2str(FirstPosition),'_Partners.csv'),PotentialBindingPartners');
    if ~isempty(PotentialBindingPartners)
        SecondPositions = PotentialBindingPartners;
        %for j = PotentialBindingPartners
        CoVarAnalysis_BasedOnTriVarData16_02_15_ComputeKDs(FirstPosition,SecondPositions,FirstPos,LastPos,ExperimentalPairs,Outfolder)
        NrPairs = NrPairs+length(SecondPositions);
        %end
    end
end
dlmwrite(strcat(Outfolder,'NumberPairs.csv'),NrPairs);


