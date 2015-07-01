clear all
close all

%-----load data------
filename1 = './PositionWiseKDEstimates.csv';
filename2 = './PositionWiseKDDifferences.csv';
filename3 = './Data_CoVariation/RefSeq.txt';

FirstPosition = 220;
LastPosition = 340;
AllPositions = 1:535;
NrAllPositions = length(AllPositions);

M1 = dlmread(filename1,',',1,0);
M2 = dlmread(filename2,',',1,0);
M3 = dlmread(filename3);
RefSeq = M3(:,2);

startIdx = find(M1(:,1)==FirstPosition);
endIdx = find(M1(:,1)==LastPosition);
M1 = M1(startIdx:endIdx,:);

startIdx = find(M2(:,1)==FirstPosition);
endIdx = find(M2(:,1)==LastPosition);
if isempty(startIdx) || isempty(endIdx)
    disp('Defined Range exceeds the Input Table')
    return
end
M2 = M2(startIdx:endIdx,:);

Positions = M1(:,1);
Wildtypes = M1(:,2);
NrPositions = length(Positions);
%---------------------
bases = 1:4;                                
StructureDestroyingMutations = [2 4; % A:U -> C:U or U:U
                                1 3; % C:G -> A:G or G:G
                                2 4; % G:U -> C:U or U:U; G:C -> X (all destroying)
                                1 3];% U:G -> A:G or G:G; U:A -> X (all destroying)
StructurePreservingMutations = [3; % A:U -> G:U
                                4; % C:G -> U:G
                                1; % G:U -> A:U
                                2];% U:G -> C:G
PotentialBindingPartner =      [4 nan; %A:U
                                3 nan; %C:G
                                4 2; %G:U, G:C
                                3 1];%U:G, U:A
TestsPreservingVsDestroying = [2 6;%G:U; A:U < A:C & A:U < U:U
                               1 8;%U:G; C:G < A:G & C:G < G:G
                               7 12;%A:U; G:U < C:U & G:U < U:U
                               5 11]+2;%C:G; U:G < A:G & U:G < G:G
                            
OutputMatrix = nan(NrAllPositions,4);%pos, wt, green/yellow/blue, potential partner 
OutputMatrix(:,1) = AllPositions;
OutputMatrix(:,2) = RefSeq;
OutputMatrix(:,3) = 0;

 for i = 1:NrPositions
     wt = Wildtypes(i);
     mbases = bases(bases~=wt);
     pos = Positions(i);
     %OutputMatrix(i,1) = pos;
     %OutputMatrix(i,2) = wt;
     m1 = StructureDestroyingMutations(wt,1);
     m2 = StructureDestroyingMutations(wt,2);
     m3 = StructurePreservingMutations(wt);
     Kd1 = M1(i,(m1-1)*7 + 4);
     p1 = M1(i,(m1-1)*7 + 5);
     Kd2 = M1(i,(m2-1)*7 + 4);
     p2 =  M1(i,(m2-1)*7 + 5);
     Kd3 = M1(i,(m3-1)*7 + 4);
     p3 =  M1(i,(m3-1)*7 + 5);
     %1) 
      % a) structure destroying mutations significantly increase Kd.
    if Kd1 > 1 && p1 < 0.05 && Kd2 >1 && p2 < 0.05 && (wt == 1 || wt == 2)
         p4 = M2(i,TestsPreservingVsDestroying(m3,1));
         p5 = M2(i,TestsPreservingVsDestroying(m3,2));
         % Test if structure, or structure and sequence matter
         % i.e. if WC base
         if p4 < 0.05 &&  p5 < 0.05  % (Watson-Crick pairs A:U, C:G)
             OutputMatrix(pos,3) = 1;
             %OutputMatrix(i,4) = nan;
             OutputMatrix(pos,4) = PotentialBindingPartner(wt,1);
         elseif p4 < 0.05 ||  p5 < 0.05 % (Watson-Crick pairs A:U, C:G)
             %OutputMatrix(i,3) = nan;
             OutputMatrix(pos,3) = 2;
             OutputMatrix(pos,4) = PotentialBindingPartner(wt,1);
             % all mutations matter (seq. specific)
         elseif p3 < 0.05 && Kd3 >1
             OutputMatrix(pos,3) = 5;
         end 
         % b) structure destroying mutations significantly increase Kd. and
         % transition to WC improves binding
         % (Wobble pairs G:U, U:G)
     elseif Kd1 > 1 && p1 < 0.05 && Kd2 >1 && p2 < 0.05 && Kd3 < 1 && (wt == 3 || wt == 4)% && p3 < 0.05 % && ((p3 >= 0.05 && Kd3 >1)|| Kd3 < 1) && (wt == 3 || wt == 4)
         p4 = M2(i,TestsPreservingVsDestroying(m3,1));
         p5 = M2(i,TestsPreservingVsDestroying(m3,2));
         if p4 < 0.05 &&  p5 < 0.05
             OutputMatrix(pos,3) = 3;
             %OutputMatrix(i,4) = nan;
             OutputMatrix(pos,4) = PotentialBindingPartner(wt,1);
         elseif p4 < 0.05 ||  p5 < 0.05
             %OutputMatrix(i,3) = nan;
             OutputMatrix(pos,3) = 3;
             OutputMatrix(pos,4) = PotentialBindingPartner(wt,1);
         end 
         % G or U base and all mutations increase Kd (Watson-Crick G:C or U:A or seq. specific)
      elseif Kd1 > 1 && p1 < 0.05 && Kd2 >1 && p2 < 0.05 && p3 < 0.05 && Kd3 >1 && (wt == 3 || wt == 4)
             OutputMatrix(pos,3) = 4;
         % two mutations including the structure preserving one increase Kd
         % (seq. matters?)
      elseif Kd1 > 1 && Kd2 > 1 && Kd3 > 1 && ((p1 < 0.05  && p3 < 0.05)|| (p2 < 0.05 && p3 < 0.05) || p2 < 0.05 && p1 < 0.05)
            OutputMatrix(pos,3) = 6;
%       elseif (Kd1 > 1 && p1 < 0.05  && p3 < 0.05 && Kd3 >1)|| (Kd2 >1 && p2 < 0.05 && p3 < 0.05 && Kd3 >1)
%             OutputMatrix(pos,3) = 6;
     end   
     
 end
 
filename = strcat('./Figure4_Data.csv');
try,
    delete(filename);
catch
end
%generate Header
Header = {'Pos.';'wt Base';'0 = nothing 1 = green (WC strong) 2= yellow (WC+seq.) 3 = blue (wobble) 4 = red (G or U WC or seq. important [strong]) 5 = pink (A or C seq. important [strong]) 6 = purple all mutations increase Kd but only 2 significantly';'potential binding partner1'};%;;'potential binding partner2'};
fid = fopen(filename,'w');
for i = 1:length(Header)-1
    fprintf(fid,'%s,',char(Header(i)));
end
fprintf(fid,'%s\n',char(Header(end)));
fclose(fid);
%write Output values
dlmwrite(filename,OutputMatrix,'-append')

