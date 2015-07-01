
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

close all
clear all
%with Watson-Crick and Wobble base- vs. non Watson-Crick base pairing
%           AAA   AAC   AAG   AAT   ACA   ACC   ACG   ACT  ....   TTA   TTC   TTG   TTT

%% Make transition Rules for expanding (index changes)
% Transitionsrules for entries: 
% Rule 1&2 are for completing/expanding where pos 1 < pos 2 hold

% Rule 3 is for creating the mirror matrix (lower left triangular), where
% pos 1 > pos 2 holds

% Rule 1: [x1 x2 x3] -> [x1 x3 x2]
% Rule 2: [x1 x2 x3] -> [x2 x3 x1]
% Rule 3: [y1 y2 y3] -> [y2 y1 y3]

TransitionRule1 =  zeros(1,70);
TransitionRule1(1:6) = [1 2 5 6 3 4];

TransitionRule2 =  zeros(1,70);
TransitionRule2(1:6) = [5 6 1 2 3 4];

TransitionRule3 =  zeros(1,70);
TransitionRule3(1:6) = [3 4 1 2 5 6];

for i = 1:4
    for j = 1:4
        for z = 1:4
            oldidx = (i-1)*16 + (j-1)*4 + z + 6;
            %Rule 1
            newidx = (i-1)*16 + (z-1)*4 + j + 6;
            TransitionRule1(oldidx) = newidx;
            %Rule 2
            newidx = (j-1)*16 + (z-1)*4 + i + 6;
            TransitionRule2(oldidx) = newidx;
            %Rule 3
            newidx = (j-1)*16 + (i-1)*4 + z + 6;
            TransitionRule3(oldidx) = newidx;           
        end
    end
end

%% Start File reading and Expansion (generates new files)

disp('Read RefSeq..')
Tmp = dlmread('../Data_CoVariation/RefSeq.txt');
RefSeq = Tmp(:,2);

FirstPosition = 200;

for counter = 8:9 % over all experiments
    % a) ----READ DATA ----
    %read data for experiment 'counter'
    Experiment = counter;
    mkdir(strcat('./',num2str(Experiment)));
    disp(strcat('---Experiment: ',num2str(Experiment),' ---'))
    folder = '';
    filename1 = strcat('./',folder,num2str(Experiment),'.txt');

    fid = fopen(filename1);
    C = textscan(fid, '%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u','HeaderLines',1);
    NrCollums = length(C);
    NrEntries = length(C{1});
    M = zeros(NrEntries,NrCollums+3);
    %Assign nt Position and wt Sequence for positions
    M(:,1) = C{1}(:,1)+FirstPosition;%pos 1;
    M(:,3) = C{2}(:,1)+FirstPosition;%pos 2;
    M(:,5) = C{3}(:,1)+FirstPosition;%pos 3;
    M(:,2) = RefSeq(M(:,1));% nt of pos 1
    M(:,4) = RefSeq(M(:,3));% nt of pos 2
    M(:,6) = RefSeq(M(:,5));% nt of pos 3
    %Assign contents (nt counts)
    for j = 4:67;
        M(1:NrEntries,j+3) = C{j}(:,1);
    end
     
    fclose(fid);
    % b) ----EXPAND----
    
    %detect when pos 1 changes + add the last position (i.e. the last idx when pos1 is still the same)
    LastPositions = [find(M(2:end,1) - M(1:end-1,1)); NrEntries]; 
    %the first idx of pos1 in M
    FirstPositions = [1;LastPositions(1:end-1)+1];

    Position_name1 = M(LastPositions,1);
    NrPositions1 = length(LastPositions);

    for i = 1:NrPositions1% over all positions1 1...x

        %disp(strcat('Position1: ',num2str(Position_name1(i)),' ..'))
        M1 = M(FirstPositions(i):LastPositions(i),:); % all with pos1 = X, pos2 >X, pos3 > pos2

        NrEntries2 = LastPositions(i)-FirstPositions(i)+1;
        %detect when pos 2 changes + add the last position (i.e. the last idx when pos1 is still the same)
        LastPositions2 = [find(M1(2:end,3) - M1(1:end-1,3)); NrEntries2] ;
        %the first idx of pos2
        FirstPositions2 = [1;LastPositions2(1:end-1)+1];

        Position_name2 = M1(LastPositions2,3);
        NrPositions2 = length(LastPositions2);
         for j = 1:NrPositions2  %over all positions 1...y| pos1 = x   
            %write the ordered data pos1 < pos2 < pos3
            disp(strcat('Position1: ',num2str(Position_name1(i)),'Position2: ',num2str(Position_name2(j)),' ..'))
            filename_i = strcat('./',num2str(Experiment),'/',num2str(Position_name1(i)),'_',num2str(Position_name2(j)),'.txt');
            M2 = M1(FirstPositions2(j):LastPositions2(j),:); %all with pos2 = Y | pos1 = X
            dlmwrite(filename_i,M2,'-append');
            n = LastPositions2(j)-FirstPositions2(j)+1;
            for z = 1:n
                %permutation rule 1<--------also permute entries-----------
                Perm1 = zeros(1,70);
                Perm1(TransitionRule1) = M2(z,:);
                filename_z = strcat('./',num2str(Experiment),'/',num2str(Perm1(1)),'_',num2str(Perm1(3)),'.txt');
                dlmwrite(filename_z,Perm1,'-append');
                clear Perm1
                %permutation rule 2 <--------also permute
                %entries-----------
                Perm2 = zeros(1,70);
                Perm2(TransitionRule2) = M2(z,:);
                filename_z = strcat('./',num2str(Experiment),'/',num2str(Perm2(1)),'_',num2str(Perm2(3)),'.txt');
                dlmwrite(filename_z,Perm2,'-append');
            end
         end
    end   
end

