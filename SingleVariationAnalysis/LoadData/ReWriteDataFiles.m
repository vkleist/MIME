
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
%           AA   AC   AG   AT   CA   CC   CG   CT   GA   GC   GG   GT   TA   TC   TG   TT
markers = {'ob';'or';'ok';'pb';'or';'og';'pr';'om';'ok';'pr';'oc';'sb';'pb';'om';'sb';'oy'};
FaceColors = {'b';'r';'k';'b';'r';'g';'r';'m';'k';'r';'b';'b';'b';'m';'b';'y'};
%markers = {'ob';'or';'ok';'og';'sb';'sr';'sk';'sg';'^b';'^r';'^k';'^g';'vb';'vr';'vk';'vg'};

disp('Read RefSeq..')
Tmp = dlmread('./RefSeq.txt');
RefSeq = Tmp(:,2);

K = [1,3,5:20];
%Experiment = 1;
for counter = 1:1
    Experiment = counter;
    mkdir(strcat('./',num2str(Experiment)));
    disp(strcat('---Experiment: ',num2str(Experiment),' ---'))
    folder = '';
    filename1 = strcat('./',folder,'tg',num2str(Experiment),'.txt');

    fid = fopen(filename1);
    C = textscan(fid, '%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u','HeaderLines',1);
    fclose(fid);

    NrCollums = length(C);
    NrEntries = length(C{1});

    M = zeros(NrEntries,NrCollums+2);


    for j = [1:18]
        M(:,K(j)) = C{j};
    end
    
    %clear residues > 535
    I = M(:,1)<=length(RefSeq);
    M = M(I,:);
    
    %clear residues > 535
    I = M(:,3)<=length(RefSeq);
    M = M(I,:);
    
    M(:,2) = RefSeq(M(:,1));
    M(:,4) = RefSeq(M(:,3));
    
   [NrEntries,NrCollums] = size(M);
    
    %% Generate new Files and assign Data

    LastPositions = [find(M(2:end,1) - M(1:end-1,1)); NrEntries];
    FirstPositions = [1;LastPositions(1:end-1)+1];

    Position_name = M(LastPositions,1);
    NrPositions = length(LastPositions);

    for i = 1:NrPositions
        disp(strcat('Position: ',num2str(Position_name(i)),' ..'))
        filename_i = strcat('./',num2str(Experiment),'/',num2str(Experiment),'_',num2str(Position_name(i)),'.txt');
        M1 = M(FirstPositions(i):LastPositions(i),:);
        dlmwrite(filename_i,M1,'-append');
        [n,m] = size(M1);
        for j = 1:n
            filename_j = strcat('./',num2str(Experiment),'/',num2str(Experiment),'_',num2str(M1(j,3)),'.txt');
            M2 = M1(j,:);
            dlmwrite(filename_j,M2,'-append');
        end
    end
end

