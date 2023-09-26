%This script computes Procrustes to estimate motor distance between subjects.
%Then, Mantel test is used to compare the solution obtained with different trials used as references.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A1=This is the type of action made by A1 (i.e. the confederate). 
%   Allowed values are 1 for precision grip, 2 for full hand grip

%A2=This is the type of action made by A2 (i.e. the subject). 
%   Allowed values are 1 for precision grip, 2 for full hand grip

%session=This is the kind of task performed.
%   Allowed values are 0 for individual condition (for A2 only)
%                      1 for parallel condition
%                      2 for joint condition
%                      3 for joint+objective condition

%r1,r2,r3,r4=These are the trials to be considered as references for the
%            Procrustes rotation.

%The output 'Mantel procrustes' contains the results of the Mantel test
%among the various references for A1 and A2.

%Example:
%p=mantel_procrustes_paper_script(1,1,1,1,5,10,13)
%   This is to consider precision grip for both A1 and A2 in parallel
%   condition, comparing trial 1,5,10 and 15 as references.

%Coded by Jordi Manuello (M'N'B Lab, University of Torino) in November 2022.
%This code is distributed under the terms of the GNU General Public License
%GPL-3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function p=mantel_procrustes_paper_script(A1,A2,session,r1,r2,r3,r4)

reference=[r1 r2 r3 r4];
%First load the dataset
load '/your/path/to/AllDataset.mat';

%Remove practice trials
k=find(allData.TrialType==0);
allData(k,:)=[];

%Remove actions not congruent with the instruction
k=find(allData.A2_actionACC==0);
allData(k,:)=[];

%Select session and type of prehension
k=find(allData.A1_actiontype==A1 & allData.A2_actiontype==A2 & allData.Session==session); 
allData=allData(k,:);

%Count how many subjects are in the database
subj=unique(allData.SubjID);

%Count how many trials exists for each subject
for i=1:length(subj)
    n_trial(i,1)=length(find(allData.SubjID==subj(i)));
end

%Find the minimum number of trials per subject
min_trial=min(n_trial);

%Check there are enough trials for the selected reference
if any(reference>min_trial)
    disp(strcat('One or more selected references exceed the number of available trials. The maximum value allowed for input "r*" is ',num2str(min_trial)))
    p=0;
    return
end

m_dist_A1_ref=cell(length(reference),1); %This is to store the procustes distance for each reference
m_dist_A2_ref=cell(length(reference),1); %This is to store the procustes distance for each reference

for rr=1:length(reference) %This is to manage different references
   
    %Now work on A1
    Z_all_A1=zeros(10,3,length(subj)); %This is to initialize the matrix storing the last transformation for each subject
    
    for i=1:length(subj)
        [tr_r,tr_c,tr_v]=find(allData.SubjID==subj(i));
        
        %Take  velocity, acceleration and jerk for reference trial
        v_s1=allData(tr_r(reference(rr)),127:136);
        a_s1=allData(tr_r(reference(rr)),137:146);
        j_s1=allData(tr_r(reference(rr)),147:156);
        
        %Build the matrix for s1 that is the reference trial
        s1=[table2array(v_s1)',table2array(a_s1)',table2array(j_s1)'];
        
        Z_trial=zeros(10,3,length(tr_r)); %This is to initialize the matrix storing the transformation for each trial
        
        for j=1:length(tr_r)
            if j~=reference(rr)
                %Take  velocity, acceleration and jerk for comparison trial (s2)
                v_s2=allData(tr_r(j),127:136);
                a_s2=allData(tr_r(j),137:146);
                j_s2=allData(tr_r(j),147:156);
                
                %Build the matrix for s2 that is the comparison trial
                s2=[table2array(v_s2)',table2array(a_s2)',table2array(j_s2)'];
                
                %Compute procrustes for s1,s2
                [d_trial(i,j),Z]=procrustes(s1,s2);
                
                Z_trial(:,:,j)=Z; %This retains each transformation for the subject
            end
        end
        
        %Now average across the rotated trials, after removing the first page
        %that is all zeros
        Z_trial(:,:,1)=[];
        tr_mean=mean(Z_trial,3);
        
        %Store the mean trial
        Z_all_A1(:,:,i)=tr_mean;
    end
    
    %Now compute procustes between subjects and build the distance matrix
    
    m_dist_A1=zeros(length(subj),length(subj)); %This creates the distance matrix
    sz=[length(subj),length(subj)]; %This must match the size of the distance matrix
    
    for i=1:numel(m_dist_A1)
        [r,c]=ind2sub(sz,i);
        s1=Z_all_A1(:,:,r);
        s2=Z_all_A1(:,:,c);
        [m_dist_A1(i),Z]=procrustes(s1,s2);
    end
    clims=[0 1];
    f1=imagesc(m_dist_A1,clims);
    colorbar
    title(strcat('procrustes distance A1; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),'Ref',num2str(reference(rr))))
    saveas(f1,strcat('procrustes distance A1; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),'Ref',num2str(reference(rr)),'.jpg'));
    
    m_dist_A1_ref{rr}=m_dist_A1;
    clear d_trial
    
    %Repeate everything on A2
    Z_all_A2=zeros(10,3,length(subj)); %This is to initialize the matrix storing the last transformation for each subject
    
    for i=1:length(subj)
        [tr_r,tr_c,tr_v]=find(allData.SubjID==subj(i));
        
        %Take  velocity, acceleration and jerk for reference trial
        v_s1=allData(tr_r(reference(rr)),31:40);
        a_s1=allData(tr_r(reference(rr)),41:50);
        j_s1=allData(tr_r(reference(rr)),51:60);
        
        %Build the matrix for s1 that is the reference trial
        s1=[table2array(v_s1)',table2array(a_s1)',table2array(j_s1)'];
        
        Z_trial=zeros(10,3,length(tr_r)); %This is to initialize the matrix storing the transformation for each trial
        
        for j=1:length(tr_r)
            if j~=reference(rr)
                %Take  velocity, acceleration and jerk for comparison trial (s2)
                v_s2=allData(tr_r(j),31:40);
                a_s2=allData(tr_r(j),41:50);
                j_s2=allData(tr_r(j),51:60);
                
                %Build the matrix for s2 that is the comparison trial
                s2=[table2array(v_s2)',table2array(a_s2)',table2array(j_s2)'];
                
                %Compute procrustes for s1,s2
                [d_trial(i,j),Z]=procrustes(s1,s2);
                
                Z_trial(:,:,j)=Z; %This retains each transformation for the subject
            end
        end
        
        %Now average across the rotated trials, after removing the first page
        %that is all zeros
        Z_trial(:,:,1)=[];
        tr_mean=mean(Z_trial,3);
        
        %Store the mean trial
        Z_all_A2(:,:,i)=tr_mean;
    end
    
    %Now compute procustes between subjects and build the distance matrix
    
    m_dist_A2=zeros(length(subj),length(subj)); %This creates the distance matrix
    sz=[length(subj),length(subj)]; %This must match the size of the distance matrix
    
    for i=1:numel(m_dist_A2)
        [r,c]=ind2sub(sz,i);
        s1=Z_all_A2(:,:,r);
        s2=Z_all_A2(:,:,c);
        [m_dist_A2(i),Z]=procrustes(s1,s2);
    end
    clims=[0 1];
    f1=imagesc(m_dist_A2,clims);
    colorbar
    title(strcat('procrustes distance A2; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),'Ref',num2str(reference(rr))))
    saveas(f1,strcat('procrustes distance A2; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),'Ref',num2str(reference(rr)),'.jpg'));
    
    m_dist_A2_ref{rr}=m_dist_A2;
end

%Now prepare for Mantel on A1

%Make 0 in the diagonal
m_dist1_0 = m_dist_A1_ref{1,1} - diag(diag(m_dist_A1_ref{1,1}));
m_dist2_0 = m_dist_A1_ref{2,1} - diag(diag(m_dist_A1_ref{2,1}));
m_dist3_0 = m_dist_A1_ref{3,1} - diag(diag(m_dist_A1_ref{3,1}));
m_dist4_0 = m_dist_A1_ref{4,1} - diag(diag(m_dist_A1_ref{4,1}));

%Build the permutation matrix
m_perm=zeros(4,4);
p_perm=zeros(4,4);
sz=[size(m_perm,1),size(m_perm,1)];

%Concatenate the examples
m_dist_all=cat(3,m_dist1_0,m_dist2_0,m_dist3_0,m_dist4_0);

%Run Mantel
for i=1:numel(m_perm)
    [r,c]=ind2sub(sz,i);
    [m_perm(i) p_perm(i)]=bramila_mantel(m_dist_all(:,:,r),m_dist_all(:,:,c),5000,'pearson');
end

m_perm_A1=m_perm;
p_perm_A1=p_perm;

clear m_perm;
clear p_perm;

%Now prepare for Mantel on A2

%Make 0 in the diagonal
m_dist1_0 = m_dist_A2_ref{1,1} - diag(diag(m_dist_A2_ref{1,1}));
m_dist2_0 = m_dist_A2_ref{2,1} - diag(diag(m_dist_A2_ref{2,1}));
m_dist3_0 = m_dist_A2_ref{3,1} - diag(diag(m_dist_A2_ref{3,1}));
m_dist4_0 = m_dist_A2_ref{4,1} - diag(diag(m_dist_A2_ref{4,1}));

%Build the permutation matrix
m_perm=zeros(4,4);
p_perm=zeros(4,4);
sz=[size(m_perm,1),size(m_perm,1)];

%Concatenate the examples
m_dist_all=cat(3,m_dist1_0,m_dist2_0,m_dist3_0,m_dist4_0);

%Run Mantel
for i=1:numel(m_perm)
    [r,c]=ind2sub(sz,i);
    [m_perm(i) p_perm(i)]=bramila_mantel(m_dist_all(:,:,r),m_dist_all(:,:,c),5000,'pearson');
end

m_perm_A2=m_perm;
p_perm_A2=p_perm;

save(strcat('Mantel procrustes; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),' Ref',num2str(reference(1)),';',num2str(reference(2)),';',num2str(reference(3)),';',num2str(reference(4))),'m_perm_A1','p_perm_A1','m_perm_A2','p_perm_A2');
p=1;
return
