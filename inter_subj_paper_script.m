%This script computes Procrustes to estimate motor distance between subjects.
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

%r=These is the trial to be considered as references for the
%            Procrustes rotation.

%The output 'inter-subject distance' contains the inter-subject distance matrix based on Procrustes distance 

%Example:
%p=mantel_procrustes_paper_script(1,1,1,1)
%   This is to consider precision grip for both A1 and A2 in parallel
%   condition, using trial 1 as references.

%Coded by Jordi Manuello (M'N'B Lab, University of Torino) in December 2022.
%This code is distributed under the terms of the GNU General Public License
%GPL-3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=inter_subj_paper_script(A1,A2,session,r)

reference=r;
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
    disp(strcat('The selected references exceed the number of available trials. The maximum value allowed for input "r" is ',num2str(min_trial)))
    p=0;
    return
end

%Now work on A1
    Z_all_A1=zeros(10,3,length(subj)); %This is to initialize the matrix storing the last transformation for each subject
    
    for i=1:length(subj)
        [tr_r,tr_c,tr_v]=find(allData.SubjID==subj(i));
        
        %Take  velocity, acceleration and jerk for reference trial
        v_s1=allData(tr_r(r),127:136);
        a_s1=allData(tr_r(r),137:146);
        j_s1=allData(tr_r(r),147:156);
        
        %Build the matrix for s1 that is the reference trial
        s1=[table2array(v_s1)',table2array(a_s1)',table2array(j_s1)'];
        
        Z_trial=zeros(10,3,length(tr_r)); %This is to initialize the matrix storing the transformation for each trial
        
        for j=1:length(tr_r)
            if j~=r
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
    
    clear d_trial
    
    %Repeate everything on A2
    Z_all_A2=zeros(10,3,length(subj)); %This is to initialize the matrix storing the last transformation for each subject
    
    for i=1:length(subj)
        [tr_r,tr_c,tr_v]=find(allData.SubjID==subj(i));
        
        %Take  velocity, acceleration and jerk for reference trial
        v_s1=allData(tr_r(r),31:40);
        a_s1=allData(tr_r(r),41:50);
        j_s1=allData(tr_r(r),51:60);
        
        %Build the matrix for s1 that is the reference trial
        s1=[table2array(v_s1)',table2array(a_s1)',table2array(j_s1)'];
        
        Z_trial=zeros(10,3,length(tr_r)); %This is to initialize the matrix storing the transformation for each trial
        
        for j=1:length(tr_r)
            if j~=r
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

    %Concatenate the two matrices
Z_A1_A2=cat(3,Z_all_A1,Z_all_A2);

%Now compute procustes between subjects and build the distance matrix

m_dist=zeros(size(Z_A1_A2,3),size(Z_A1_A2,3)); %This creates the distance matrix
sz=[size(Z_A1_A2,3),size(Z_A1_A2,3)]; %This must match the size of the distance matrix

for i=1:numel(m_dist)
    [r,c]=ind2sub(sz,i);
    s1=Z_A1_A2(:,:,r);
    s2=Z_A1_A2(:,:,c);
    [m_dist(i),Z]=procrustes(s1,s2);
end
    
    clims=[0 1];
    f1=imagesc(m_dist,clims);
    colorbar
    title(strcat('inter-subject distance; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),' Ref',num2str(reference)))
    saveas(f1,strcat('inter-subject distance; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),'Ref',num2str(reference),'.jpg'));
    
    save(strcat('inter-subject distance; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),' Ref',num2str(reference)),'m_dist')
    
    %This is to prepare the distance matrix to be used as input for the MDS in Orange
m_dist_orange = m_dist - diag(diag(m_dist));
m_dist_orange=tril(m_dist_orange);

dlmwrite(strcat('inter-subject distance; A1=',num2str(A1),' A2=',num2str(A2),' Session ',num2str(session),'Ref',num2str(reference),'.txt'),m_dist_orange,'delimiter','\t','precision',4)

p=1;
return