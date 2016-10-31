%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               REMOVE UNWANTED VARIATION (i.e. BATCH EFFECTS) BY USING OPTIMAL 
%               ORTHONORMAL SYSTEM FOR DISCRIMINANT ANALYSIS. 
%
%               This example is performed on data batched in four groups (c1-c4). 
%               The d is set to 10 as default, but this need to be optimixed for 
%               each dataset specifically using the parameter_tuning.m
%
%               NOTE: The samples need to be sorted batch/group-wise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear()
close()

[Y,drugID,raw]  = xlsread('C:\Users\Milos\Documents\Klinisk farm\MCF-10\Set pipeline\corr79\batchsorted.xlsx');
[loadings,score,eigenvalues] = pca(Y','Centered',false);
no = length(drugID);
m = mean(Y')';

% REDUCE DATA TO PREVENT SINGULARITY
for (i=1:no)
X(:,i)=loadings(:,1:no-8)'*(Y(:,i)-m); % <-- reduce further if singularity is encountered during OOS-DA
end


% PERFORM THE OPTIMAL ORTHONORMAL SYSTEM FOR DISCRIMINANT ANALYSIS
displayFlag=0;
d=10; % <----- set d

% SET NUMBER OF MEMEBERS IN EACH CLASS/BATCH
c1 = 29;
c2 = 30;
c3 = 36;
c4 = 32;

classIndices=[ones(1,c1),2*ones(1,c2),3*ones(1,c3),4*ones(1,c4)];
[BasisVectorsAsColumns,ScoreValuesAsColumns,rank] = oosDA(X,classIndices,d,displayFlag);


% REMOVE BATCH EFFECTS OR OTHER UNWANTED VARIATION

Xreduced=X-BasisVectorsAsColumns*ScoreValuesAsColumns;

for (i=1:no)
Yhat(:,i)=m+loadings(:,1:no-8)*Xreduced(:,i); 
end


%xlswrite('C:\Users\Milos\Desktop\MCF-10\Set pipeline\corr79\batcheremoved_d10.xlsx',Yhat)

% VISUALIZE THE CORRECTED DATA IN 3D PLOT
figure
K=d-2;
plot3(ScoreValuesAsColumns(K,:),ScoreValuesAsColumns(K+1,:),ScoreValuesAsColumns(K+2,:),'.w');grid on
hold on
plot3(ScoreValuesAsColumns(K,1:c1),ScoreValuesAsColumns(K+1,1:c1),ScoreValuesAsColumns(K+2,1:c1),'*r')
plot3(ScoreValuesAsColumns(K,c1+1:c1+c2),ScoreValuesAsColumns(K+1,c1+1:c1+c2),ScoreValuesAsColumns(K+2,c1+1:c1+c2),'*b')
plot3(ScoreValuesAsColumns(K,59+1:59+c3),ScoreValuesAsColumns(K+1,59+1:59+c3),ScoreValuesAsColumns(K+2,59+1:59+c3),'*g')
plot3(ScoreValuesAsColumns(K,95+1:95+c4),ScoreValuesAsColumns(K+1,95+1:95+c4),ScoreValuesAsColumns(K+2,95+1:95+c4),'*m')
