%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               TUNE THE PARAMETER D FOR OPTIMAL ORTHONORMAL SYSTEM FOR DISCRIMINANT 
%               ANALYSIS
%
%               This example is performed on data batched in two groups (c1 and c2). 
%               OOS-DA modeling will be performed using d ranging from 1-50 and a 
%               normalised distance score between the groups is calculated for each 
%               iteration. The loss in distance between each step is then further 
%               calculated. If the loss is strictly less than 0.01, the current d is 
%               picked. Two plot are generated to visualise the change in score and 
%               variance between the groups.
%
%               NOTE: The samples need to be sorted batch/group-wise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear()
close()
[Y,drugID,raw]  = xlsread('C:\path\to\data.xlsx');




[loadings,score,eigenvalues] = pca(Y','Centered',false);
no = length(drugID);
m=mean(Y')';

% REDUCE DATA TO PREVENT SINGULARITY
for (i=1:no)
    X(:,i)=loadings(:,1:no-2)'*(Y(:,i)-m); % <-- reduce further if singularity is encountered during OOS-DA
end

% SET NUMBER OF MEMEBERS IN EACH CLASS/BATCH
c1 = 55;
c2 = 72;

classIndices=[ones(1,c1),2*ones(1,c2)];
between_variance = zeros(1,50);
scores = zeros(1,50);
p = ProgressBar(50);
% PERFORM OOS-DA MODELING WITH D FROM 1 TO 50
for (j=1:50)
    p.progress;
    % PERFORM THE OPTIMAL ORTHONORMAL SYSTEM FOR DISCRIMINANT ANALYSIS
    displayFlag=0;
    d=j;
    
    [BasisVectorsAsColumns,ScoreValuesAsColumns,rank] = oosDA(X,classIndices,d,displayFlag);
    
    % FIND INDICES OF QC (MEBENDAZOLE) SAMPLES
    QCindex = strfind(drugID,'Meb');
    QCindex = find(~cellfun(@isempty,QCindex));
    
    G1 = QCindex(QCindex<c1); 
    G2 = QCindex(QCindex>c1);

    % REMOVE BIAS EFFECTS

    Xreduced=X-BasisVectorsAsColumns*ScoreValuesAsColumns;

    mG1 = mean(Xreduced(:,G1)');
    mG2 = mean(Xreduced(:,G2)');

    m = mean([mG1;mG2]);

    b = norm(mG1-m)+norm(mG2-m);
    b = b/2;
    between_variance(j) = b;
    
    wG1 = 0;
    for (i=1:length(G1)) 
        wG1 = wG1 + norm(Xreduced(:,G1(i))-mG1'); 
    end
    wG1 = wG1/length(G1);

    wG2 = 0;
    for (i=1:length(G2)) 
        wG2 = wG2 + norm(Xreduced(:,G2(i))-mG2'); 
    end
    wG2 = wG2/length(G2);
    
    w = (wG1+wG2)/2;
    scores(j) = b/w;
end

figure
plot(cumsum(ones(1,50)),scores); 
ylabel('score')

figure
plot(cumsum(ones(1,50)),between_variance);
ylabel('b')
xlabel('d')

% FIND PARAMETER D 

for (i=1:50)
    t = abs(scores(i+1)-scores(i));
    gain = t/scores(i);
    if (gain < 0.01) 
        fprintf('Dimensionality for OOS-DA modeling: %d \n',i)
        break
    end
end


