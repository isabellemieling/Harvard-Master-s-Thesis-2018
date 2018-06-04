clear -regexp ^[^D]([\w]+)?;
close all; tic;
addpath(genpath('\\argon-cifs\chembio_datasets\csdev\IMM'));
wf = '\\argon-cifs\chembio_datasets\csdev\IMM\Data\SIDER';
cd(wf);

% read in SE and targets and compounds
DS = DFread('pubchem_se_siderPT_DB_geneName.txt',wf,[],[],true);
[xpDS,xvDS,xiDS] = DFindex(DS,{'pubchem_ID','side_effect'});
mDS = sparse(xiDS(:,1),xiDS(:,2),true,numel(xvDS.pubchem_ID),numel(xvDS.side_effect));

DT = DFread('drug_target_pairings_geneName.txt',wf,[],[],true);
[xpDT,xvDT,xiDT] = DFindex(DT,{'pubchem_ID','target_gene_name'});
mDT = sparse(xiDT(:,1),xiDT(:,2),true,numel(xvDT.pubchem_ID),numel(xvDT.target_gene_name));

xx = sum(mDS,1)./sum(sum(mDS,1));
xy = sum(mDS,2)./sum(sum(mDS,2));
x = xy*xx;
x0 = max(1/numel(mDS),mDS);
ef = log10(x0)-log10(x);
pp = pdist(ef,'cityblock');
mDlnk = linkage(pp,'complete');

filename = '\\argon-cifs\chembio_datasets\csdev\IMM\Data\Final_dataset\pubchem_names.txt';
A = tdfread(filename);

mDsqf = squareform(pp);

% redo for each result 

% RESULT 1 
shortlist = [53232 ;54454 ;446156; 4915];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 2 
shortlist = [2479 ;4737 ;5193];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 3 
shortlist = [2662; 3883; 4679; 5029; 119607];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 4 
shortlist = [3715; 3826; 4594; 55891; 5362129; 5484727;];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 5 
shortlist = [4595; 5206; 42113];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 6 
shortlist = [57363; 6918296];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 7 
shortlist = [445643; 446541; 5281078];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 8 
shortlist = [9871419; 9875401; 10182969];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 9 
shortlist = [20279; 119182];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 10 
shortlist = [5245; 60852];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 11 
shortlist = [9444; 451668];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 12 
shortlist = [3333; 3784];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

% RESULT 13 
shortlist = [5743; 5754; 6741];
subset =  ismember(xvDS.pubchem_ID, shortlist);
result = mDsqf(subset, subset);
Result = squareform(result); % do squareform again to get unsquareformed 
resLink = linkage(Result,'complete');

sumbset1 = find(ismember(xvDS.pubchem_ID, shortlist));
labels = A.Name(sumbset1,:);

figure(1);
H = dendrogram(resLink,0,'Orientation', 'right', 'Labels', labels);
set(H, 'LineWidth', 3);
set(gca,'FontSize',35);

