clear -regexp ^[^D]([\w]+)?;
close all; tic;
addpath(genpath('\\argon-cifs\chembio_datasets\csdev\IMM'));
wf = '\\argon-cifs\chembio_datasets\csdev\IMM\Data\SIDER';
cd(wf);

% read SE and PTs, make matrices
% rerun with uniprot names converted to gene name - lost 5 pubchemIDs
%DS = DFread('pubchem_se_siderPT_DB.txt',wf,[],[],true);
DS = DFread('pubchem_se_siderPT_DB_geneName.txt',wf,[],[],true);
[xpDS,xvDS,xiDS] = DFindex(DS,{'pubchem_ID','side_effect'});
mDS = sparse(xiDS(:,1),xiDS(:,2),true,numel(xvDS.pubchem_ID),numel(xvDS.side_effect));
%DT = DFread('drug_target_pairings.txt',wf,[],[],true);
DT = DFread('drug_target_pairings_geneName.txt',wf,[],[],true);
%[xpDT,xvDT,xiDT] = DFindex(DT,{'pubchem_ID','uniprot_target'});
[xpDT,xvDT,xiDT] = DFindex(DT,{'pubchem_ID','target_gene_name'});
%mDT = sparse(xiDT(:,1),xiDT(:,2),true,numel(xvDT.pubchem_ID),numel(xvDT.uniprot_target));
mDT = sparse(xiDT(:,1),xiDT(:,2),true,numel(xvDT.pubchem_ID),numel(xvDT.target_gene_name));

xx = sum(mDS,1)./sum(sum(mDS,1));
xy = sum(mDS,2)./sum(sum(mDS,2));
x = xy*xx;
x0 = max(1/numel(mDS),mDS);
ef = log10(x0)-log10(x);
pp = pdist(ef,'cityblock');
mDlnk = linkage(pp,'complete');
figure();

% to make heatmap have better colors
c1 = make1cmap(1, false);

[~,~,mDprm] = dendrogram(mDlnk,0);
mDsqf = squareform(pp);
mDsqf = mDsqf(mDprm,mDprm);
figure();
imagesc(mDsqf);

mDclid = clinkid(mDlnk);
mDobj = obj2node(mDclid);
[Nlod,Npv,Nn1,Nn2,Npur,Ncnf,Nov] = nodenfast(mDobj,mDT);

% new stuff 

% filtering based on p-values first 
% f = (Npv<0.05); % gets matrix of p-values significant according
% to this value 
% makes matrix of all NaN the size of the Npv matrix that holds p-values 
Nqv = nan(size(Npv)); 

% pass this function DataMatrix object containing p-values for each feature in a data set
% Estimate false discovery rate (FDR) for multiple hypothesis testing
q = mafdr(Npv(:));
% this copies what we have in q back into Nqv 
Nqv(:) = q;

% now that we have prepared and gotten what we want - the false discovery
% rate for each of these p-values 
% now we filter based on the p-values taking only false discovery rates
% less than 0.25 and purity and confidence greater than 0.5 and the overlap greater than 1  
f = (Nqv<0.25&Npur>=0.5&Ncnf>=0.5&Nov>1);
nnz(f); % this tells us how many values 
% we get 27! 

[out.fx,out.fy]=find(f);  % this gives us the 25 row and columns numbers 
out.q = Nqv(f); % this goes and gives us the actual Q values for these 27
out.p = Npv(f); % actual p values for the 27
out.pur = Npur(f); % actual purity values for the 27
out.cnf = Ncnf(f); % actual confidence values for the 27
out.ov = Nov(f);  % actual overlap values for the 27

% sorting based on Q : smalled values first --> most significant 
[~,sort_order] = sort(out.q,'ascend');
out = DFkeeprow(out, sort_order) ; % use this to sort and filter 
% keeprow so it keeps all the rows and sorts it 

% now we want to get the values that sit at the row and column values
% need to make for loop to loop through all of these 
%mDobj(:,out.fx(1)); % gives the number of compound IDs that correspond to
% the specific column(cluster) 

% this gives the list of pubchemIDs from the indices we got for result 1
%xvDS.pubchem_ID(mDobj(:,out.fx(1)));
%xvDT.uniprot_target(out.fy(1)); % protein ID for result 1
%xvDT.target_gene_name(out.fy(1)); % gene name for result 1 (just one)
% list of SE in 1/2 or more of the set of result 1 
%xvDS.side_effect(mean(mDS(mDobj(:,out.fx(1)),:))>=0.5); 

% make empty df with fields we want ot have in it
outputs_27 = struct('groupID', [], 'target', [], 'cpd', []);
% make a for loop to go through 27 outputs 
clusters = 27;
for c = 1:clusters
    temp = xvDS.pubchem_ID(mDobj(:,out.fx(c)));
    num = numel(temp); % number of elements (pubchemID)
    outputs_27.cpd = [outputs_27.cpd;temp]; % semicolon means append down and space means append across 
    outputs_27.groupID = [outputs_27.groupID;repmat(c, num, 1)]; % column vector which is 4 copies of c
    outputs_27.target = [outputs_27.target;repmat(xvDT.target_gene_name(out.fy(c)), num, 1)];
end

DFwrite(outputs_27, 'outputs_27_2018_04_22.txt');

% do the same thing to get the list of SEs 

outputs_27_SE = struct('groupID', [], 'SE', []);
for c = 1:clusters
    temp = xvDS.side_effect(mean(mDS(mDobj(:,out.fx(c)),:))>=0.5); 
    num = numel(temp);
    outputs_27_SE.SE = [outputs_27_SE.SE;temp];
    outputs_27_SE.groupID = [outputs_27_SE.groupID;repmat(c, num, 1)]; 
end

DFwrite(outputs_27_SE, 'outputs_27_SE_2018_04_22.txt');

toc;

% notes

% will have nested trees of significance, among the 77
% we think there may be 25 different results  

% pruning tree: if parents has better p-value than child ,keep that 
% want to keep the bedst think in those non-independent paths
% we have a child that is significant and then a parent may or may not have
% the same target but could still be significant, so we will get a lot of
% non-independen results 

% need to figure out hwo to tabulate the 25 outputs 
% 3 outputs: protein & its biology, 2+ compounds, what caused them to be
% called similar in the first plce
% we clustered based on SE but that can go in different ways 
% small # of SE or could be a lot of SE in common --> need to know what
% these are and think about how to present this information
% investiage string concatenation functions on my own and figure out how
% best to do it 

% think about interpretation of purity and confidence 

% Hypothesis is that the target of a compound is causing its SE or soething
% like this 

% now we are looking at the 
% imaein recording the set of SE that represent the union across all these
% compounds or the intersection or anything in between 
% 

% we now have a complete path to ask questions 
% we have these 3 outputs for each of these 25 significant findings 
% since we see so many SE in common, this tells us that maybe when we tried
% to account for the fact that cpds have many SE, they are catgorized as
% similary, event though we tried to account for htis
% might see many of the same compounds or SE 
% how do we export these in a way that we study these individually and
% study them as a group as well
% do e see trends emerging, study them as a group 
% 

% put this into R 
