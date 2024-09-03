% This script loads the human GEM Human1, integrates the model with
% thermodynamic properties of metabolic reactions predicted by different
% models, and analyzes the relationship between reaction thermodynamics and
% network topology features.
% This script was written by Ziwei Dai (daizw@sustech.edu.cn)

%% Read model from the .mat file
load Human-GEM.mat

% Read real reaction Gibbs free energy from file
prediction_model_names = ["eQuilibrator","dGbyG"];
n_rxns = length(ihuman.rxns);
feature_names = ["Average dGr","Maximal dGr","Minimal dGr", ...
    "Standard dGr", "Reaction type", "Maximal absolute flux"];
n_features = length(feature_names);
thermodynamic_data = cell(1,2);
flux_data = cell(1,2);
for i = 1:2
    filename_real = sprintf("Human1_real_dGrs_distribution_%s.xlsx",prediction_model_names{i});

    % Read fluxes from files
    flux_table = readtable(filename_real,"Sheet","Vs","ReadRowNames",true,...
        "ReadVariableNames",true);
    flux_data{i} = flux_table;

    % Read model-predicted reaction deltaGs from files
    real_dG_table = readtable(filename_real,"Sheet","real dGrs","ReadRowNames", ...
        true,"ReadVariableNames",true);
    average_dGr = mean(real_dG_table{:,:},2);
    max_dGr = max(real_dG_table{:,:},[],2); %Real dGs that are the most "not negative"
    min_dGr = min(real_dG_table{:,:},[],2); %Real dGs that are the most negative
    rxn_type = repmat("NTDR",n_rxns,1);
    max_flux_abs = max(abs(flux_table{:,:}),[],2);
    rxn_type(max_dGr < -200 & max_flux_abs > 0) = "TDR";
    % Read model-predicted standard deltaGs from files
    % Read model-predicted standard deltaGs from files
    filename_standard = sprintf("Human1_standard_dGr_%s.csv",prediction_model_names{i});
    std_dGr_table = readtable(filename_standard,"ReadRowNames",true,"ReadVariableNames",true);
    std_dGr = std_dGr_table{:,1};
    thermodynamic_data{i} = table(average_dGr,max_dGr,min_dGr,std_dGr,...
        rxn_type,max_flux_abs);
    thermodynamic_data{i}.Properties.VariableNames = feature_names;
end

%% Process information about reaction subsystems
rxn_subsystems = cell(n_rxns,1);
for i = 1:n_rxns
    rxn_subsystems{i} = ihuman.subSystems{i}{1};    
end

clear filename_real real_dG_table average_dGr max_dGr min_dGr TDR_labels NTDR_labels

% Read the CCLE protein abundance data and compute average expression level
% and CV/PCA-based heterogeneity scores
% Load proteomics data and compare between TDRs and others
%% Analyze the expression level of genes
expression_matrix_raw = readtable("human_ensembl_expression.csv","ReadRowNames",true);
CCLE_proteomics_matrix = table2array(expression_matrix_raw);
CCLE_proteomics_rxns = ihuman.rxnGeneMat * CCLE_proteomics_matrix;
average_protein_abundance = mean(CCLE_proteomics_rxns,2);
CCLE_proteomics_CV = std(CCLE_proteomics_rxns,[],2)./mean(CCLE_proteomics_rxns,2);
CCLE_proteomics_normalized = zscore(CCLE_proteomics_rxns');
CCLE_proteomics_normalized(isnan(CCLE_proteomics_normalized)) = 0;

[coef,~,~,~,explained,~] = pca(CCLE_proteomics_normalized);
n_pcs_kept = find(cumsum(explained)>90,1);
heterogeneity_score = sum(coef(:,1:n_pcs_kept).^2,2);
heterogeneity_score(average_protein_abundance==0) = NaN;
heterogeneity_score(isnan(average_protein_abundance)) = NaN;
rxn_proteomics_features = table(average_protein_abundance,...
    CCLE_proteomics_CV,heterogeneity_score);

%% Compute network topological features of reactions
degree_mets = sum((ihuman.S~=0),2);
met_hubs = readtable("met_hubs.csv");
met_hub_names = met_hubs{met_hubs{:,3}==1,"Var1"};
hub_mets_idx = find(ismember(ihuman.metNames,met_hub_names));
S_rm_hubs = ihuman.S;
S_rm_hubs(hub_mets_idx,:) = [];
rxn_adjacency = [];
for i = 1:size(S_rm_hubs,1)
    rxn_connected = find(S_rm_hubs(i,:)~=0);
    if length(rxn_connected)>1
        rxn_adjacency = [rxn_adjacency;nchoosek(rxn_connected,2)];
    end
end
rxn_adjacency = unique(rxn_adjacency,'rows');
graph_rxns = graph(rxn_adjacency(:,1),rxn_adjacency(:,2));
degree_rxns = degree(graph_rxns);

% Compute topological metrics for each reaction
closeness = centrality(graph_rxns,'closeness'); 
% Closeness centrality
betweenness = centrality(graph_rxns,'betweenness'); 
% Betweenness centrality
pagerank = centrality(graph_rxns,'pagerank');
% Pagerank centrality
eigenvector = centrality(graph_rxns,'eigenvector');
% Eigenvector centrality

%{
uptake_rxn_list = ["Exchange of D-Glucose", "Exchange of L-Glutamine",...
    "Exchange of L-Homoserine",  "Exchange of L-Isoleucine", ...
    "Exchange of L-Leucine",  "Exchange of L-Lysine", ... 
    "Exchange of L-Methionine",  "Exchange of L-Phenylalanine", ...
    "Exchange of L-Serine", "Exchange of L-Threonine", ...
    "Exchange of L-Tryptophan", "Exchange of L-Tyrosine",...
    "Exchange of L-Valine"];
%}
uptake_rxn_IDs = [7595 7621 7718 7600 7601 7602 7603 7604 7627 7605 7606 7622 7607];
distance_uptake = distances(graph_rxns, uptake_rxn_IDs);
% Shortest path distance to extracellular transport rxns

rxn_topological_features = table(degree_rxns,closeness,betweenness,min(distance_uptake)');
topo_feature_names = {'Degree','Closeness','Betweenness','Distance to nutrient uptake'};
rxn_topological_features.Properties.VariableNames = topo_feature_names;

%% Compare the sets of TDRs and NTDRs predicted by dGbyG and eQuilibrator
rxn_type_eQ = thermodynamic_data{1}.("Reaction type");
rxn_type_dGbyG = thermodynamic_data{2}.("Reaction type");
uniq_rxn_types = ["TDR","NTDR"];
sankey_links = cell(4,3);
n = 0;
for i = 1:2
    for j = 1:2
        n = n + 1;
        node_start = uniq_rxn_types(i);
        node_end = uniq_rxn_types(j);
        count = sum(ismember(rxn_type_dGbyG,node_start) & ismember(rxn_type_eQ,node_end));
        sankey_links{n,1} = convertStringsToChars(strcat(node_start,".dGbyG"));
        sankey_links{n,2} = convertStringsToChars(strcat(node_end,".eQuilibrator"));
        sankey_links{n,3} = count;
    end
end
sankey_links(4,:) = [];
SK=SSankey(sankey_links(:,1),sankey_links(:,2),sankey_links(:,3));
SK.ColorList = [0.9451 0.7412 0.2353;0.3412 0.4980 0.8588;...
    0.6400 0.2949 0.4110;0.2328 0.1871 0.2785];
SK.draw();

%% Extract information about what reactions are artificial, 
% exchange, and transport reactions
idx_exchange_demand = find(ismember(rxn_subsystems,{'Exchange/demand reactions'}));
idx_transport = find(contains(rxn_subsystems,'Transport'));
idx_misc = find(ismember(rxn_subsystems,{'Miscellaneous'}));
idx_true_reaction = setdiff(1:n_rxns,[idx_exchange_demand;idx_transport;idx_misc]);
max_abs_flux = zeros(n_rxns,2);
for i = 1:2
    max_abs_flux(:,i) = max(abs(flux_data{i}{:,:}),[],2);
end
is_nonzero_flux = (max_abs_flux > 0);

%% Train a random forest model to predict TDR/NTDRs
rxn_features_comb = [rxn_proteomics_features rxn_topological_features];
x = table2array(rxn_topological_features);
n_topo_features = size(x,2);
rxn_topo_features_offset = zeros(1,n_topo_features);
for i = 1:n_topo_features
    if min(x(:,i))==0
        rxn_topo_features_offset(i) = min(x(x(:,i)>0,i));
    end
end
x(isinf(x)) = 100; % Distance between two disconnected nodes assumed to be 10
x = x(idx_true_reaction,:);
x_log_transform = log(x + rxn_topo_features_offset);

%x(isnan(x)) = 0;
figure;
hold on;
line_colors = [241 189 60;87 127 219]/255;
for n = 1:2
    x = x_log_transform;
    x(is_nonzero_flux(idx_true_reaction,n)==0,:) = NaN;
    y = double(ismember(thermodynamic_data{n}.("Reaction type"),"TDR"));
    y = y(idx_true_reaction);

    y = y(~isnan(x(:,1)));
    x = x(~isnan(x(:,1)),:);
    cv_label = randi(5,length(y),1);
    ypred = zeros(size(y));
    for i = 1:5
        idx_test = find(cv_label == i);
        idx_train = find(cv_label ~= i);
        xtrain = x(idx_train,:);
        ytrain = y(idx_train);
        xtest = x(idx_test,:);
        rfmodel = fitrensemble(xtrain,ytrain);
        ypred(idx_test) = predict(rfmodel,xtest);
    end
    [xx,yy,~,auc] = perfcurve(y,ypred,1);
    %subplot(1,2,n);
    plot(xx,yy,'Color',line_colors(n,:));
    hold on;
    plot([0 1],[0 1],':','Color',[0 0 0]);
    xlabel('False positive rate');
    ylabel('True positive rate');
    annotation('textbox',[n*0.3 0.6 0.3 0.1],'String',...
        strcat(prediction_model_names{n},sprintf(",AUC=%.2f",auc)));   
end
box on;
title("ROC for the random-forest model predicting TDR/NTDR from network topology");

%% Plot distribution of network topology features in different types of rxns
% Do not plot distance to uptake fluxes in this step
figure;
m = n_topo_features - 1;
for n = 1:2
    for i = 1:m
        subplot(2,m,i+(n-1)*m);
        data = rxn_topological_features{:,i};
        data(isinf(data)) = NaN;

        data(is_nonzero_flux(:,n)==0) = NaN; %Discard rxns with zero fluxes

        labels = thermodynamic_data{n}.("Reaction type");
        data = data(idx_true_reaction);

        data = log(data + rxn_topo_features_offset(i)); %Log-transform the data

        labels = labels(idx_true_reaction);
        violinplot(data,labels,'ViolinColor',[87 127 219]/255);
        title(strcat(prediction_model_names{n},'.',topo_feature_names{i}));
        p = ranksum(data(ismember(labels,"TDR")),data(ismember(labels,"NTDR")));
        annotation('textbox',[i*0.2 (3-n)*0.4 0.1 0.1],'String', ...
            sprintf("Wilcoxon's rank-sum p = %.2e",p), 'EdgeColor','none');
        box on;
    end
end

%% Plot distribution of rxn deltaG vs distance to uptake flux
figure;
labels = rxn_topological_features.("Distance to nutrient uptake");
labels(labels>7) = 7;
labels(labels<2) = 2;
labels_full = labels;
labels_str = ["<=2","3","4","5","6",">=7"];
median_std_dG = zeros(2,6);
median_rxn_dG = zeros(2,6);

for n = 1:2
    subplot(1,2,n);
    data = thermodynamic_data{n}.("Average dGr");
    
    flux_table = flux_data{n};
    flux_matrix_tr = table2array(flux_table)';
    idx_nonzero = find(max(abs(flux_matrix_tr)) > 0);
    data = data(intersect(idx_nonzero, idx_true_reaction));
    labels = labels_full(intersect(idx_nonzero, idx_true_reaction));
    violinplot(data(~isnan(data)),labels(~isnan(data)),...
        'ViolinColor',[87 127 219]/255);
    box on;
    title(prediction_model_names{n});
    ylabel("Average \Delta_rG");
    xticklabels(labels_str);
    ylim([-500 300]);
    xlabel("Distance to nutrient uptake reaction");
    for i = 1:6
        x = data(labels == i+1);
        median_rxn_dG(n,i) = nanmedian(x);
    end
end

figure;
for n = 1:2
    subplot(1,2,n);
    data = thermodynamic_data{n}.("Standard dGr");
    
    data = data(idx_true_reaction);
    labels = labels_full(idx_true_reaction);
    violinplot(data(~isnan(data)),labels(~isnan(data)),...
        'ViolinColor',[87 127 219]/255);
    box on;
    title(prediction_model_names{n});
    ylabel("\Delta_rG^o");
    xticklabels(labels_str);
    ylim([-500 300]);
    xlabel("Distance to nutrient uptake reaction");
    for i = 1:6
        x = data(labels == i+1);
        median_std_dG(n,i) = nanmedian(x);
    end
end

%% Compare proteomics-based features between TDRs and NTDRs
figure;
for n = 1:2
    subplot(1,2,n);
    labels = thermodynamic_data{n}.("Reaction type");
    data = heterogeneity_score(idx_true_reaction);

    data(is_nonzero_flux(idx_true_reaction,n)==0) = NaN; %Discard rxns with zero fluxes

    labels = labels(idx_true_reaction);
    violinplot(data,labels,'ViolinColor',[87 127 219]/255);
    box on;
    title(prediction_model_names{n});
    ylabel("PCA-based heterogeneity score");
    p = ranksum(data(ismember(labels,"TDR")),data(ismember(labels,"NTDR")));
    annotation('textbox',[n*0.3 0.6 0.1 0.1],'String', ...
        sprintf("Wilcoxon's rank-sum p = %.2e",p), 'EdgeColor','none');
end

figure;
for n = 1:2
    subplot(1,2,n);
    labels = thermodynamic_data{n}.("Reaction type");
    data = average_protein_abundance(idx_true_reaction);

    data(is_nonzero_flux(idx_true_reaction,n)==0) = NaN; %Discard rxns with zero fluxes
    data(data == 0) = NaN; % Discard rxns with zero protein expression

    labels = labels(idx_true_reaction);
    violinplot(log10(data),labels,'ViolinColor',[87 127 219]/255);
    box on;
    title(prediction_model_names{n});
    ylabel("log10(protein abundance)");
    p = ranksum(data(ismember(labels,"TDR")),data(ismember(labels,"NTDR")));
    annotation('textbox',[n*0.3 0.6 0.1 0.1],'String', ...
        sprintf("Wilcoxon's rank-sum p = %.2e",p), 'EdgeColor','none');
end
