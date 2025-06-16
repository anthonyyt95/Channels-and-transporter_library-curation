%% SBP: Processes list copied from SLC BioParadigms

data = copied;

% Acquires SLCs
loc = find(strcmp(data(:,1), 'SLC name'));
output = {};
for i = [1:length(loc(:,1))]
    start_ind = loc(i) + 1;
    
    for j = [start_ind:length(data(:,1))]
        val = data{j,1};
        
        if isempty(val)
            break
        end
        
        output = [output; data(j,:)];
    end
end

%% SBP: Removes 'pseudogenes'

data = list;

exclude = [];
for i = [1:length(data(:,1))]
    val = data{i,1};
    
    loc = strfind(lower(val), 'ensg000');
    if not(isempty(loc))
        exclude = [exclude; i];
    end
end
data(exclude,:) = [];
output = data;

clear data exclude i loc val



%% HGNC: Processes list copied from HGNC

data = list;

% Acquires SLCs
loc = find(strcmp(data(:,1), 'Gene'));
output = {};
step = 1;
for i = [1:length(loc(:,1))]
    ind = loc(i)-1;
    val = data{ind,1};
    
    loc1 = strfind(val,':');
    gene = val([1:loc1-1]);
    name = val([loc1+2:end]);
    
    output{step,1} = gene;
    output{step,2} = name;
    step = step + 1;
end

clear data loc ind val loc1 gene name


%% Removes duplicates

uniques = list;

[temp unique_ind] = unique(uniques(:,1),'stable');
comp_ind = [1:length(uniques(:,1))]';
exclude = [];
for j = [1:length(comp_ind(:,1))]
    cind = comp_ind(j);
 
    match = 0;
    for k = [1:length(unique_ind(:,1))]
        uind = unique_ind(k);
        if uind == cind
            match = 1;
            break
        end
    end
 
    if match == 0
        exclude = [exclude; j];
    end
end
uniques(exclude,:) = [];
uniques = sortrows(uniques,1,'asc');

clear temp i unique_ind comp_ind exclude cind j match uind

%% Looks for genes with an interest term in the name

data = list;
%query = ' channel'
%query = ' cation'
%query = ' anion'
%query = ' ion'
%query = 'transporter'
%query = 'sodium'
%query = 'na+'
%query = 'k+'
%query = 'h+'
%query = 'proton'
%query = 'potassium'
%query = 'calcium'
%query = 'ca2+'
%query = 'copper'
%query = 'magnesium'
%query = 'mg2+'
%query = 'chloride'
%query = 'cl-'
%query = 'uniporter'
%query = 'synport'
%query = 'antiport'
%query = 'exchanger'
%query = 'voltage'
%query = 'pore'
%query = 'hyperpolariz'
%query = 'depolariz'
%query = 'potential'
%query = 'carrier'

ind = [];
for i = [1:length(data(:,1))]
    val = lower(data{i,2});
    
    loc = strfind(val, lower(query));
    if not(isempty(loc))
        ind = [ind; i];
    end
end
contains = data(ind,:);
data(ind,:) = [];
does_not_contain = data;

clear data i ind loc val

%% Removes overlapping items (from list1)

data1 = list1;
data2 = list2;

data2 = lower(data2(:,1));
exclude = [];
for i = [1:length(data1(:,1))]
    val2 = lower(data1{i,1});
    
    loc = find(strcmp(val2, data2));
    if not(isempty(loc))
        exclude = [exclude; i];
    end
end
data1(exclude,:) = [];
output = data1;

clear data1 data2 exclude val2 loc i



%% Searches each gene from UniProt & acquires brief function info
tic
data = list;
threshold = 10;

output = {};
step = 1;
for i = [1:length(data)]
    gene_item = data{i,1};
    i
    url_parent = {'https://www.uniprot.org/uniprot/?query=';'&sort=score'};
    url = [string(url_parent{1}) + string(gene_item) + string(url_parent{2})];
    source = urlread(char(url));
    permanent = source;
    
    % Trims source file so it begins with the 1st entry
    try
        op_1 = strfind(source, 'class=" entry');
        source = source([op_1(1)+10:end]);
    catch
        continue
    end
    
    % Looks for the option corresponding to Homo Sapiens
    cont = 0;
    wstep = 1;
    try
        while cont == 0
            if not(wstep == 1)
                op_1 = strfind(source, 'class=" entry');
                source = source([op_1(1)+2:end]);
            end
            
            loc = strfind(source, 'href="/taxonomy');
            taxo = source([loc(1):loc(1)+1000]);
            loc1 = strfind(taxo,'>');
            loc2 = strfind(taxo,'<');
            taxo = taxo([loc1(1)+1:loc2(1)-1]);

            if not(lower(string(taxo)) == 'homo sapiens (human)')
                if wstep == threshold
                    continue
                else
                    wstep = wstep + 1;
                end
            else
                cont = 1;
            end
        end
    catch
        continue
    end
    
    % Acquires the 1st symbol shown in bold
    try
        loc1 = strfind(source, '<span class="shortName">');
        first_op = source([loc1(1)+24:end]);
        loc1 = strfind(first_op, '<strong>');
        loc2 = strfind(first_op, '</strong>');
        first_op = first_op([loc1(1)+8:loc2(1)-1]);
    catch
        continue
    end
    
    % Acquires the full gene name
    try
        loc1 = strfind(source, '<div class="protein_names">');
        gname = source([loc1(1):end]);
        loc1 = strfind(gname, 'title=');
        gname = gname([loc1(1):loc1(1)+1000]);
        loc1 = strfind(gname,'"');
        gname = gname([loc1(1)+1:loc1(2)-1]);
    catch
        continue
    end
   
    
    if lower(string(gene_item)) == lower(string(first_op))
        output{step,1} = gene_item;
        output{step,2} = gname;
        step = step + 1;
    else
        continue
    end 
end
toc

clear data first_op gene_item gname i loc1 loc2 permanent source step
clear url url_parent taxo

%% Looks for genes with an interest term in the name

data = list;
%query = 'zinc finger';
%query = 'transcription factor';
%query = 'actin'
%query = 'peptidase'
%query = 'deaminase'
%query = 'g protein-coupled receptor'
%query = 'apolipoprotein'
%query = 'open reading frame'
%query = 'complement'Z
%query = 'collagen'
%query = 'chemokine ligand'
%query = 'centromere protein'
%query = 'centrosomal protein'
%query = 'helicase'
%query = 'cytochrome p450'
%query = 'dynein'
%query = 'extracellular'
%query = 'translation initiation factor'
%query = 'coagulation factor'
%query = 'growth factor'
%query = 'forkhead box'
%query = 'histone'
%query = 'major histocompatibility complex'
%query = 'ribonucleoprotein'
%query = 'interleukin'
%query = 'kinesin'
%query = 'keratin'
%query = 'LDL receptor'
%query = 'microtubule'
%query = 'ribosomal protein'
%query = 'myosin'
%query = 'nuclear receptor'
%query = 'nucleoporin'
%query = 'olfactory receptor'
%query = 'phosphodiesterase'
%query = 'phospholipase'
%query = 'polymerase'
%query = 'proteasome'
%query = 'ras oncogene'
%query = 'rna binding'
%query = 'ring finger protein'
%query = 'tata-box'
%query = 'taste'
%query = 'tubulin'
%query = 'ubiquitin'

ind = [];
for i = [1:length(data(:,1))]
    val = lower(data{i,2});
    
    loc = strfind(val, lower(query));
    if not(isempty(loc))
        ind = [ind; i];
    end
end
contains = data(ind,:);
data(ind,:) = [];
does_not_contain = data;

clear data i ind loc val query


%%
[o2,m2] = getGeneInfo(list(1:500));
xlswrite('a1.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(501:1000));
xlswrite('a2.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(1001:1500));
xlswrite('a3.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(1501:2000));
xlswrite('a4.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(2001:2500));
xlswrite('a5.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(2501:3000));
xlswrite('a6.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(3001:3500));
xlswrite('a7.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(3501:4000));
xlswrite('a8.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(4001:4500));
xlswrite('a9.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(4501:5000));
xlswrite('a10.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(5001:5500));
xlswrite('a11.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(5501:6000));
xlswrite('a12.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(6001:6500));
xlswrite('a13.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(6501:7000));
xlswrite('a14.xlsx',string(o2));

[o2,m2] = getGeneInfo(list(7001:7199));
xlswrite('a15.xlsx',string(o2));


%% Retrieves data from 'getGeneInfo'

a1 = importdata('a1.xlsx');
a2 = importdata('a2.xlsx');
a3 = importdata('a3.xlsx');
a4 = importdata('a4.xlsx');
a5 = importdata('a5.xlsx');
a6 = importdata('a6.xlsx');
a7 = importdata('a7.xlsx');
a8 = importdata('a8.xlsx');
a9 = importdata('a9.xlsx');
a10 = importdata('a10.xlsx');
a11 = importdata('a11.xlsx');
a12 = importdata('a12.xlsx');
a13 = importdata('a13.xlsx');
a14 = importdata('a14.xlsx');
a15 = importdata('a15.xlsx');

output = [a1;a2;a3;a4;a5;a6;a7;a8;a9;a10;a11;a12;a13;a14;a15];

%% Looks for genes with an interest term in the annotated description

data = list;
%query = 'channel'
%query = 'transporter'
%query = 'uniporter'
%query = 'symport'
%query = 'antiport'
%query = 'exchanger'
%query = 'voltage'
%query = 'pore'
%query = 'hyperpolariz'
%query = 'depolariz'
%query = 'sodium'
%query = 'na+'
%query = 'k+'
%query = 'potassium'
%query = 'h+'
%query = 'proton'
%query = 'chloride'
%query = 'cl-'

%query = ' cation'
%query = ' anion'
%query = ' ion'
%query = 'carrier'

ind = [];
for i = [1:length(data(:,1))]
    val = lower(data{i,2});
    
    loc = strfind(val, lower(query));
    if not(isempty(loc))
        ind = [ind; i];
    end
end
contains = data(ind,:);
data(ind,:) = [];
does_not_contain = data;

clear data i ind loc val query

%% Looks for genes with an interest term in the annotated GO description

data = list;
%query = 'ion transmembrane transport';
query = 'ion transport'

ind = [];
for i = [1:length(data(:,1))]
    val_molec = lower(data{i,3});
    val_bp = lower(data{i,4});
    
    loc = strfind(val_molec, lower(query));
    if not(isempty(loc))
        ind = [ind; i];
        continue
    end
    
    loc = strfind(val_bp, lower(query));
    if not(isempty(loc))
        ind = [ind; i];
        continue
    end
end
contains = data(ind,:);
data(ind,:) = [];
does_not_contain = data;

clear data i ind loc val query val_bp val_molec


%% Acquires data given a list of genes and a reference

lref = ref;
data = list;

step = 1;
for i = [1:length(data(:,1))]
    gene = data{i};
    
    loc = find(strcmp(lref(:,1),gene));
    output(step,:) = lref(loc,:);
    step = step + 1;
end

clear lref data step i gene loc 


%% Processes ICTs derived from HUGO, UniProt, & other online databases
%
%
%
%

%% Searches each gene from UniProt & acquires brief function info
tic
data = list;

output = {};
step = 1;
for i = [1:length(data)]
    gene_item = data{i,1};
    i
    url_parent = {'https://www.uniprot.org/uniprot/?query=';'&sort=score'};
    url = [string(url_parent{1}) + string(gene_item) + string(url_parent{2})];
    source = urlread(char(url));
    permanent = source;
    
    try
        loc1 = strfind(source, '<span class="shortName">');
        first_op = source([loc1(1)+24:end]);
        loc1 = strfind(first_op, '<strong>');
        loc2 = strfind(first_op, '</strong>');
        first_op = first_op([loc1(1)+8:loc2(1)-1]);
    catch
        continue
    end
        
    if lower(string(gene_item)) == lower(string(first_op))
        item_find = '<tr id=';
        loc1 = strfind(source,item_find);
        source = source([loc1:end]);

        item_find = '"';
        loc1 = strfind(source, item_find);
        gene_id = source([loc1(1)+1:loc1(2)-1]);

        url = ['https://www.uniprot.org/uniprot/' + string(gene_id)];
        source = urlread(char(url));

        loc1 = strfind(source, 'This section provides any useful information');
        try
            loc2 = strfind(source, 'GO - Molecular function');
            source = source([loc1:loc2(1)]);
        catch
            try
                loc2 = strfind(source, 'GO - Biological process');
                source = source([loc1:loc2(1)]);
            catch
                loc2 = strfind(source,'<h4>');
                source = source([loc1:loc2(1)]);
            end
        end
        
        try
            loc1 = strfind(source, '<span property="text">');
            loc1 = loc1(1) + 22;
            source = source([loc1:end]);
        catch
            output{step,1} = gene_item;
            output{step,2} = '';
            step = step + 1;
            continue
        end
        
        try
            loc2 = strfind(source, '<span class="attribution ECO269">');
            source = source([1:loc2(1)-1]);
        catch
            try
                loc2 = strfind(source, '</div>');
                source = source([1:loc2(1)-1]);
            catch 
                output{step,1} = gene_item;
                output{step,2} = "";
                step = step + 1;
                continue
            end
        end
       
        output{step,1} = gene_item;
        output{step,2} = string(extractHTMLText(source));
        step = step + 1;
    else
        continue
    end 
end
toc




%% Removes items enclosed by <> in a string

data = list;

for i = [1:length(data(:,1))]
    descrip = data{i,1};
    if isempty(descrip)
        continue
    end
    
    locL = strfind(descrip, '<');
    locR = strfind(descrip, '>');
    if isempty(locL)
        continue
    end
    
    exclude = [];
    for j = [1:length(locL)]
        try
            exclude = [exclude, [locL(j):locR(j)]];
        catch
            break
        end
    end
    descrip = char(descrip);
    descrip(exclude) = [];
    data{i,1} = regexprep(descrip,'\n+','');
    
end
    
%% Searches for keywords in the ICT gene function descriptions
%
% Input
%   'data'          X-by-2 matrix consisting of gene names (col 1) & their
%                   associated description (col 2)
%   'query'         character variable for term of interest
%
% Output
%   'output'        X-by-2 matrix containing gene names (col 1) & their
%                   associated description (col 2)
%

data = list;
query = 'ion';
col = 2;

output = {};
step = 1;
for i = [1:length(data(:,1))]
    descrip = data{i,col};
    loc = strfind(lower(string(descrip)), lower(string(query)));
    if isempty(loc)
        continue
    else
        output{step,1} = data{i,1};
        output{step,2} = descrip;
        step = step + 1;
    end
end


%% Excludes genes in list A from list B
%
% Input
%   'data'          X-by-1 list from which genes will be excluded
%   'list'          X-by-1 list containing the genes to be excluded from
%                   'data'
%
% Output
%   'output'        X-by-1 list of genes from 'data' with the genes from
%                   'list' excluded

data = ticts;
list = aicts;

exclude = [];
for i = [1:length(data(:,1))]
    data_val = data{i,1};
    for j = [1:length(list(:,1))]
        list_val = list{j,1};
        if lower(string(data_val)) == lower(string(list_val))
            exclude = [exclude; i];
        end
    end
end

data(exclude,:) = [];
output = data;

%% Searches each gene from UniProt & acquires GO info
tic
data = list;

output = {};
step = 1;
for i = [1:length(data)]
    gene_item = data{i,1};
    i
    url_parent = {'https://www.uniprot.org/uniprot/?query=';'&sort=score'};
    url = [string(url_parent{1}) + string(gene_item) + string(url_parent{2})];
    source = urlread(char(url));
    permanent = source;
    
    try
        loc1 = strfind(source, '<span class="shortName">');
        first_op = source([loc1(1)+24:end]);
        loc1 = strfind(first_op, '<strong>');
        loc2 = strfind(first_op, '</strong>');
        first_op = first_op([loc1(1)+8:loc2(1)-1]);
    catch
        continue
    end
        
    if lower(string(gene_item)) == lower(string(first_op))
        item_find = '<tr id=';
        loc1 = strfind(source,item_find);
        source = source([loc1:end]);

        item_find = '"';
        loc1 = strfind(source, item_find);
        gene_id = source([loc1(1)+1:loc1(2)-1]);

        url = ['https://www.uniprot.org/uniprot/' + string(gene_id)];
        source = urlread(char(url));
        permanent = source;
        
        try
            loc1 = strfind(permanent, 'GO - Molecular function');
            source = permanent([loc1(1):end]);
            
            loc1 = strfind(source, '<ul');
            source = source([loc1(1):end]);
            loc1 = strfind(source, '<li>');
            loc2 = strfind(source, '<div');
            source1 = source([loc1(1):loc2(1)]);
        catch
        end
    

        try
            loc1 = strfind(permanent, 'GO - Biological process');
            source = permanent([loc1(1):end]);
            
            loc1 = strfind(source, '<ul');
            source = source([loc1(1):end]);
            loc1 = strfind(source, '<li>');
            loc2 = strfind(source, '<div');
            source2 = source([loc1(1):loc2(1)]);
        catch
        end
        
        output{step,1} = gene_item;
        output{step,2} = string(source1);
        output{step,3} = string(source2);
        step = step + 1;
    else
        continue
    end 
end
toc


%% Searches each gene from UniProt & acquires GO info
tic
data = list;

output = {};
step = 1;
for i = [1:length(data)]
    gene_item = data{i,1};
    i
    url_parent = {'https://www.uniprot.org/uniprot/?query=';'&sort=score'};
    url = [string(url_parent{1}) + string(gene_item) + string(url_parent{2})];
    source = urlread(char(url));
    permanent = source;
    
    try
        loc1 = strfind(source, '<span class="shortName">');
        first_op = source([loc1(1)+24:end]);
        loc1 = strfind(first_op, '<strong>');
        loc2 = strfind(first_op, '</strong>');
        first_op = first_op([loc1(1)+8:loc2(1)-1]);
    catch
        continue
    end
        
    if lower(string(gene_item)) == lower(string(first_op))
        item_find = '<tr id=';
        loc1 = strfind(source,item_find);
        source = source([loc1:end]);

        item_find = '"';
        loc1 = strfind(source, item_find);
        gene_id = source([loc1(1)+1:loc1(2)-1]);

        url = ['https://www.uniprot.org/uniprot/' + string(gene_id)];
        source = urlread(char(url));
        permanent = source;
        
        try
            loc1 = strfind(permanent, 'GO - Molecular function');
            source = permanent([loc1(1):end]);
            
            loc1 = strfind(source, '<ul');
            source = source([loc1(1):end]);
            loc1 = strfind(source, '<li>');
            loc2 = strfind(source, '<div');
            source1 = source([loc1(1):loc2(1)]);
        catch
        end
    

        try
            loc1 = strfind(permanent, 'GO - Biological process');
            source = permanent([loc1(1):end]);
            
            loc1 = strfind(source, '<ul');
            source = source([loc1(1):end]);
            loc1 = strfind(source, '<li>');
            loc2 = strfind(source, '<div');
            source2 = source([loc1(1):loc2(1)]);
        catch
        end
        
        output{step,1} = gene_item;
        output{step,2} = string(source1);
        output{step,3} = string(source2);
        step = step + 1;
    else
        continue
    end 
end
toc

%% Acquires all aliases of a list of genes (based on UniProt)

data = list;
threshold = 20;
output = {};
step = 1;
missing = {};

for i = [1:length(data)]
    gene_item = data{i,1};
    i
    url_parent = {'https://www.uniprot.org/uniprot/?query=';'&sort=score'};
    url = [string(url_parent{1}) + string(gene_item) + string(url_parent{2})];
    source = urlread(char(url));
    
    % Trims source file so it begins with the 1st entry
    try
        op_1 = strfind(source, 'class=" entry');
        source = source([op_1(1)+10:end]);
    catch
        missing = [missing; gene_item];
        continue
    end
    
    % Looks for the option corresponding to Homo Sapiens
    cont = 0;
    wstep = 1;
    try
        while cont == 0
            if not(wstep == 1)
                op_1 = strfind(source, 'class=" entry');
                source = source([op_1(1)+2:end]);
            end
            
            loc = strfind(source, 'href="/taxonomy');
            taxo = source([loc(1):loc(1)+1000]);
            loc1 = strfind(taxo,'>');
            loc2 = strfind(taxo,'<');
            taxo = taxo([loc1(1)+1:loc2(1)-1]);

%             if not(lower(string(taxo)) == 'homo sapiens (human)')
            if not(lower(string(taxo)) == 'mus musculus (mouse)')    
                if wstep == threshold
                    continue
                else
                    wstep = wstep + 1;
                end
            else
                % Checks to determine if gene symbol is present
                try
                    loc1 = strfind(source, '<span class="shortName">');
                    loc2 = strfind(source, '</span>');
                    alias = extractHTMLText(source([loc1(1):loc2(1)+6]));
                    alias = strtrim(alias);
                    alias = strrep(alias,',','');
                    cont = 1;
                catch
                    continue
                end
            end
        end
    catch
        missing = [missing; gene_item];
        continue
    end
    
    output{step,1} = alias;
    step = step + 1;
end

clear alias cont data gene_item i loc loc1 loc2 op_1 source step taxo
clear threshold url url_parent wstep

%% Cleans up situations in which an alias of a gene corresponds to another gene

list = ICT794;

output = list;
for i = [1:length(list(:,1))]
    aliases = strsplit(list{i,1}, ' ');
    gene = aliases{1};
    for j = [1:length(list(:,1))]
        if j == i
            continue
        else
            daliases = strsplit(list{j,1}, ' ');
            loc = find(strcmp(lower(daliases), lower(gene)));
            if not(isempty(loc))
                daliases(loc) = [];
                daliases = strjoin(daliases, ' ');
                output{j,1} = daliases;
            end
        end
    end
end

%% Acquires all aliases of a list of genes (based on HGNC)

data = list;
missing = {};

tic
% Acquires HGNC ID data ('2019.12.17_HGNCgeneIDs.txt')
filename = 'C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\2019.12.17_HGNCgeneIDs.xlsx';
ref = importdata(filename);
ref = ref.textdata;

% Acquires gene name synonyms
output = {};
for i = [1:length(data(:,1))]
    gene = data{i,1};
    loc = find(strcmp(lower(ref(:,2)), lower(gene)));
    if isempty(loc)
        missing = [missing; gene];
        output{i,1} = gene;
        output{i,2} = '';
    elseif length(loc) > 1
        for j = [1:length(loc)]
            if ref{loc(j),4} == 'Approved'
                loc = loc(j);
                output{i,1} = gene;
                output{i,2} = ref{loc,6};
                break
            end
        end
    else
        output{i,1} = gene;
        output{i,2} = ref{loc,6};
    end
end
toc

clear data filename ref i gene loc a




