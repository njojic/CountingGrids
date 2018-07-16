function [counts_output,vocabulary_display] = wDictionary_lessRestrictive( counts, vocabulary, targetpath )

% Merges words with the same stem, removes special characters, etc.
% Also has one reamining comment from Alessandro that did not have bad
% language

[Z,T] = size(counts);

% Step 1. Pulisco il dizionario
vocabulary_full_cleaned = cell(1,Z);
for z=1:Z
   tmp = vocabulary{z};
   vocabulary_full_cleaned{z} = tmp( regexp( tmp, '[A-Za-z\d\-_!()*{}]'));
end
id = cellfun( @length, vocabulary_full_cleaned ) > 2;
counts_1 = counts(id,:);
vocabulary_1 = lower( {vocabulary_full_cleaned{id}} );
[Z,T] = size(counts_1);

stems = cell(Z,1);
for z=1:Z
   stems{z} = porterStemmer( vocabulary_1{z} );
end

unique_stems = unique( stems );
Zs = length( unique_stems );

l2s = zeros(1,Z);
for z=1:Z
   l2s(z) = find( strcmp( unique_stems, stems{z} )  == 1);
end

T = size(counts,2);
counts_stemmed = zeros( Zs,T );
for z=1:Z
    tmp = counts_1(z,:);
    counts_stemmed(l2s(z),:) = counts_stemmed(l2s(z),:)+tmp;
end

% Scegli 
vocabulary_display = cell(1,Zs);
for z=1:Zs
    id = find( l2s == z );
    candidate = {vocabulary_1{id}};
    len = cellfun(@length, candidate);
    % [~,id_sorted] = sort(len,'ascend');
    [~,id_sorted] = sort(len,'descend');
    candidate_chosen = [];
    
    for l=1:length( id_sorted )
        
        tmp = candidate{id_sorted(l)};
        try
            if ~strcmp( tmp(end-1:end),'ed') && ~strcmp( tmp(end-2:end), 'ing' ) && ~strcmp( tmp(end), 's' ) && ~strcmp( tmp(end-1:end), 'ty' )  && ~strcmp( tmp(end-1:end), 'ly' ) && ~strcmp( tmp(end-2:end), 'ism' )
                candidate_chosen = candidate{id_sorted(l)};
                break
            end
        catch
            candidate_chosen = tmp;
        end

    end
    
    if isempty( candidate_chosen  )
        
        for l=1:length( id_sorted )
            tmp = candidate{id_sorted(l)};
            if ~strcmp( tmp(end-1:end),'ed') && ~strcmp( tmp(end-2:end), 'ing' ) && ~strcmp( tmp(end), 's' )
                candidate_chosen = candidate{id_sorted(l)};
                break
            end
        end
        if isempty( candidate_chosen  )
            candidate_chosen  = candidate{id_sorted(1)};
        end
    end
    
    if strcmp( candidate_chosen,'powerful')
        candidate_chosen = 'power';
    end
    if ~isempty(  findstr(candidate_chosen, 'address') )
        candidate_chosen = 'address';
    end
    if strcmp( candidate_chosen,'binged')
        candidate_chosen = 'bing';
    end
        if ~isempty(  findstr(candidate_chosen, 'download') )
        candidate_chosen = 'download';
    end
       
    vocabulary_display{z} = candidate_chosen;
end

fid = fopen([targetpath,'\vocabulary.txt'],'w+');
for z=1:Zs
    tmp = [num2str(z),'\t',vocabulary_display{z},'\n'];
    fprintf( fid, tmp );
end
b = fclose( fid );


fid = fopen([targetpath,'\correspondences.txt'],'w+');
for z=1:Z
   tmp =  [vocabulary_1{z},'\t',vocabulary_display{l2s(z)},'\t',num2str(l2s(z)),'\n'];
    fprintf( fid, tmp );   
end
b = fclose( fid );

counts_output = counts_stemmed;





