function wDatabaseLayers( counts, titles, abstracts, links, images, layers, curr_target )

[Z,T] = size( counts );

if isempty( titles ); for t=1:T; titles{t} = 'na'; end; end
if isempty( titles ); for t=1:T; abstracts{t} = 'na'; end; end
if isempty( titles ); for t=1:T; links{t} = 'http://www.bing.com'; end; end
if isempty( titles ); for t=1:T; images{t} = 'no_image.jpg'; end; end


fid = fopen([curr_target, '\database.txt'],'w+');
for t=1:T
    tmp = ['id:',num2str(t),'\t','title:',titles{t},'\t','abstract:',abstracts{t},'\t',...
        'link:',links{t},'\t','image:',images{t},'\t','layer:',num2str(layers(t)),'\n']; 
    fprintf( fid, tmp );
end
b = fclose( fid );

fid = fopen([curr_target,'\words.txt'],'w+');
for z=1:Z
    tmp = ['id:',num2str(z)];
    idw = find( counts(z,:) ~= 0);
    for w=1:length(idw)
        tmp = [tmp,'\t',num2str( idw(w)),':',num2str( counts(z,idw(w)))];
    end
    tmp = [tmp,'\n'];
    fprintf( fid, tmp );
end
b = fclose(fid);

fid = fopen([curr_target,'\cooccurrences.txt'],'w+');
for z=1:Z

    sign = 1;
    idw = find( counts(z,:) ~= 0);
    words = [];
    for i=1:length( idw )
        words = [words,find( counts(:,idw(i)) ~= 0 )'];
    end
    words = unique( words );
    
    nwords = [];
    for i=1:length( idw )
        nwords = [nwords,find( counts(:,idw(i)) ~= 0 )'];
    end
    nwords = setdiff(1:Z, unique( nwords ) );
    
    if length( nwords) < length( words)
        words = nwords;
        sign = -1;
    end
    tmp = ['id:',num2str(z),'\t','sig:',num2str(sign)];
    
    for w=1:length(words)
        tmp = [tmp,'\t',num2str( words(w))];
    end
    tmp = [tmp,'\n'];
    fprintf( fid, tmp );
end
b = fclose(fid);



