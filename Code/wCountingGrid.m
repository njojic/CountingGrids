function wCountingGrid( pi, N, curr_target )

% To reduce the size of the CG which has probabilities for the entire vocab
% in each cell, we keep only the top words in each cell.
% Some vocab words may no be in the top for any cell. We add those to the
% cells where they had the max probability. THis way the entire voacabualry
% is represented in the stored CG.

[cg_size(1), cg_size(2), Z] = size(pi);

for z=1:Z
    tmp_pi = pi(:,:,z);
    conn_comp = bwconncomp( tmp_pi>1e-3, 4);
    
    tmp_sub = zeros( cg_size );
    for r=1:conn_comp.NumObjects
        tmp = zeros( conn_comp.ImageSize );
        tmp( conn_comp.PixelIdxList{r})=1;
        id_conn_comp = find( tmp==1);
        [~,idmax] = max( tmp_pi( id_conn_comp  ) );
        set_to_zero = setdiff( id_conn_comp , id_conn_comp(idmax) );
        tmp_pi(set_to_zero) = 0;
    end
    pi(:,:,z) = tmp_pi;
end

top_pi = cell(cg_size);
top_pi_out = zeros([cg_size,N]);
TW = [];
for y=1:cg_size(1)
    for x=1:cg_size(2)
        wordDist = squeeze( pi(y,x,:) );
        [~,idw] = sort( wordDist ,'descend');    
        top_pi{y,x} = idw(1:N);
        top_pi_out(y,x,:) = idw(1:N);
        TW = [TW,idw(1:N)];
    end
end
basedir=pwd;
cd(curr_target);
save top_pi_out top_pi_out
cd(basedir);

missing_words = setdiff(1:Z, unique( TW) );
for m=1:length( missing_words )
    tmp_pi = pi(:,:,missing_words(m));
    idw = find( tmp_pi > 1e-3 );
    if isempty( idw )
        [~,idw] = max( tmp_pi(:) );
    end
    
    for i=1:length( idw )
        [tmp_r,tmp_c] = ind2sub( cg_size, idw(i));
        top_pi{tmp_r,tmp_c} = [top_pi{tmp_r,tmp_c}; missing_words(m)];
    end
end


fid = fopen([curr_target ,'\top_pi.txt'],'w');
for y=1:cg_size(1)
    for x=1:cg_size(2)
        tmp = ['row:',num2str(y),'\t','col:',num2str(x)];
        for n=1:length( top_pi{y,x} )
            tmp = [tmp,'\t',num2str(top_pi{y,x}(n)-1),':',num2str(pi(y,x,top_pi{y,x}(n)))]; 
           % tmp = [tmp,'\t',num2str(idw(n)),':',num2str(pi(y,x,idw(n)))]; 
        end
        tmp = [tmp,'\n'];
        fprintf( fid, tmp );        
    end
end
b = fclose( fid );
