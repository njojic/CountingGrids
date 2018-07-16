function wCountingGridLayers( pi_la, N, curr_target )


[cg_size(1), cg_size(2), Z, L] = size(pi_la);
top_pi = cell([L,cg_size]);

for lay=1:L
    
    pi = pi_la(:,:,:,lay);
    
    for z=1:Z
        tmp_pi = pi(:,:,z);
        conn_comp = bwconncomp( tmp_pi>1e-3, 8);
        
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
    

    TW = [];
    for y=1:cg_size(1)
        for x=1:cg_size(2)
            wordDist = squeeze( pi(y,x,:) );
            [~,idw] = sort( wordDist ,'descend');
            top_pi{lay,y,x} = idw(1:N);
            TW = [TW,idw(1:N)];
        end
    end
    
    missing_words = setdiff(1:Z, unique( TW) );
    for m=1:length( missing_words )
        tmp_pi = pi(:,:,missing_words(m));
        idw = find( tmp_pi > 1e-3 );
        if isempty( idw )
            [~,idw] = max( tmp_pi(:) );
        end
        
        for i=1:length( idw )
            [tmp_r,tmp_c] = ind2sub( cg_size, idw(i));
            top_pi{lay,tmp_r,tmp_c} = [top_pi{lay,tmp_r,tmp_c}; missing_words(m)];
        end
    end
end


    
fid = fopen([curr_target,'\top_pi_layers.txt'],'w');
for l=1:L
    for y=1:cg_size(1)
        for x=1:cg_size(2)
            tmp = ['layer:',num2str(l),'\t','row:',num2str(y),'\t','col:',num2str(x)];
            for n=1:length( top_pi{l,y,x} )
                tmp = [tmp,'\t',num2str(top_pi{l,y,x}(n)-1),':',num2str(pi_la(y,x,top_pi{l,y,x}(n),l) )];
                % tmp = [tmp,'\t',num2str(idw(n)),':',num2str(pi(y,x,idw(n)))];
            end
            tmp = [tmp,'\n'];
            fprintf( fid, tmp );
        end
    end
end
b = fclose( fid );

