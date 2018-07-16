function wMappingsLayers( q, W, idLayer, curr_target )

% The mapping of documents to CG locations is stored by this routine

[cg_size(1),cg_size(2),T]= size(q);
mask = padarray( ones(W),[cg_size-W],0,'post');
mappings = cell(T,1);
for t=1:T
    thr = 1e-2;
    tmp_q = q(:,:,t);
    conn_comp = bwconncomp( tmp_q>thr, 4);
    
    for r=1:conn_comp.NumObjects
        tmp = zeros( conn_comp.ImageSize );
        tmp( conn_comp.PixelIdxList{r})=1;
        id_conn_comp = find( tmp==1);
        [~,idmax] = max( tmp_q( id_conn_comp  ) );
        set_to_zero = setdiff( id_conn_comp , id_conn_comp(idmax) );
        tmp_q(set_to_zero) = 0;
    end
    tmp_q( tmp_q<=thr )=0;
    idm = find( tmp_q ~=0, 1);
    if isempty( idm )
        tmp_q2 = q(:,:,t);
        [~,idm] = max( tmp_q2(:) );
        tmp_q = zeros( cg_size );
        tmp_q( idm ) = tmp_q2(idm);
    end
    
    tmp_qs = real( ifft2( fft2( tmp_q).*fft2( mask )));
    tmp_qs = double( tmp_qs  > thr ).*tmp_qs;
    idm = find( tmp_qs~=0);
    if isempty( idm )
        [~,idm] = max( tmp_qs(:) );
    end
    for m=1:length( idm )
        [tmp_r,tmp_c] = ind2sub( cg_size, idm(m));
        mappings{t} = [mappings{t};[tmp_r,tmp_c,tmp_qs( tmp_r,tmp_c )]];
    end

end

posmap = cell( cg_size );
for y=1:cg_size(1)
    for x=1:cg_size(2)
        posmap{y,x} = ['row:',num2str(y),'\t','col:',num2str(x)];
    end
end

for t=1:T
    for i=1:size( mappings{t},1 )
        posmap{ mappings{t}(i,1), mappings{t}(i,2) } = [posmap{ mappings{t}(i,1),mappings{t}(i,2) },'\t',num2str(t),':',num2str( abs( mappings{t}(i,3) ) ),':',num2str(idLayer(t))];
    end
end

fid = fopen([curr_target,'\docmap.txt'],'w+');
for y=1:cg_size(1)
    for x=1:cg_size(2)
        fprintf( fid, [posmap{y,x},'\n']);
    end
end
b=fclose(fid);