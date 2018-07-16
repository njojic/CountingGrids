function [qla,ql,pi_la,dirip,plal] = cg_layers( WD, pi, ql, W, L )

T = size(WD,2);
[cg_size(1),cg_size(2),Z]= size(pi);

pi_la = zeros([size(pi),L]);
h_la = zeros([size(pi),L]);
for l =1:L
    pi_la(:,:,:,l) = pi + rand(size(pi));
    pi_la(:,:,:,l) = bsxfun(@rdivide, pi_la(:,:,:,l) , sum(pi_la(:,:,:,l) ,3));
    PI = padarray( ...
        permute(cumsum( permute( cumsum( ...
        padarray(pi_la(:,:,:,l),W,'circular','post' )),[2 1 3]) ),[2 1 3] ), ...
        [1 1],0,'pre');
    tmp = compute_h_noLoopFull( PI, W(2)+1, W(1)+1);
    h_la(:,:,:,l) = bsxfun( @rdivide, tmp(1:end-1,1:end-1,:), sum( tmp(1:end-1,1:end-1,:),3 ));
end

alpha = 1e-10;
miter = 1;
P = prod( cg_size );
lqla = zeros([L,T]);
pseudocounts =  mean( sum(WD) / prod(cg_size) )  / 2.5;
plal= ones([L,P]) ./ L;
dirip = ones([L,Z]);

qlsm = zeros( size( ql));
mask = padarray( ones(W),[cg_size-W],0,'post');
qlsm = real( ifft2( bsxfun(@times, fft2( ql ),fft2( mask) ))) ./ prod(W);
every_iter = 1;
nmax = 50;
loglik = zeros(1,nmax);
start_ql = 2;
minp = 1e-10;
for iter = 1:nmax
    
    if iter >= start_ql
        
        if iter == start_ql;
            save('Model\before_iterating_ql','pi_la','qla');
%             qla = qla + rand( size( qla));
%             qla = bsxfun(@rdivide, qla, sum( qla ));
        end
        
        lql = zeros( [P,T] );
        for l=1:L
            tmp = reshape( log( eps + h_la(:,:,:,l)),[P,Z]);
            lql = lql + bsxfun( @times, tmp*WD, qla(l,:));
        end
        
        Lq = reshape( bsxfun( @minus, bsxfun(@minus, lql, max(lql) ),  log( sum( exp( bsxfun(@minus, lql, max(lql) ) ) ))), [cg_size,T]);
        ql = exp( Lq );
        qlsm = real( ifft2( bsxfun(@times, fft2( ql ),fft2( mask) ))) ./ prod(W);
        %tmp = exp( Lq ) ;  tmp( tmp< minp ) = minp;
        %ql = bsxfun(@rdivide, tmp, sum( sum(tmp)) );
    end
   
    
    tmpq = reshape( ql, [P,T]);
    for l=1:L
        tmp = reshape( log( eps + h_la(:,:,:,l)),[P,Z]);
        % lqla(l,:) = sum( tmpq.*(tmp*WD),1);
        lqla(l,:) = sum( tmpq.*( bsxfun(@plus, (tmp*WD), log( plal(l,:)') ) ),1);
    end
    qla = exp( bsxfun( @minus, bsxfun( @minus, lqla, max( lqla)), log( sum( exp(   bsxfun( @minus, lqla, max( lqla)) ))) ) );
    
    for l=1:L
        
        tmpdirip = dirip(l,:)-1;
        tmpdirip( isnan( tmpdirip ) ) = 0;
        
        nrm = bsxfun(@plus, reshape(tmpdirip,[1 1 Z]), ...
            reshape(  reshape( padarray( ql, W, 'circular','pre'), [prod(cg_size+W),T])*bsxfun( @times, WD, qla(l,:))', [ cg_size+W,Z ]));
        
        % nrm = reshape(  reshape( padarray( ql, W, 'circular','pre'), [prod(cg_size+W),T])*bsxfun( @times, WD, qla(l,:))', [ cg_size+W,Z ]);
        
        QH = permute( cumsum( permute( cumsum(...
            bsxfun( @rdivide, nrm,   padarray(h_la(:,:,:,l)+prod(W)*alpha ,[W,0],'circular','pre')) ), ...
            [2 1 3]) ),[2 1 3]);
        QH = compute_h_noLoopFull(QH,W(2)+1,W(1)+1);
        QH(QH<0) = 0;
        
        
        un_pi =   pseudocounts  + QH.*(pi_la(:,:,:,l)+alpha);
        mask = sum(un_pi,3) ~= 0;
        
        pi_la(:,:,:,l) = bsxfun(@times, bsxfun(@rdivide, un_pi, sum(un_pi,3)), double( mask ) ) ...
            + bsxfun(@times, 1/Z*ones([cg_size,Z]), double( ~mask ) );
        
%         % Smooth pi_la
%         tmppeak = reshape( pi_la(:,:,:,l), [P,Z]);
%         tmppeak = tmppeak.^(2);
%         tmppeak = bsxfun(@rdivide, tmppeak, sum( tmppeak,2));
%         pi_la(:,:,:,l) = reshape( tmppeak,[cg_size,Z]);
        

        PI = padarray( ...
            permute(cumsum( permute( cumsum( ...
            padarray(pi_la(:,:,:,l),W,'circular','post' )),[2 1 3]) ),[2 1 3] ), ...
            [1 1],0,'pre');
        tmp = compute_h_noLoopFull( PI, W(2)+1, W(1)+1);
        h_la(:,:,:,l) = bsxfun( @rdivide, tmp(1:end-1,1:end-1,:), sum( tmp(1:end-1,1:end-1,:),3 ));
        
%         if mod(iter, every_iter ) ==1
%             input = reshape( pi_la(:,:,:,l), [P,Z]);
%             % input = input.^(1.2);
%             input = bsxfun(@rdivide, input, sum( input,2));
%             dirip(l,:) = dirichlet_fit_simple(input);
%         end
        
    end

    plal =  bsxfun(@rdivide, reshape( squeeze( sum( bsxfun(@times, qlsm, reshape( qla',[1 1 T L]) ),3)), [P,L]), ...
        reshape( sum( qlsm,3),[1,P])')';
    plal( plal < 1e-100 ) = 1e-100;
    plal = bsxfun( @rdivide, plal, sum( plal,1));
    
%     oL = zeros(1,3);
%     for l=1:L
%         tmp = reshape( log(  h_la(:,:,:,l)),[P,Z]);
%         % lqla(l,:) = sum( tmpq.*(tmp*WD),1);
%         oL(l) = qla(l,:)*sum( tmpq.*( bsxfun(@plus, (tmp*WD), log( plal(l,:)') ) ),1)'; 
%     end
%     pL = sum( sum(sum( bsxfun( @times, reshape( bsxfun( @times, ql, reshape( qla',[1 1 T L])),[P,T,L]), reshape( log( plal)', [P,1,L])) )));
%     loglik(iter) = pL + sum( oL );
%     plot(1:iter,loglik(1:iter),'r'); drawnow;
%     
    
end
save('Model\after_iterating_ql','pi_la','qla');
end

function h = compute_h_noLoopFull( H, xW, yW )
        h = H(yW:end,xW:end,:,:,:) - H(1:end-yW+1,xW:end,:,:,:) ...
            - H(yW:end,1:end-xW+1,:,:,:) + H(1:end-yW+1,1:end-xW+1,:,:,:);    
end