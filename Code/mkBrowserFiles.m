
% This script runs the CG learning then make the files for the CG browser
%
% titles : cell array with titles of the documents 
% abstracts : cell array with abstracts of the documents 
% links : cell array with the links for documents
% WD : Word document matrix
% VOC : Vocabulary (cell array) that refers to WD
% pi: the counting grid consisting of word probabilities in a grid


browser_data_dir = ['Browser\Data\',output_file_name]; 

%wd_size = [5 5]; % window size
%cg_size = [32,32]; % grid size

mkdir(browser_data_dir );

for t=1:size(counts,2)
    images{t}='http://www.bing.com/';
    % this is for a web version that serves an image for each abstract
    % Note enabled here
end

% Merge words with the same root; Can lead to suprising results, like using
% the word visualize for all of these words: vision, visual, visualize
[WD,VOC] = wDictionary_lessRestrictive(counts, VOC,browser_data_dir );


%%%%%%%%%%%%%%% COUNTING GRID LEARNING
% Coarse-to-fine approach here may be an overakill; results are often the
% same if CG routine is called only once. Here, the window is initially
% large, than reduced to the desired size wd_size. 100 iterations is also
% likely an overkill most of the time
clear options;
options.max_iter=15; % coarse to fine training: windows go form bigger to smaller=
[pi,pl,Lq,loglikelihood_samples] = cg( WD, cg_size, [10 10], options);

options.pi=pi;
[pi,pl,Lq,loglikelihood_samples] = cg( WD, cg_size, [6 6], options);

options.max_iter=100;
options.pi=pi;
[pi,pl,Lq,loglikelihood_samples] = cg( WD, cg_size, wd_size, options);



[Z,T] = size(WD);
% compute idf weights
idf = log( size(WD,2) ./ sum( WD ~= 0,2 )');
idf(isinf(idf))=0;

%%%%%%%%%%%%% LAYER ESTIMATION
% Layers speed up the search in the browser and also allow visualziation of
% the breaks in the topic drifts in the grid; See the KDD paper
L = 4; % number of layers
ql = exp(Lq);
qln = ql + 0.25*rand(size( ql)); % add noise
qln = bsxfun( @rdivide, qln, sum( qln,3));

[qla,ql2,pi_la,dirip,plal] = cg_layers( WD, pi, qln, wd_size, L );
[~,id_layer] = max( qla);

pi_idf = bsxfun(@times, pi, reshape( idf, [1 1 Z]));
pi_la_idf = zeros( size( pi_la));
for l=1:L
    idf_layer = log( size(WD(:,id_layer==l),2) ./ sum( WD(:,id_layer==l) ~= 0,2 )');
    idf_layer( isinf(idf_layer)) = 0;   
    %idf_layer=ones(size(idf_layer));
    pi_la_idf(:,:,:,l) = bsxfun(@times, pi_la(:,:,:,l), reshape( idf_layer, [1 1 Z]));
end


mask = padarray( ones( wd_size),cg_size-wd_size,0,'post');
wg = zeros( [cg_size,L]);
for l=1:L
    idl = find( id_layer == l);
    for t=1:length( idl )
        wg(:,:,l) =   wg(:,:,l)+ real( ifft2( fft2( mask ).*fft2( ql2(:,:, idl(t)))));
    end
end
wg = bsxfun( @rdivide, wg, (eps+ sum( wg,3)));

wg_lay=wg;

pi2_idf = squeeze( sum( bsxfun(@times, pi_la_idf, reshape( wg, [cg_size,1,L])),4));

% Clean the titles and abstracts from special characters if any
titlesC = clean_text_single( titles );
abstractsC = clean_text_single( abstracts );

% Make browser files
% These are all text files, and are thus big
% This limits the browser to datasets of up to 30-40k documents
% Larger datasets will be slow to write down on disk, or load into memory
% (The browser loads everything into memory)
wCountingGrid( pi2_idf, 80, browser_data_dir );
wDatabaseLayers( WD, titlesC, abstractsC, links, images, id_layer, browser_data_dir );
wMappingsLayers( ql2+eps , wd_size, id_layer, browser_data_dir );
wCountingGridLayers( pi_la_idf, 80, browser_data_dir );


% Compute and save a different color for each cell.
% If you do not want colors, just dont make this file (or erase it)
layer_colors=[1 1 0; 0 0.6 1; 1 0 0.6; 0.3 1 0.4    ]';

fid = fopen([browser_data_dir '\colors_browser.txt'],'wb');
for i = 1:cg_size(1)
    for j=1:cg_size(2)
        w=squeeze(wg_lay(i,j,:));
        w=w.^1.3; w=w/sum(w);
        color=layer_colors*w;
        line = [num2str(i),'\t',num2str(j),'\t',num2str( color(1)),'\t',num2str( color(2)),'\t',num2str( color(3)),'\n'];
        fprintf(fid,line);
    end
end

fclose(fid);

% Instead of coloring layers, we could color the CG cells in any other way.
% E.g. if we had labels for documents, indicating three document categories
% we can do something like this:
% colors_code = [1,0.2,0.1; 0.1,1,0.4; 0.3 0.2 1]; % palette for 3 labels
%     C = 3; % for R G B
%     mask = padarray( ones( wd_size),cg_size-wd_size,0,'post');
%     colors = zeros( [cg_size,C]);
%     den = zeros(cg_size);
%     
%     for t=1:length(abstracts)
%         wg = sum( labels == labels(t)) / length( labels );
%         colors(:,:, labels(t)) = colors(:,:,labels(t)) + real( ifft2( fft2( ql2(:,:,t)).*fft2(mask)))/wg;
%         den = den + real( ifft2( fft2( ql2(:,:,t)).*fft2(mask)))/wg;
%     end
%     colors = bsxfun(@rdivide, colors, den );
%     colors = colors.^1.5; % Sometimes I want peakier colors
%     colors = bsxfun(@rdivide, colors, sum(colors,3) );
%     colors = squeeze( sum( bsxfun(@times, colors, reshape( colors_code,[1,1,C,3] )),3) );
%     fid = fopen([browser_data_dir '\colors_browser.txt'],'wb');
%     for i = 1:cg_size(1)
%         for j=1:cg_size(2)
%             line = [num2str(i),'\t',num2str(j),'\t',num2str( colors(i,j,1)),'\t',num2str( colors(i,j,2)),'\t',num2str( colors(i,j,3)),'\n'];
%             fprintf(fid,line);
%         end
%     end
%     fclose(fid);

save Model\pi pi
save Model\pi_la pi_la

