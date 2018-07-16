function [FO,WO] = show_text_in_grid(  pi, vocabulary, options )


if ~exist( 'options', 'var'); options = []; end
if ~isfield( options,'minf'); options.minf = 8; end
if ~isfield( options,'maxf'); options.maxf =16; end
if ~isfield( options,'freq'); options.freq = ones(1,length( vocabulary )); end
if ~isfield( options,'no_words'); options.no_words = 3; end
if isfield( options,'fig'); 
    try 
        close(options.fig); 
    end
    handle=figure(options.fig); 
    %options.handle = handle; 
end
if ~isfield( options,'fig'); 
    handle=figure(1); 
    %options.handle = handle; 
end

load('font_to_pixels.mat');
%handle = options.handle;

[E(1),E(2),Z]=size(pi);
axis on
%% Ingrandisci la figura - enlarge figure
% TODO: Come metterla sul secondo schermo? - How to put it on the second screen?
set(handle, 'Position', get(0,'Screensize'));
set( get(handle,'Children'),'Position',get( get(handle,'Children'),'OuterPosition' ) );
       
axis off

%% Plotta le linee - Plot the lines
pos = get( handle, 'Position' );

lx = pos(3);
ly = pos(4);
axis([0 lx 0 ly]);
limx = linspace(0,lx,E(2)+1); limx=limx(2:end-1);
limy = linspace(0,ly,E(1)+1); limy=limy(2:end-1);
axis off
for l=1:length( limx )
    line( ones(1,ly).*limx(l),1:ly,'Color',[.7 .7 .7] );
end
for l=1:length( limy )
    line( 1:lx,ones(1,lx).*limy(l),'Color',[.7 .7 .7] );
end

%% Calcola midpoint and size di ogni cella (serve per il testo dopo) - Calculates the midpoint and size of each cell (serves for the text afterwards)
limx = linspace(0,lx,E(2)+1);
limy = linspace(0,ly,E(1)+1);

size_cell = zeros([E,2]);
midpoint_cell = zeros([E,2]);
for y=1:E(1)
    for x=1:E(2)
        size_cell(y,x,:) = [ limy(y+1)-limy(y), limx(x+1)-limx(x)];
        midpoint_cell(y,x,:) = [limy(y),limx(x)] + squeeze(size_cell(y,x,:))'./2;
    end
end


different_fonts = 10;
font_sizes = round( linspace( options.minf,options.maxf,different_fonts) );

prob = zeros([E,options.no_words]);
words_in_grid = zeros([E,options.no_words]);
used = [];
for y=1:E(1)
    for x=1:E(2)
        [values_sorted,indeces_sorted] = sort( squeeze( pi(y,x,:) )'.*options.freq ,'descend');
        words_to_show = {vocabulary{indeces_sorted(1:options.no_words)}};
        prob(y,x,:) = values_sorted(1:options.no_words);
        for w=1:options.no_words
            words_in_grid(y,x,w) = indeces_sorted(w);
        end
        used = cat(1,used,{words_to_show{:}}');
    end
end
used = unique( used );
U = length(used);

% BW = zeros([E,U]);
% for u=1:U
%     id = find( strcmp( vocabulary, used{u} )==1);
%     BW(:,:,u) = imregionalmax(pi(:,:,id),8);
% end

mask_word = zeros([E,options.no_words,U]);

for u=1:U
    id = find( strcmp( vocabulary, used{u} )==1);
    id = id(1);
    for w=1:options.no_words
        bw = double(words_in_grid(:,:,w)==id);
        % Trova il numero di regioni - Find the number of regions
        conn_comp = bwconncomp(bw,8);
        
        for r=1:conn_comp.NumObjects
            tmp = zeros( conn_comp.ImageSize );
            tmp( conn_comp.PixelIdxList{r})=1;
            stat = regionprops(tmp,'Centroid');
            xc = stat.Centroid(1);
            yc = stat.Centroid(2);
            [yp,xp] = ind2sub( E,conn_comp.PixelIdxList{r});
            P = [xp,yp]';
            C = [xc,yc]';
            distance = slmetric_pw(P,C, 'eucdist');
            [~,chosen] = min( distance );
            
            tmp = zeros( conn_comp.ImageSize );
            tmp( conn_comp.PixelIdxList{r}(chosen) ) =1;
            
            mask_word(:,:,w,u) = mask_word(:,:,w,u) + tmp;
        end
    end
end

limits = prctile(prob(:),[10 40 60 70 75 90 95 97 99]);
font_sizes(1) = 0;
% OK 
%limits = prctile(prob(:),[40 55 70 85 90 98]);
idfond = sum( bsxfun( @ge, prob(:), limits ),2)+1;
fonts_out = font_sizes( reshape( idfond, [E,options.no_words]) );

space_between_words = 1;
WO = cell(E);
FO = zeros([E,options.no_words]);
for y=1:E(1)
    for x=1:E(2)
        
        [values_sorted,indeces_sorted] = sort( squeeze( pi(y,x,:) )'.*options.freq ,'descend');
        words_to_show = {vocabulary{indeces_sorted(1:options.no_words)}};
        fonts = squeeze( fonts_out(y,x,:) );
        WO{y,x} = words_to_show;
        FO(y,x,:) = fonts; 
        
        pixels_fonts = zeros(1,options.no_words );
        for w=1:options.no_words
            id = find( font_to_pixels(:,1) == fonts(w)) ;
            if ~isempty( id )
                pixels_fonts(w) = font_to_pixels(id,2);
            else
                gamma = 10;
                average =  exp( -gamma*abs( font_to_pixels(:,1) - fonts(w) ) ) ./ sum( exp( -gamma*abs( font_to_pixels(:,1) - fonts(w) ) ) );
                pixels_fonts(w) = average'*font_to_pixels(:,2);
            end
        end
        
        height_text_box = sum(pixels_fonts)+space_between_words*(options.no_words-1); % C
        upper_left = [limx(x),limy(y)]; % in x-y 
        offset = ( size_cell(y,x,1) - height_text_box ) ./ 2 ; 
        
        for ww=1:options.no_words
            
            % Conrolla se la parola e' stata usata in un intorno... - Check if the word was used in a circle ...
            index_us = find( strcmp( used, words_to_show{ww} )==1);
            
            % wgt = 1./(exp( 0*values_sorted(ww) ));
            if mask_word(y,x,ww,index_us) ~=0
                
                % Questa sotto funziona! - This works under!
                % text( midpoint_cell(y,x,2),ly-midpoint_cell(y,x,1), words_to_show{ww}, 'HorizontalAlignment','center','FontSize',fonts_out(y,x) );

                % Sotto per piu' di una cella! - Below for more than one cell!
                if fonts(ww) > 0
                    offset_word = sum(pixels_fonts(1:ww-1))+ space_between_words*(ww-1);
                    text( midpoint_cell(y,x,2), ly-upper_left(2)-offset-offset_word, words_to_show{ww}, 'HorizontalAlignment','center','VerticalAlignment','Cap','FontSize', fonts(ww) );
                    word_out{y,x} = words_to_show{1};
                end
            end
            
        end
        
    end
end


