function [text_clean] = clean_text_single( text )
dbstop if error

T = length( text );

text_clean = cell( size(text));
for t=1:T
    if ~isempty(text{t})
    id = regexp( text{t},'\n');
    text{t}(id) = ' ';
    id = regexp( text{t},'•');
    text{t}(id) = '';    
    id = regexp( text{t},'\t');
    text{t}(id) = '';    
%     Cm = bwconncomp( text{t} == ' ');    
%     modid = find( cellfun( @length, Cm.PixelIdxList) > 1);
%     
%     for m=1:length( modid )
%         text{t}( Cm.PixelIdxList{modid(m)}(2:end)) = '';
%     end
%     
    end
    if isempty (text{t} )
        text_clean{t} = '';
    else
            tmp = text{t};
        id = regexp( tmp, '[a-zA-Z_0-9\-_(){}!*#,:/\s\.]' );
        tmp = regexprep( tmp(id), '&','and');
        text_clean{t} = tmp;
    end
end