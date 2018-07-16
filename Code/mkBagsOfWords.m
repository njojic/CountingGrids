% This scripts turns the input file into bags of words
% A document is considered a combination of title, abstract and a weblink
% The processed data is stored in alldata.mat

%%%%%%%%%%%%%%%% Hyper parameters
%includeTitles=1;    % Include the words from titles in the bags of words
%includeAbstracts=1; % Include the words from abstracts in the bags

%minWordFreq=3;  % Limit the vocab to words that occur at least this many times
%minDocLength=5; % Only documents with at least this many words will included 
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Load Data 
loadInput;     % input.txt has 3 tab delimited columns; 
               % loaded into cell array input (titles, abstracts, links)
%%%%%%%%%%%%%%%%

titles=input(:,1);
if size(input,2)>1
    abstracts=input(:,2);
else
    abstracts=cell(size(input,1),1);
end
if size(input,2)>2
links=input(:,3);
else
    links=cell(size(input,1),1);
end

Ndocs=length(abstracts);
absmeta=abstracts; % will include titles, abstracts, or both

for t=1:Ndocs
    if includeTitles
        absmeta{t}=titles{t};
    end
    if includeAbstracts
        absmeta{t}=[ absmeta{t},' ',abstracts{t},' ###docend### '];
    end 
    absmeta{t}(absmeta{t}<33)=' '; absmeta{t}(absmeta{t}>126)=' ';
end

tmp=char(absmeta)';
tmp(end+1,:)=32;
unwData=lower(strsplit(char(tmp(:)')));
[VOCbig,~,tokenData] = unique(unwData);

ind_docend=find(strcmp(unwData,'###docend###'));
countsBig=sparse(length(VOCbig),Ndocs);
dBegin=1;
for t=1:Ndocs
    ind=tokenData(dBegin:ind_docend(t));
    dBegin=ind_docend(t)+1;
    countsBig(ind,t)=1;
end
    
tmp=sum(countsBig,2);
[~, fr_ind]=sort(tmp,'Descend');
counts=countsBig(fr_ind(2:end),:); %get rid of docend
VOC=VOCbig(fr_ind(2:end));


% VOC contains the vocabulary, and counts contains the word counts
% The following code cleans these up

nhFlag=zeros(length(VOC),1);
for i=1:length(VOC)
    if length(VOC{i}>0)
    nhFlag(i)=VOC{i}(1)=='@'; %depending on the CG browser, @ may break it
    end
end

ind=find(nhFlag==0);
counts=counts(ind,:);
VOC=VOC(ind);

sz=[0 0];
while sz(1)~=size(counts,1) || sz(2)~=size(counts,2)
    sz=size(counts);
    ind=find(sum(counts,2)>=minWordFreq);
    counts=counts(ind,:);
    VOC=VOC(ind);

    ind=find(sum(counts,1)>=minDocLength);
    abstracts=abstracts(ind);
    links=links(ind);
    titles=titles(ind);
    counts=counts(:,ind);
end

save -v7.3 Model\alldata abstracts titles links counts  VOC
% counts is a words X document matrix of word counts;
% VOC is the vocabulary; size(counts,1) is equal to length(VOC)

