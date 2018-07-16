% Main script
addpath('./Code')

input_file_name='dictionaryVerySmallSample.txt';
output_file_name='dictionaryVerySmallSample';

%%%%%%%%%%%%%%%% Hyper parameters
includeTitles=1;    % Include the words from titles in the bags of words
includeAbstracts=1; % Include the words from abstracts in the bags

minWordFreq=1;  % Limit the vocab to words that occur at least this many times
minDocLength=3; % Only documents with at least this many words will included 

wd_size = [5 5]; % window size
cg_size = [32,32]; % grid size
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Tokenization
% The input file in Input/<input_file_name>
% consists of three tab delimited text columns: titles, abstracts and links
% E.g.: 
% titles{1}='Counting grids'; 
% abstracts{1}='Our model embeds bags of words into a grid of word counts'
% links{1}='https://www.microsoft.com/en-us/research/people/jojic/'

mkBagsOfWords; 

%%%%%%%%%%%%%%%%%%% CG learning and browser file computation
% The intermediate files are stored in the Model dir. 
% The final browser files are stored in the Data subdir of the Browser
% directiory

mkBrowserFiles; 
