function [maxs,locs]=maxNvalues(A,n,Unique)
% this function finds the top n values in a multidimension array
% A, matrix to be searched
% n, number of desired values
% Unique, (Optional) true if only unique maximums are desired (first
% location, as defined in column format, so it may not be what you think)
% max are the values returned and locs are their locations in A

% based on and expanded from
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/238315
% but I though it would be handy, so I put it in my `personal' toolbox 

if isempty(A)
    error('Hey, how am I supposed to find a max in an empty matrix?')
end

if nargin<=2||Unique==false
    if nargin<2 % just return one value
        [maxs,temp_1]=max(A(:));
        temp_2=cell(length(size(A)),1);
        [temp_2{:}]=ind2sub(size(A),temp_1);
        locs=cat(1,temp_2{:}).';
        return
    end
if numel(A)<n
    error('Not enough entries in data for outputs requested.')
end
[a,ix] = sort(A(:),'descend');
maxs = a(1:n);
ix = ix(1:n);
temp_2=cell(length(size(A)),1);
[temp_2{:}]=ind2sub(size(A),ix);
locs=cat(2,temp_2{:});
else %Unique must exist and is true
[a,temp_3,dummy_1]=unique(A(:),'first');%#ok<NASGU>,3rd output is required
if numel(a)<n
    error('Not enough unique entries in data for outputs requested.')
end
% unique sorts in ascending order
maxs = a(end:-1:(end-n));
temp_2=cell(length(size(A)),1);
[temp_2{:}]=ind2sub(size(A),temp_3(end:-1:(end-n)));
locs=cat(2,temp_2{:});
end



