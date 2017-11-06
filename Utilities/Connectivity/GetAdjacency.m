%GETADJACENCY computes the adjacency matrix
%   G = GETADJACENCY(lambda) computes the Adjacency matrix based on the
%   connectivity of a kinematic structure as described by lambda.
%   No circular chains are allowed in this function. Also, lambda should
%   contain the base with identifier 0 because the base will be included in
%   the adjacency matrix.
%
%   Made by Erik Vlasblom
%   Last modified: 01-04-2015
function G = GetAdjacency(lambda)

G = zeros(length(lambda)+1); % add one more for the base with id 0
for ii = 1:length(lambda)
    jj = lambda{ii};
    G(ii+1,jj+1) = 1;
    G(jj+1,ii+1) = 1;
end

end