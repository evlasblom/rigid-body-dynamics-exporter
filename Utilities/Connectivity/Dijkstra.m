%% DIJKSTRA Computes the shortest path in a graph
% path = DIJKSTRA(G,pair)
%   The Dijkstra algorithm calculates the shortest path from node pair(1)
%   to node pair(2) in the graph described by the Adjecency matrix G.
%   Note that, because a 0 index does not exist, pair can not contain a 0.
%
%   E. Vlasblom, April 2015
%
%   Based on:
%   https://gist.github.com/jcchurch/930280
%   http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm#Pseudocode
%   accessed 01-04-2015
function path = Dijkstra(G, pair)
%#codegen

source = pair(1);
target = pair(2);

vl = size(G,1); % number of vertices
dist = inf(1,vl); % distance from source for each node
prev = -ones(1,vl); % previous node in optimal path for each node
Q = ones(1,vl); % unvisited set of vertices

dist(source) = 0; 
not_seen = vl; 
while not_seen > 0
    [distance, ii] = min(dist .* Q); % least distance in unvisited set
    Q(ii) = inf; % mark visited
    not_seen = not_seen - 1;
    
    if ii == target; break; end % target reached
    if distance == inf; break; end
    
    for jj = 1:vl
        if Q(jj) == 1 && G(ii,jj) > 0 % if unseen neighbour
            alt = distance + G(ii, jj);
            if alt < dist(jj)
                dist(jj) = alt;
                prev(jj) = ii;
            end
        end
    end 
end

% construct shortest path by reverse iteration
path = zeros(1,vl); 
num = 1;
jj = target;
while prev(jj) > 0
    path(num) = jj;
    jj = prev(jj);
    num = num+1;
end
path(num) = source;


end