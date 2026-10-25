function [valid_vectors, mean_edge_length] = valid_edges(X, DT, validity_range)
%> @file valid_edges.m
%> @brief Classifies simplex edges by an explicit relative-length threshold.
%> @param X N-by-D vertex/sample matrix; the first simplex-dimension columns are coordinates.
%> @param DT M-by-K simplex vertex-index matrix.
%> @param validity_range Optional multiplier of mean edge length; default 3.
%> @return valid_vectors M-by-nchoosek(K,2) logical edge-validity matrix.
%> @return mean_edge_length Root-mean-square simplex edge length in coordinate units.
%>
%> The function owns mesh quality/topology classification only. It neither
%> extracts level sets nor performs file/UI operations.
%% GET_VALID_TRIANGULATION_VECTORS
% Identifies "valid" edge vectors within a triangulation.
%
%   valid_vectors = GET_VALID_TRIANGULATION_VECTORS(X, DT)
%   Returns a logical matrix indicating which edge vectors in the triangulation
%   are considered "valid" based on their squared lengths being below a threshold.
%   By default, the threshold is 3 times the mean edge length.
%
%   valid_vectors = GET_VALID_TRIANGULATION_VECTORS(X, DT, validity_range)
%   Allows specification of a custom threshold multiplier `validity_range`.
%
% PARAMETERS:
%   X              - [NxD] matrix of vertex coordinates.
%                    Only the first D-1 columns are used for the distance computation,
%                    corresponding to the coordinates used for the triangulation.
%
%   DT             - [MxK] matrix of indices defining M simplices with K vertices each.
%                    Typically obtained from a Delaunay triangulation.
%
%   validity_range - (Optional) Scalar multiplier defining the threshold for
%                    valid edge lengths. The threshold is set to
%                    `validity_range^2 * mean(edge_length)^2`.
%                    Default is 3.
%
% RETURNS:
%   valid_vectors  - [MxE] logical matrix, where E is the number of edges per simplex.
%                    Each element is `true` if the squared length of the corresponding
%                    edge is less than the squared threshold.
%
% EXAMPLE:
%   X = randn(100, 3);                              
%   DT = delaunayn(X(:,1:2));                              
%   valid_vectors = backend.numeric.geometry.valid_edges(X, DT);


    % Set default threshold if not provided
    if nargin < 3 
        validity_range = 3; 
    end 

    N_dim = size(DT,2);               % Number of vertices per simplex
    vertex_pairs = nchoosek(1:N_dim, 2);  % All unique edge pairs in each simplex
    d = N_dim - 1;                    % Dimensionality of space used (assumes embedded in d+1 space)
    num_simplices = size(DT,1);      % Number of simplices in triangulation
    num_edges_per_simplex = size(vertex_pairs,1); % Number of edges per simplex

    % Preallocate matrix for squared edge lengths
    edge_lengths_squared = zeros(num_simplices, num_edges_per_simplex);

    % Loop over all simplices
    for k = 1:num_simplices
        simplex_pts = X(DT(k,:), 1:d);  % Get coordinates of simplex vertices (only first d dims)
        for e = 1:num_edges_per_simplex
            i = vertex_pairs(e,1);     % Index of first vertex of edge
            j = vertex_pairs(e,2);     % Index of second vertex of edge
            vec = simplex_pts(i,:) - simplex_pts(j,:);  % Edge vector
            edge_lengths_squared(k,e) = vec * vec';     % Squared length
        end
    end

    % Determine validity by comparing to scaled mean edge length squared    
    mean_edge_length = sqrt(mean(edge_lengths_squared(:)));
    threshold = (validity_range * mean_edge_length).^2;
    valid_vectors = edge_lengths_squared < threshold;
end
