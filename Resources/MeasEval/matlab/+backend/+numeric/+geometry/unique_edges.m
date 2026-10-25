function edges = unique_edges(T, X)
%> @file unique_edges.m
%> @brief Extracts unique undirected edges from simplex topology.
%> @param T K-by-(D+1) simplex vertex-index matrix.
%> @param X Optional N-by-D vertex coordinates used for the established length filter.
%> @return M-by-2 sorted unique vertex-index pairs.
%>
%> When X is supplied, edges with squared length at least nine times the mean
%> squared edge length are removed. The function owns topology normalization;
%> it does not perform iso-/level-set intersection and has no side effects.
%GET_UNIQUE_EDGES_FROM_DELAUNAYN Extracts unique undirected edges from a Delaunay triangulation.
%
%   This function returns a list of unique, undirected edges from a Delaunay
%   triangulation. Each edge is defined by two point indices that reference the
%   original set of coordinates X used for triangulation. Edges are treated as 
%   undirected, so the pair (i, j) is considered the same as (j, i).
%
%   The function supports triangulations in arbitrary dimensions, i.e., simplices
%   in 2D (triangles), 3D (tetrahedra), or higher.
%
%   If the point coordinates X are provided, the function removes edges that are
%   considered too long compared to the average edge length, which helps in 
%   filtering out invalid or degenerate elements.
%
% INPUT:
%   T : K ?? (D+1) matrix containing indices of the D-dimensional simplices.
%       Each row defines one simplex with (D+1) point indices.
%
%   X : (Optional) M ?? D matrix of point coordinates. If provided, edges that 
%       exceed a multiple of the average length are excluded.
%
% OUTPUT:
%   edges : N ?? 2 matrix of unique, undirected edges (pairs of point indices).
%
% EXAMPLE:
%   X = rand(30, 3);            % Random 3D points
%   T = delaunayn(X);           % Delaunay triangulation (tetrahedra)
%   edges = backend.numeric.geometry.unique_edges(T, X);
%
%   figure;
%   tetramesh(T, X, 'FaceAlpha', 0.01); % Transparent tetrahedra
%   hold on;
%   scatter3(X(:,1), X(:,2), X(:,3), 8, 'k', 'filled'); % Points
%   axis equal;
%   xlabel('X'); ylabel('Y'); zlabel('Z');
%   title('3D Delaunay Tetrahedral Mesh with Unique Edges');
%
%   % Plot edges
%   for k = 1:size(edges, 1)
%       i = edges(k, 1);
%       j = edges(k, 2);
%       plot3([X(i,1), X(j,1)], [X(i,2), X(j,2)], [X(i,3), X(j,3)], 'r-.');
%   end

    if nargin == 2
        remove_invalid_edges = true;
        validity_range = 9;  % Threshold multiplier for length filtering
    else
        remove_invalid_edges = false;
    end

    % Number of points per simplex (e.g., 4 for tetrahedra)
    num_points_per_simplex = size(T, 2);

    % All combinations of 2 vertices per simplex (i.e., all potential edges)
    combs = nchoosek(1:num_points_per_simplex, 2);

    % Collect all edges across all simplices
    all_edges = [];
    for i = 1:size(combs, 1)
        edge_pair = T(:, combs(i, :)); % Edge indices from all simplices
        all_edges = [all_edges; edge_pair];
    end

    % Sort each edge so (i,j) == (j,i), i.e., edges are undirected
    all_edges = sort(all_edges, 2);

    % Keep only unique edges
    edges = unique(all_edges, 'rows');

    % If coordinates are provided, filter out edges that are too long
    if remove_invalid_edges
        start_points = X(edges(:,1), :);
        end_points   = X(edges(:,2), :);
        edge_vectors = end_points - start_points;
        length_squared = sum(edge_vectors.^2, 2);
        mean_length_squared = mean(length_squared);
        is_valid = length_squared < validity_range * mean_length_squared;
        edges = edges(is_valid, :);
    end
end
