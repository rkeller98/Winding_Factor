function iso_points = edge_intersections(X, search_dim, iso_val, DT_Edges)
%> @file edge_intersections.m
%> @brief Intersects a scalar level with supplied triangulation edges.
%> @param X N-by-D sample matrix; all columns are linearly interpolated.
%> @param search_dim One-based column containing the scalar level field.
%> @param iso_val Requested scalar iso-value, in the search column's unit.
%> @param DT_Edges M-by-2 vertex-index pairs supplied by triangulation code.
%> @return K-by-D intersection samples with column search_dim equal to iso_val.
%>
%> Endpoints bracket the level when the linear interpolation factor is in
%> [0,1]. This function owns level-set extraction only; mesh topology belongs
%> to backend.numeric.geometry. It has no side effects.
%% backend.numeric.geometry.edge_intersections
% Computes the intersection points of a scalar isovalue along the edges of a
% Delaunay triangulation in arbitrary dimensions.
%
% This function searches all edges of a given Delaunay triangulation and
% performs linear interpolation to find the coordinates where the scalar
% field (represented by one column of the input matrix `X`) crosses a
% specified isovalue.
%
% PARAMETERS:
%   X             - [NxD_total] array of coordinates of N points. The first d columns
%                   are used to define the D-dimensional Delaunay triangulation.
%                   Additional columns may contain scalar or vector fields and are
%                   interpolated accordingly along the edges.
%
%   search_dim    - Scalar index indicating the column in `X` along which the
%                   isovalue should be searched (e.g., 3 for the Z-coordinate).
%                   Must be an integer ??? 1 and ??? size(X,2). Typically, `search_dim`
%                   refers to a column beyond the first d coordinates (i.e., the
%                   part of `X` not used for the triangulation).
%
%   iso_val       - Scalar value specifying the isocontour (e.g., a specific
%                   value of a scalar field along `search_dim`).
%
%   DT_Edges      - [Mx2] array of indices into `X`, where each row defines
%                   an edge (i, j). These edges are typically obtained via the
%                   function `get_unique_edges_from_delaunay`, which extracts the
%                   unique edges of a Delaunay triangulation.
%
% RETURNS:
%   iso_points - [KxD_total] array of interpolated intersection points. Each
%                row corresponds to a location along an edge where the scalar field
%                (column `search_dim`) crosses `iso_val`. The value in `search_dim`
%                is fixed to `iso_val`; all other columns are linearly interpolated.
%
% EXAMPLE USAGE:
%   X = rand(100, 5);                
%   DT = delaunayn(X(:,1:2));        % 2D triangulation based on first two columns
%   edges = get_unique_edges_from_delaunay(DT);
%   iso_pts = backend.numeric.geometry.edge_intersections(X, 3, 0.5, edges);  % Search isocontour in column 3


    % Extract scalar values along the selected dimension for each edge
    target_values = X(:, search_dim);                   % Values along search dimension
    target_points = target_values(DT_Edges);            % Edge endpoints [M x 2]

    % Compute interpolation factors lambda for each edge:
    % lambda = (iso_val - value_at_start) / (value_at_end - value_at_start)
    numerator_diffs   = iso_val - target_points(:, 1);  % ??numerator
    denominator_diffs = target_points(:, 2) - target_points(:, 1);  % ??denominator
    lambda = numerator_diffs ./ denominator_diffs;

    % Identify edges where interpolation yields a valid intersection
    valid_lambda = (lambda >= 0) & (lambda <= 1);

    % Count number of valid intersection points
    num_valid = nnz(valid_lambda);

    % Initialize output matrix for valid iso-points (NaNs by default)
    N_rows = size(X,2);                   % Dimensionality of space (e.g., 3D)
    iso_points = NaN(num_valid, N_rows); % Preallocate with NaNs
    iso_points(:, search_dim) = iso_val; % Assign fixed isovalue in search dimension

    % Determine which other dimensions to interpolate
    dim_vec = 1:N_rows;
    dim_vec(dim_vec == search_dim) = [];  % All dimensions except the search dimension

    % Interpolate remaining dimensions along valid edges
    for dim_idx = dim_vec
        dim_values = X(:, dim_idx);                % All values in current dimension
        dim_points = dim_values(DT_Edges);         % Start and end values per edge
        start_points = dim_points(:, 1);           % Start point of each edge
        direction_vectors = dim_points(:, 2) - start_points;  % Vector from start to end
        interpolated_points = start_points + lambda .* direction_vectors;  % Linear interpolation

        % Keep only interpolated values from valid edges
        iso_points(:, dim_idx) = interpolated_points(valid_lambda);
    end
end
