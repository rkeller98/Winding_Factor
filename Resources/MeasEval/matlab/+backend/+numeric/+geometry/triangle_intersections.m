
function [Result_vec] = triangle_intersections(arg1, arg2, arg3, arg4, arg5, T)
%> @file triangle_intersections.m
%> @brief Intersects a plane with supplied triangular topology.
%> @param arg1 N-by-M samples; columns 1:3 are geometry and later columns are interpolants.
%> @param arg2 Plane coefficient a in a*x+b*y+c*z=d.
%> @param arg3 Plane coefficient b in a*x+b*y+c*z=d.
%> @param arg4 Plane coefficient c in a*x+b*y+c*z=d.
%> @param arg5 Plane offset d, in the combined coordinate units.
%> @param T K-by-3 vertex-index matrix supplied by triangulation code.
%> @return Unique plane-intersection samples with additional columns interpolated.
%>
%> The routine rejects anomalously long triangle vectors using its established
%> three-times-mean heuristic. It owns level-set/plane extraction, not topology,
%> and has no external side effects.
% INPUTS
%   arg1   n x m matrix with measured values that should be interpolated
%   n: amount of measurements performed
%   m: amount of different input variables, the first two columns are the  x and y
%   values of the DoE plane, the thrid column contains the variable where
%   that defines plane intersections (equipotential lines), all additional
%   columns are measured values that should be interpolated with respect to
%   the calculated equipotential lines (amount can vary)

%   arg2   a of parameter plane equation, default = 0
%   arg3   b of parameter plane equation, default = 0
%   arg4   c of parameter plane equation, default = 1
%   arg5   d of parameter plane equation, default = value of equipotential
%   line
%   T:     Delaunay Matrix

%% Check for Validity
validity_range = 3;

% Length of T
amount_tri = length(T(:,1));

% plane form ax + by + cz = d
a = arg2;
b = arg3;
c = arg4;
d = arg5;

a_vector = arg2 * ones(amount_tri,1);
b_vector = arg3 * ones(amount_tri,1);
c_vector = arg4 * ones(amount_tri,1);
d_vector = arg5 * ones(amount_tri,1);


P1_matrix = [arg1(T(:,1),1) arg1(T(:,1),2) arg1(T(:,1),3)];
P2_matrix = [arg1(T(:,2),1) arg1(T(:,2),2) arg1(T(:,2),3)];
P3_matrix = [arg1(T(:,3),1) arg1(T(:,3),2) arg1(T(:,3),3)];


% vector 1: P1 --> P2
Vector_1 = P2_matrix - P1_matrix;
% vector 2: P2 --> P3
Vector_2 = P3_matrix - P2_matrix;
% vector 3: P3 --> P1
Vector_3 = P1_matrix - P3_matrix;

% Aux_P12: Half Vector 1 + P_1
Aux_P_12 = P1_matrix + 0.5 * Vector_1;
% Aux_P23: Half Vector 2 + P_2
Aux_P_23 = P2_matrix + 0.5 * Vector_2;
% Aux_P31: Half Vector 3 + P_3
Aux_P_31 = P3_matrix + 0.5 * Vector_3;

Aux_Vector_12_3 = P3_matrix - Aux_P_12;
Aux_Vector_23_1 = P1_matrix - Aux_P_23;
Aux_Vector_31_2 = P2_matrix - Aux_P_31;

Aux_Vector_12_23 =  Aux_P_23 - Aux_P_12;
Aux_Vector_23_31 =  Aux_P_31 - Aux_P_23;
Aux_Vector_31_12 =  Aux_P_12 - Aux_P_31;

norm_vector_1           = sqrt(Vector_1(:,1).^2 + Vector_1(:,2).^2 + Vector_1(:,3).^2);
norm_vector_2           = sqrt(Vector_2(:,1).^2 + Vector_2(:,2).^2 + Vector_2(:,3).^2);
norm_vector_3           = sqrt(Vector_3(:,1).^2 + Vector_3(:,2).^2 + Vector_3(:,3).^2);
norm_vector_A_12_3      = sqrt(Aux_Vector_12_3(:,1).^2 + Aux_Vector_12_3(:,2).^2 + Aux_Vector_12_3(:,3).^2);
norm_vector_A_23_1      = sqrt(Aux_Vector_23_1(:,1).^2 + Aux_Vector_23_1(:,2).^2 + Aux_Vector_23_1(:,3).^2);
norm_vector_A_31_2      = sqrt(Aux_Vector_31_2(:,1).^2 + Aux_Vector_31_2(:,2).^2 + Aux_Vector_31_2(:,3).^2);
norm_vector_A_12_23     = sqrt(Aux_Vector_12_23(:,1).^2 + Aux_Vector_12_23(:,2).^2 + Aux_Vector_12_23(:,3).^2);
norm_vector_A_23_31     = sqrt(Aux_Vector_23_31(:,1).^2 + Aux_Vector_23_31(:,2).^2 + Aux_Vector_23_31(:,3).^2);
norm_vector_A_31_12     = sqrt(Aux_Vector_31_12(:,1).^2 + Aux_Vector_31_12(:,2).^2 + Aux_Vector_31_12(:,3).^2);

index_norm_vector = find(norm_vector_1 >...
    validity_range * mean(norm_vector_1));
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_2 > validity_range * mean(norm_vector_2))];
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_3 > validity_range * mean(norm_vector_3))];
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_A_12_3 > validity_range * mean(norm_vector_A_12_3))];
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_A_23_1 > validity_range * mean(norm_vector_A_23_1))];
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_A_31_2 > validity_range * mean(norm_vector_A_31_2))];
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_A_12_23 > validity_range * mean(norm_vector_A_12_23))];
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_A_23_31 > validity_range * mean(norm_vector_A_23_31))];
index_norm_vector = [index_norm_vector; ...
    find(norm_vector_A_31_12 > validity_range * mean(norm_vector_A_31_12))];

index_norm_vector = unique(index_norm_vector);

% Delete Hypersimplex that are not valid in all data points
T(index_norm_vector,:) = []; 


%% INIT
% Length of T
amount_tri = length(T(:,1));

% plane form ax + by + cz = d
a = arg2;
b = arg3;
c = arg4;
d = arg5;

a_vector = arg2 * ones(amount_tri,1);
b_vector = arg3 * ones(amount_tri,1);
c_vector = arg4 * ones(amount_tri,1);
d_vector = arg5 * ones(amount_tri,1);


P1_matrix = [arg1(T(:,1),1) arg1(T(:,1),2) arg1(T(:,1),3)];
P2_matrix = [arg1(T(:,2),1) arg1(T(:,2),2) arg1(T(:,2),3)];
P3_matrix = [arg1(T(:,3),1) arg1(T(:,3),2) arg1(T(:,3),3)];


% vector 1: P1 --> P2
Vector_1 = P2_matrix - P1_matrix;
% vector 2: P2 --> P3
Vector_2 = P3_matrix - P2_matrix;
% vector 3: P3 --> P1
Vector_3 = P1_matrix - P3_matrix;

% Aux_P12: Half Vector 1 + P_1
Aux_P_12 = P1_matrix + 0.5 * Vector_1;
% Aux_P23: Half Vector 2 + P_2
Aux_P_23 = P2_matrix + 0.5 * Vector_2;
% Aux_P31: Half Vector 3 + P_3
Aux_P_31 = P3_matrix + 0.5 * Vector_3;

Aux_Vector_12_3 = P3_matrix - Aux_P_12;
Aux_Vector_23_1 = P1_matrix - Aux_P_23;
Aux_Vector_31_2 = P2_matrix - Aux_P_31;

Aux_Vector_12_23 =  Aux_P_23 - Aux_P_12;
Aux_Vector_23_31 =  Aux_P_31 - Aux_P_23;
Aux_Vector_31_12 =  Aux_P_12 - Aux_P_31;


lambda_matrix(:,1) = (d_vector - (P1_matrix(:,1).* a_vector + ...
    P1_matrix(:,2).* b_vector + P1_matrix(:,3).* c_vector)) ./ (...
    Vector_1(:,1).* a_vector + ...
    Vector_1(:,2).* b_vector + Vector_1(:,3).* c_vector);

lambda_matrix(:,2) = (d_vector - (P2_matrix(:,1).* a_vector + ...
    P2_matrix(:,2).* b_vector + P2_matrix(:,3).* c_vector)) ./ (...
    Vector_2(:,1).* a_vector + ...
    Vector_2(:,2).* b_vector + Vector_2(:,3).* c_vector);

lambda_matrix(:,3) = (d_vector - (P3_matrix(:,1).* a_vector + ...
    P3_matrix(:,2).* b_vector + P3_matrix(:,3).* c_vector)) ./ (...
    Vector_3(:,1).* a_vector + ...
    Vector_3(:,2).* b_vector + Vector_3(:,3).* c_vector);

lambda_matrix(:,4) = (d_vector - (Aux_P_12(:,1).* a_vector + ...
    Aux_P_12(:,2).* b_vector + Aux_P_12(:,3).* c_vector)) ./ (...
    Aux_Vector_12_3(:,1).* a_vector + ...
    Aux_Vector_12_3(:,2).* b_vector + Aux_Vector_12_3(:,3).* c_vector);

lambda_matrix(:,5) = (d_vector - (Aux_P_23(:,1).* a_vector + ...
    Aux_P_23(:,2).* b_vector + Aux_P_23(:,3).* c_vector)) ./ (...
    Aux_Vector_23_1(:,1).* a_vector + ...
    Aux_Vector_23_1(:,2).* b_vector + Aux_Vector_23_1(:,3).* c_vector);

lambda_matrix(:,6) = (d_vector - (Aux_P_31(:,1).* a_vector + ...
    Aux_P_31(:,2).* b_vector + Aux_P_31(:,3).* c_vector)) ./ (...
    Aux_Vector_31_2(:,1).* a_vector + ...
    Aux_Vector_31_2(:,2).* b_vector + Aux_Vector_31_2(:,3).* c_vector);

lambda_matrix(:,7) = (d_vector - (Aux_P_12(:,1).* a_vector + ...
    Aux_P_12(:,2).* b_vector + Aux_P_12(:,3).* c_vector)) ./ (...
    Aux_Vector_12_23(:,1).* a_vector + ...
    Aux_Vector_12_23(:,2).* b_vector + Aux_Vector_12_23(:,3).* c_vector);

lambda_matrix(:,8) = (d_vector - (Aux_P_23(:,1).* a_vector + ...
    Aux_P_23(:,2).* b_vector + Aux_P_23(:,3).* c_vector)) ./ (...
    Aux_Vector_23_31(:,1).* a_vector + ...
    Aux_Vector_23_31(:,2).* b_vector + Aux_Vector_23_31(:,3).* c_vector);

lambda_matrix(:,9) = (d_vector - (Aux_P_31(:,1).* a_vector + ...
    Aux_P_31(:,2).* b_vector + Aux_P_31(:,3).* c_vector)) ./ (...
    Aux_Vector_31_12(:,1).* a_vector + ...
    Aux_Vector_31_12(:,2).* b_vector + Aux_Vector_31_12(:,3).* c_vector);


Intersection_matrix_1 = P1_matrix + lambda_matrix(:,1).* Vector_1;
Intersection_matrix_2 = P2_matrix + lambda_matrix(:,2).* Vector_2;
Intersection_matrix_3 = P3_matrix + lambda_matrix(:,3).* Vector_3;
Intersection_matrix_4 = Aux_P_12 + lambda_matrix(:,4).* Aux_Vector_12_3;
Intersection_matrix_5 = Aux_P_23 + lambda_matrix(:,5).* Aux_Vector_23_1;
Intersection_matrix_6 = Aux_P_31 + lambda_matrix(:,6).* Aux_Vector_31_2;
Intersection_matrix_7 = Aux_P_12 + lambda_matrix(:,7).* Aux_Vector_12_23;
Intersection_matrix_8 = Aux_P_23 + lambda_matrix(:,8).* Aux_Vector_23_31;
Intersection_matrix_9 = Aux_P_31 + lambda_matrix(:,9).* Aux_Vector_31_12;

% Evaluate if intersection Point is within Delaunay triangle (Interpolation
% can be applied if lambda is between 0 and 1)
Result_index_1 = find(lambda_matrix(:,1) <=1 & lambda_matrix(:,1) >= 0);
Result_index_2 = find(lambda_matrix(:,2) <=1 & lambda_matrix(:,2) >= 0);
Result_index_3 = find(lambda_matrix(:,3) <=1 & lambda_matrix(:,3) >= 0);
Result_index_4 = find(lambda_matrix(:,4) <=1 & lambda_matrix(:,4) >= 0);
Result_index_5 = find(lambda_matrix(:,5) <=1 & lambda_matrix(:,5) >= 0);
Result_index_6 = find(lambda_matrix(:,6) <=1 & lambda_matrix(:,6) >= 0);
Result_index_7 = find(lambda_matrix(:,7) <=1 & lambda_matrix(:,7) >= 0);
Result_index_8 = find(lambda_matrix(:,8) <=1 & lambda_matrix(:,8) >= 0);
Result_index_9 = find(lambda_matrix(:,9) <=1 & lambda_matrix(:,9) >= 0);

% Assign results

Result_vec = Intersection_matrix_1(Result_index_1,:,:);
Result_vec = [Result_vec; Intersection_matrix_2(Result_index_2,:,:)];
Result_vec = [Result_vec; Intersection_matrix_3(Result_index_3,:,:)];
Result_vec = [Result_vec; Intersection_matrix_4(Result_index_4,:,:)];
Result_vec = [Result_vec; Intersection_matrix_5(Result_index_5,:,:)];
Result_vec = [Result_vec; Intersection_matrix_6(Result_index_6,:,:)];
Result_vec = [Result_vec; Intersection_matrix_7(Result_index_7,:,:)];
Result_vec = [Result_vec; Intersection_matrix_8(Result_index_8,:,:)];
Result_vec = [Result_vec; Intersection_matrix_9(Result_index_9,:,:)];



%% Interpolation of additional values
% The interpolation of these values has to be performed but that are not
% subject of the search condition, first thee colums of arg1 are x,y and z
% values, the additional colums are values that should be interpolated for
% the values that fullfill z = arg5

if size(arg1,2) > 3
    for count_interpolants = 4:size(arg1,2)
        
        P1_matrix = [arg1(T(:,1),1) arg1(T(:,1),2) arg1(T(:,1), ...
            count_interpolants)];
        P2_matrix = [arg1(T(:,2),1) arg1(T(:,2),2) arg1(T(:,2), ...
            count_interpolants)];
        P3_matrix = [arg1(T(:,3),1) arg1(T(:,3),2) arg1(T(:,3), ...
            count_interpolants)];
        
        
        % vector 1: P1 --> P2
        Vector_1 = P2_matrix - P1_matrix;
        % vector 2: P2 --> P3
        Vector_2 = P3_matrix - P2_matrix;
        % vector 3: P3 --> P1
        Vector_3 = P1_matrix - P3_matrix;
        
        % Aux_P12: Half Vector 1 + P_1
        Aux_P_12 = P1_matrix + 0.5 * Vector_1;
        % Aux_P23: Half Vector 2 + P_2
        Aux_P_23 = P2_matrix + 0.5 * Vector_2;
        % Aux_P31: Half Vector 3 + P_3
        Aux_P_31 = P3_matrix + 0.5 * Vector_3;
        
        Aux_Vector_12_3 = P3_matrix - Aux_P_12;
        Aux_Vector_23_1 = P1_matrix - Aux_P_23;
        Aux_Vector_31_2 = P2_matrix - Aux_P_31;
        
        Aux_Vector_12_23 =  Aux_P_23 - Aux_P_12;
        Aux_Vector_23_31 =  Aux_P_31 - Aux_P_23;
        Aux_Vector_31_12 =  Aux_P_12 - Aux_P_31;
        
        % Use the lambdas that are calculated in the previous part to interpolate
        % the values at the point where the plane intersection are
        Intersection_matrix_1 = P1_matrix + lambda_matrix(:,1).* ...
            Vector_1;
        Intersection_matrix_2 = P2_matrix + lambda_matrix(:,2).* ...
            Vector_2;
        Intersection_matrix_3 = P3_matrix + lambda_matrix(:,3).* ...
            Vector_3;
        Intersection_matrix_4 = Aux_P_12 + lambda_matrix(:,4).* ...
            Aux_Vector_12_3;
        Intersection_matrix_5 = Aux_P_23 + lambda_matrix(:,5).* ...
            Aux_Vector_23_1;
        Intersection_matrix_6 = Aux_P_31 + lambda_matrix(:,6).* ...
            Aux_Vector_31_2;
        Intersection_matrix_7 = Aux_P_12 + lambda_matrix(:,7).* ...
            Aux_Vector_12_23;
        Intersection_matrix_8 = Aux_P_23 + lambda_matrix(:,8).* ...
            Aux_Vector_23_31;
        Intersection_matrix_9 = Aux_P_31 + lambda_matrix(:,9).* ...
            Aux_Vector_31_12;

% Get results        
Temp_result_vec = Intersection_matrix_1(Result_index_1,3,:);
Temp_result_vec = [Temp_result_vec; Intersection_matrix_2(Result_index_2,3,:)];
Temp_result_vec = [Temp_result_vec; Intersection_matrix_3(Result_index_3,3,:)];
Temp_result_vec = [Temp_result_vec; Intersection_matrix_4(Result_index_4,3,:)];
Temp_result_vec = [Temp_result_vec; Intersection_matrix_5(Result_index_5,3,:)];
Temp_result_vec = [Temp_result_vec; Intersection_matrix_6(Result_index_6,3,:)];
Temp_result_vec = [Temp_result_vec; Intersection_matrix_7(Result_index_7,3,:)];
Temp_result_vec = [Temp_result_vec; Intersection_matrix_8(Result_index_8,3,:)];
Temp_result_vec = [Temp_result_vec; Intersection_matrix_9(Result_index_9,3,:)];  

% Assign Results to Ouput vector
Result_vec = [Result_vec Temp_result_vec];
        
    end
end

Result_vec = unique(Result_vec,'Rows');



