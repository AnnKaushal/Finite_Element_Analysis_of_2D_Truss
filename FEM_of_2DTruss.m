% -------------------------------------------------------------------------
% Project Title: Finite Element Analysis of a 2D Truss Structure
%
% Description: This MATLAB project demonstrates finite element analysis on
% a 2D truss structure. The analysis involves discretization of
% the truss into nodes and elements, calculation of the global stiffness
% matrix, solving for nodal displacements, computing element forces and
% support reactions, and finally visualizing the results.
%
% Authors: Anandita Kaushal
% Contact: 20bce034@nith.ac.in
% -------------------------------------------------------------------------

% Get started on a clean command window 
clc

% Step 1: Discretization of a 2D Truss Structure

% In this step, the user is prompted to provide input for the number of nodes
% and elements in the truss. The node coordinates, element connectivity, and
% element properties are obtained from the user. The nodes and elements are
% then plotted to visualize the truss geometry.

% Ask the user for input
numNodes = input('Enter the number of nodes: ');
numElements = input('Enter the number of elements: ');

% Initialize arrays to store node coordinates, element connectivity & element properties
nodeCoordinates = zeros(numNodes, 2);
elementConnectivity = zeros(numElements, 2);
elementProperties = zeros(numElements, 2);  

% Input node coordinates
disp('Enter the coordinates of nodes (x, y):');
for i = 1:numNodes
    fprintf('Node %d:\n', i);
    nodeCoordinates(i, 1) = input('x-coordinate: ');
    nodeCoordinates(i, 2) = input('y-coordinate: ');
end

% Input element connectivity and properties
disp('Enter the connectivity of elements and their properties:');
for i = 1:numElements
    fprintf('Element %d:\n', i);
    elementConnectivity(i, 1) = input('First node: ');
    elementConnectivity(i, 2) = input('Second node: ');
    elementProperties(i, 1) = input('Cross-sectional area: ');
    elementProperties(i, 2) = input('Modulus of elasticity: ');
end
% Display the node coordinates and element connectivity side by side
disp('Node Coordinates:');
disp(nodeCoordinates);

disp('Element Connectivity:');
disp(elementConnectivity);

% Plot the nodes and elements
figure;
hold on;
title('2D Truss Discretization');
xlabel('X-coordinate');
ylabel('Y-coordinate');

% Plot nodes & show them in blue colour
plot(nodeCoordinates(:, 1), nodeCoordinates(:, 2), 'o', 'MarkerFaceColor', 'b');

% Plot elements & show them in red colour
for i = 1:numElements
    startNode = elementConnectivity(i, 1);
    endNode = elementConnectivity(i, 2);
    plot(nodeCoordinates([startNode, endNode], 1), nodeCoordinates([startNode, endNode], 2), 'r');
end

% Set axis limits and display grid
axis equal;
grid on;
hold off;

% Step 2: Computing Global Stiffness Matrix and External Force Vector 

% The external loads at each node are input by the user. The global stiffness
% matrix and external force vector are computed based on element properties
% and nodal displacements. These matrices and vectors are used in the later
% steps to solve for displacements and analyze forces.

% Initialize arrays to store element properties and external loads
% External loads at each node [Fx, Fy]
externalLoads = zeros(numNodes, 2);         

% Input external loads at each node
disp('Enter the external loads at each node (Fx, Fy):');
for i = 1:numNodes
    fprintf('Node %d:\n', i);
    externalLoads(i, 1) = input('Fx: ');
    externalLoads(i, 2) = input('Fy: ');
end

% Calculate total number of degrees of freedom (DOFs)
totalDOFs = 2 * numNodes;

% Initialize arrays for stiffness matrix and force vector
K = zeros(totalDOFs);    % Global stiffness matrix
F = zeros(totalDOFs, 1);  % External force vector

% Loop over each element
for i = 1:numElements
    startNode = elementConnectivity(i, 1);
    endNode = elementConnectivity(i, 2);
    
    % Calculate element length
    dx = nodeCoordinates(endNode, 1) - nodeCoordinates(startNode, 1);
    dy = nodeCoordinates(endNode, 2) - nodeCoordinates(startNode, 2);
    L = sqrt(dx^2 + dy^2);
    
    % Calculate cosine and sine of element angle
    c = dx / L;
    s = dy / L;
    
    % Calculate element stiffness matrix
    EA_over_L = elementProperties(1) * elementProperties(2) / L;
    k = [c^2, c*s, -c^2, -c*s;
         c*s, s^2, -c*s, -s^2;
         -c^2, -c*s, c^2, c*s;
         -c*s, -s^2, c*s, s^2] * EA_over_L;
    
    % Assemble element stiffness into global stiffness matrix
    dofIndices = [2*startNode-1, 2*startNode, 2*endNode-1, 2*endNode];
    for j = 1:4
        for l = 1:4
            K(dofIndices(j), dofIndices(l)) = K(dofIndices(j), dofIndices(l)) + k(j, l);
        end
    end
    
    % Calculate and add external forces to force vector
    F(dofIndices) = F(dofIndices) + [externalLoads(startNode, 1); externalLoads(startNode, 2); externalLoads(endNode, 1); externalLoads(endNode, 2)];
end

% Display the stiffness matrix and force vector
disp('Global Stiffness Matrix K:');
disp(K);

disp('External Force Vector F:');
disp(F);

% Step 3: Computing Axial Displacements, elemnent forces & support reactions

% Using the global stiffness matrix and external force vector, the axial
% displacements of nodes are calculated. Element forces and support reactions
% are also determined. Stress in each element is computed, and the modified
% displacement matrix is displayed, accounting for missing values.

% Solve for displacements using K matrix
displacementsFilled = K ./ F; % Using direct solver to solve linear system

% Fill missing values with zeros

% Replace NaN and Inf values with placeholders
displacementsFilled(isnan(displacementsFilled) | isinf(displacementsFilled)) = 0;

% Display the modified displacement matrix
disp('Displacements:');
disp(displacementsFilled);

% Calculate element forces and support reactions
elementForces = zeros(numElements, 1);
supportReactions = zeros(totalDOFs, 1);

% Initialize arrays for stress in each element and support reactions
Stress = zeros(numElements, 1);
R = zeros(2 * size(nodeCoordinates, 1), 1);

% Loop over each element
for i = 1:numElements
    startNode = elementConnectivity(i, 1);
    endNode = elementConnectivity(i, 2);
    
    dx = nodeCoordinates(endNode, 1) - nodeCoordinates(startNode, 1);
    dy = nodeCoordinates(endNode, 2) - nodeCoordinates(startNode, 2);
    L = sqrt(dx^2 + dy^2);
    
    c = dx / L;
    s = dy / L;
    
    BJ = startNode; % Joint number at the start of the element
    EJ = endNode;   % Joint number at the end of the element
    
    % Calculate stress in the element
    Stress(i) = (elementProperties(2) / L) * (-c * (displacementsFilled(2*BJ-1) - displacementsFilled(2*EJ-1)) - s * (displacementsFilled(2*BJ) - displacementsFilled(2*EJ)));

    % Calculate support reactions
    R = R - K(:, (2 * startNode - 1):(2 * startNode)) * [displacementsFilled(2 * startNode - 1); displacementsFilled(2 * startNode)];
    R = R - K(:, (2 * endNode - 1):(2 * endNode)) * [displacementsFilled(2 * endNode - 1); displacementsFilled(2 * endNode)];
end

% Display stress in each element
disp('Stress in Each Element:');
disp(Stress');

% Display support reactions
disp('Support Reactions:');
disp(R);

% Step 4: Plotting the results

% The project concludes with visualizations of the analysis results. The global
% stiffness matrix and external force vector are plotted. Stress distribution
% in each element and support reactions are also displayed. The original and
% deformed truss structures are visualized side by side, showing the effects
% of nodal displacements.

% Plot the global stiffness matrix K
figure;
subplot(1, 2, 1);
imagesc(K);
title('Global Stiffness Matrix');
xlabel('DOFs');
ylabel('DOFs');
colorbar;

% Plot the external force vector F
subplot(1, 2, 2);
bar(F);
title('Force Vector');
xlabel('DOFs');
ylabel('Force');
grid on;

% Adjust figure layout
set(gcf, 'Position', [100, 100, 1000, 500]);

% Plot stress distribution
figure;
bar(Stress);
title('Element Stress Distribution');
xlabel('Element');
ylabel('Stress');
grid on;

% Plot the original and the deformed truss

% Calculate scaled displacements
scaleFactor = 500;
scaledDisplacements = scaleFactor * displacementsFilled;

% Original truss
figure;
subplot(1, 2, 1);
hold on;
for i = 1:numElements
    startNode = elementConnectivity(i, 1);
    endNode = elementConnectivity(i, 2);
    plot([nodeCoordinates(startNode, 1), nodeCoordinates(endNode, 1)], ...
         [nodeCoordinates(startNode, 2), nodeCoordinates(endNode, 2)], 'b');
end
title('Original Truss');
xlabel('X');
ylabel('Y');
axis equal;

% Deformed truss with nodal displacements
subplot(1, 2, 2);
hold on;

for i = 1:numElements
    startNode = elementConnectivity(i, 1);
    endNode = elementConnectivity(i, 2);
    def_start_x = nodeCoordinates(startNode, 1) + scaleFactor * displacementsFilled(2 * startNode - 1);
    def_start_y = nodeCoordinates(startNode, 2) + scaleFactor * displacementsFilled(2 * startNode);
    def_end_x = nodeCoordinates(endNode, 1) + scaleFactor * displacementsFilled(2 * endNode - 1);
    def_end_y = nodeCoordinates(endNode, 2) + scaleFactor * displacementsFilled(2 * endNode);
    plot([def_start_x, def_end_x], [def_start_y, def_end_y], 'r');
end

title('Deformed Truss with Nodal Displacements');
xlabel('X');
ylabel('Y');
axis equal;

% Adjust figure layout
set(gcf, 'Position', [100, 100, 1000, 600]);
