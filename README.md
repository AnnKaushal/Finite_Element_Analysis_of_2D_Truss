# Finite_Element_Analysis_of_2D_Truss

## Description

This MATLAB project demonstrates finite element analysis on a 2D truss structure. The analysis involves discretization of the truss into nodes and elements, calculation of the global stiffness matrix, solving for nodal displacements, computing element forces and support reactions, and finally visualizing the results.

## Author

- **Author:** Anandita Kaushal
- **Contact:** 20bce034@nith.ac.in

---

## Getting Started

To run this project, follow the steps below:

1. Open MATLAB and create a new script or function file.

2. Copy and paste the code from this README into the script or function file.

3. Run the script or function file to perform finite element analysis on your 2D truss structure.

---

## Step 1: Discretization of a 2D Truss Structure

In this step, the user is prompted to provide input for the number of nodes and elements in the truss. The node coordinates, element connectivity, and element properties are obtained from the user. The nodes and elements are then plotted to visualize the truss geometry.

**Usage:**

- Input the number of nodes.
- Input the number of elements.
- Enter the coordinates of each node (x, y).
- Enter the connectivity of each element (start node, end node).
- Provide element properties (cross-sectional area, modulus of elasticity).

---

## Step 2: Computing Global Stiffness Matrix and External Force Vector

The external loads at each node are input by the user. The global stiffness matrix and external force vector are computed based on element properties and nodal displacements. These matrices and vectors are used in the later steps to solve for displacements and analyze forces.

**Usage:**

- Input the external loads at each node (Fx, Fy).

---

## Step 3: Computing Axial Displacements, Element Forces & Support Reactions

Using the global stiffness matrix and external force vector, the axial displacements of nodes are calculated. Element forces and support reactions are also determined. Stress in each element is computed, and the modified displacement matrix is displayed, accounting for missing values.

**Usage:**

- No additional user input required.

---

## Step 4: Plotting the Results

The project concludes with visualizations of the analysis results. The global stiffness matrix and external force vector are plotted. Stress distribution in each element and support reactions are also displayed. The original and deformed truss structures are visualized side by side, showing the effects of nodal displacements.

**Usage:**

- No additional user input required.

---

## Results

- Global Stiffness Matrix K
- External Force Vector F
- Stress in Each Element
- Support Reactions
- Displacements

---

## References

- No external references are required for this project as it is a self-contained MATLAB script for finite element analysis of 2D truss structures.

---

**Note:** Please ensure that you have MATLAB installed to run this project successfully.
