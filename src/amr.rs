//! `Adaptive Mesh Refinement Module`
//!
//! This is a simple adaptive mesh refinement module. Nothing fancy, but
//! functional. It does not support ghost cells in each branch, nor does it do
//! zone decomposition for MPI parallel computing. For such complexity we
//! refer to e.g. the PARAMESH or CHOMBO library. The objective of the present
//! AMR module is to provide a simple mesh routine which allows to refine
//! anywhere at will by dividing a cell into 2x2x2 (for 3-D) subcells.  In the
//! end this module delivers a base grid in which each grid cell can be either
//! a true cell (leaf) or a tree of sub-cells. All the branches are linked via
//! pointers (a pointer to the parent and pointers to the children) and
//! neighbors are also linked via pointers.  Moreover one can let the module
//! produce a 1-D array of pointers to all cells (leafs), so that referencing
//! to cells can be done in a do-loop rather than requiring the necessity of
//! recursion. A similar thing can be done to all branches (i.e. leafs AND
//! nodes).
//!
//!  The objective of this module is to provide a meshing for e.g.  radiative
//!  transfer or for diffusion-type equations. It is not quite efficient,
//!  however, for the explicit solving of hydrodynamics equations, since
//!  the module is not adapted for parallelization. For that we refer to
//!  the PARAMESH library or CHOMBO library.
//!
//! --------------------------------------------------------------------------
//!  BASE GRID AND TREE STRUCTURED REFINEMENTS
//!
//!  The base grid (forest) with possibly hierarchically structured
//!  refinements (tree). The end-members of the branches are called leafs.
//!  The leafs are the "true" cells! Each symbol (o or *) is a branch
//!  (the amr_branch type). A branch can be a node (consisting of children)
//!  or a leaf (end-member = AMR grid cell).
//!
//!   Base grid (level=0)        o----------*----------o----------o----------o
//!                                       /||\
//!                                     / |  | \
//!   First level refinement           *  o  o  o
//!                                  /||\\
//!                                / |  | \
//!   Second level refinement     o  o  o  o
//!
//!  Two types of branches:
//!    o = leaf
//!    * = node
//!
//!  NOTE: Compared to previous versions of amr_module, the base level is now
//!        0 instead of 1. This is in better comparison to the levels of the
//!        layers, and is in many respects more consistent.
//!
//! --------------------------------------------------------------------------
//!  HOW THE REFEINEMENT WORKS IN 2-D (EXAMPLE)
//!
//!   |-----------|-----------|-----------------------|
//!   |           |           |                       |
//!   |           |           |                       |
//!   |           |           |                       |
//!   |-----|-----|-----------|                       |
//!   |     |     |           |                       |
//!   |-----|-----|           |                       |
//!   |     |     |           |                       |
//!   |-----|-----|-----------|-----------------------|
//!   |                       |                       |
//!   |                       |                       |
//!   |                       |                       |
//!   |                       |                       |
//!   |                       |                       |
//!   |                       |                       |
//!   |                       |                       |
//!   |-----------------------|-----------------------|
//!
//! --------------------------------------------------------------------------
//!  Author: C.P. Dullemond
//!  Date:   07.09.07
//! ==========================================================================
