# pHLA–TCR Docking Angle Definition and Computation

## 1. Purpose and Scope

This document defines the **standardized geometric angles** used to characterize the **global docking geometry** of pHLA–TCR complexes in this software.

The goal of this module is **not** to predict immunogenicity or signaling outcomes, but to provide a **robust, reproducible, and physically interpretable representation of the overall TCR–pHLA docking pose**, suitable for:

- structural quality control of modeled complexes  
- statistical analysis of docking geometry distributions  
- comparison between modeled structures and reference crystal / MD ensembles  
- downstream machine learning feature construction  

Only **global rigid-body–level angles** are considered.

---

## 2. Conceptual Overview

The pHLA–TCR binding geometry is represented using **two complementary angles**:

1. **Crossing angle** — describes the *lateral orientation* of the TCR relative to the pHLA peptide-binding groove  
2. **Incident angle** — describes the *tilt* of the TCR relative to the pHLA groove plane  

Together, these angles define the **global docking pose** of the TCR over the pHLA complex.

This definition follows the **implicit consensus used in structural immunology** and is consistent with commonly adopted approaches (e.g., TCR3d-style geometry).

---

## 3. Design Principles

The angle definitions strictly follow the principles below.

### 3.1 Coordinate-System Independence
All angles are defined using **intrinsic structural vectors and planes**, not external coordinate axes.

### 3.2 Coarse-Grained Robustness
Only **structurally conserved, rigid elements** are used.  
Flexible components such as CDR loops or peptide atoms are intentionally excluded.

### 3.3 Reproducibility and Interpretability
All geometric objects are explicitly defined and can be visualized and validated in standard molecular viewers (e.g., PyMOL, ChimeraX).

---

## 4. Geometric Definitions

### 4.1 TCR Reference Axis

#### Definition

The **TCR axis** is defined as the vector connecting the conserved disulfide bonds of the TCR variable domains:

\[
\vec{v}_{\mathrm{TCR}} =
\mathrm{COM}(\text{V}_\beta\ \text{disulfide}) -
\mathrm{COM}(\text{V}_\alpha\ \text{disulfide})
\]

where:

- each disulfide bond center is computed as the geometric center of the two cysteine sulfur atoms  
- Vα and Vβ domains are identified via chain annotation or numbering schemes (e.g., ANARCI-compatible)

#### Rationale

- Disulfide bond positions in Vα and Vβ are **highly conserved**
- The resulting vector represents the **global orientation of the TCR variable module**
- The definition is insensitive to CDR length and sequence variability

---

### 4.2 pHLA Groove Axis

#### Definition

The **pHLA groove axis** is defined along the peptide-binding groove using the MHC α1 and α2 helices:

\[
\vec{v}_{\mathrm{groove}} =
\mathrm{COM}_{\mathrm{C\text{-}end}} -
\mathrm{COM}_{\mathrm{N\text{-}end}}
\]

where:

- N-terminal and C-terminal regions correspond to predefined residue ranges on the α1/α2 helices  
- the center of mass (COM) is computed using backbone atoms

#### Rationale

- The peptide-binding groove direction is structurally conserved
- This axis defines the **longitudinal reference direction** of pHLA

---

### 4.3 pHLA Groove Plane and Normal

#### Definition

A **groove plane** is fitted using backbone atoms of the α1 and α2 helices via least-squares plane fitting.

The **groove normal vector** is defined as:

\[
\vec{n}_{\mathrm{groove}} = \text{unit normal of groove plane}
\]

#### Rationale

- The groove plane approximates the pHLA binding surface
- The normal vector provides a stable reference for measuring TCR tilt

---

## 5. Angle Definitions

### 5.1 Crossing Angle

\[
\theta_{\mathrm{cross}} =
\angle\left(
\vec{v}_{\mathrm{TCR}},\;
\vec{v}_{\mathrm{groove}}
\right)
\]

**Interpretation**

The crossing angle quantifies the **lateral orientation** of the TCR as it spans across the pHLA peptide-binding groove.

---

### 5.2 Incident Angle

\[
\theta_{\mathrm{inc}} =
\angle\left(
\vec{v}_{\mathrm{TCR}},\;
\vec{n}_{\mathrm{groove}}
\right)
\]

**Interpretation**

The incident angle quantifies the **tilt or inclination** of the TCR relative to the pHLA binding surface.

---

## 6. Algorithmic Workflow

### Step 0. Structural Alignment (Mandatory)

All structures must first be aligned using the **MHC α1/α2 backbone** to ensure a consistent reference frame.

### Step 1. Feature Extraction

- Identify TCR Vα and Vβ disulfide bonds  
- Compute geometric centers of disulfide bonds  
- Extract MHC α1/α2 backbone atoms  
- Construct groove axis and groove plane  

### Step 2. Angle Computation

- Normalize all vectors  
- Compute angles using standard vector dot products  
- Output angles in degrees, within the range \([0^\circ, 180^\circ]\)

### Step 3. Validation (Strongly Recommended)

- Randomly select representative structures  
- Visualize vectors and planes in PyMOL or ChimeraX  
- Confirm geometric consistency with intuitive docking orientation  

---

## 7. Output Specification

For each structure or trajectory frame, the module outputs:

| Field | Description |
|---|---|
| `crossing_angle` | TCR–pHLA groove crossing angle (degrees) |
| `incident_angle` | TCR–pHLA groove incident angle (degrees) |
| `tcr_axis_vector` | Unit vector defining TCR orientation |
| `groove_axis_vector` | Unit vector defining pHLA groove direction |
| `groove_normal_vector` | Unit normal vector of pHLA groove plane |

---

## 8. Limitations and Intended Use

- These angles **do not describe** local binding chemistry or signaling activation  
- They represent a **coarse-grained docking geometry**, not a full reaction coordinate  
- Resulting distributions should be interpreted as **effective, environment-agnostic geometric preferences**, not instantaneous immunological states  

---

## 9. Recommended Applications

- Structural quality control for modeled pHLA–TCR complexes  
- Statistical comparison against crystal or MD-derived reference ensembles  
- Feature construction for machine learning models with explicit geometric priors  
- Identification and classification of docking modes  

---

## 10. Explicit Non-Goals

This module is **not intended** for:

- predicting T cell activation or immunogenicity  
- ranking agonist versus antagonist ligands  
- replacing detailed energetic or allosteric analyses  

---

## 11. Summary

This docking angle definition provides a **minimal, robust, and interpretable geometric representation** of pHLA–TCR binding poses.

It is designed to be **reproducible across large datasets**, **stable under MD sampling**, and **defensible in both physical and biological contexts**.
