# VMD Chain Contact Visualization Guide

## Overview

This guide shows how to visualize and analyze dynamic chain-chain interactions in multi-chain complexes using VMD and quantitative analysis.

## System Requirements

- VMD (Visual Molecular Dynamics)
- GROMACS
- Python with matplotlib

## Two-Step Workflow

### Step 1: Visual Exploration (VMD)

Use VMD to interactively explore chain contacts during trajectory playback.

### Step 2: Quantitative Analysis (Python)

Calculate contact distances, frequencies, and generate publication-quality plots.

---

## Step 1: VMD Visualization

### Basic Usage

```bash
# Navigate to your trajectory directory
cd data/7n1e_1/prod/pbc_comparison_nojump/

# Launch VMD with contact visualization script
vmd -e ../../../../scripts/vmd_chain_contacts.tcl
```

### What the Script Does

1. **Loads trajectory** with PBC-corrected coordinates
2. **Colors chains** for easy identification:
   - Chain A: Red (protein 1)
   - Chain B: Blue (protein 2)
   - Chain C: Green (protein 3)
   - Chain D: Yellow (peptide 1)
   - Chain E: Orange (peptide 2)

3. **Displays dynamic contacts** (distance < 3.5 Å):
   - A-B contacts: Pink bonds
   - A-C contacts: Cyan bonds
   - B-C contacts: Purple bonds
   - A-D contacts: Lime bonds
   - A-E contacts: Silver bonds

### Interactive Controls

In VMD:

1. **Play trajectory**: Use Graphics → Representations → Trajectory
2. **Rotate view**: Left mouse drag
3. **Zoom**: Scroll wheel or right mouse drag
4. **Toggle representations**: Check/uncheck in Graphics window

### Customization

Edit `vmd_chain_contacts.tcl` to adjust:

```tcl
# Change contact distance cutoff
set cutoff 4.0  # Increase to 4.0 Angstrom

# Add new chain pair
mol representation DynamicBonds $cutoff 0.3
mol color ColorID 14
mol selection "protein and ((chain B and within $cutoff of chain D) or (chain D and within $cutoff of chain B))"
mol material Transparent
mol addrep top
```

### Tips for Better Visualization

1. **Simplify view for large systems**:
   ```tcl
   # Show only backbone
   mol representation NewCartoon
   mol selection "protein and backbone"
   ```

2. **Add transparency to see buried contacts**:
   ```tcl
   mol material Transparent
   ```

3. **Show specific residues in contact**:
   ```tcl
   mol representation Licorice
   mol selection "protein and ((chain A and resid 10 to 20) or (chain B and resid 50 to 60))"
   ```

---

## Step 2: Quantitative Analysis

### Basic Usage

```bash
# Run contact analysis script
python scripts/analyze_chain_contacts.py
```

### What the Script Calculates

1. **Minimum distance** between each chain pair at every frame
2. **Contact frequency**: % of frames where distance < 3.5 Å
3. **Distance statistics**: min, mean, max for each pair
4. **Time evolution**: How contacts change over trajectory

### Output Files

```
data/7n1e_1/prod/chain_contact_analysis/
├── mindist_A-B.xvg                    # Raw distance data
├── mindist_A-C.xvg
├── mindist_B-C.xvg
├── mindist_A-D_Protein-Peptide.xvg
├── mindist_A-E_Protein-Peptide.xvg
├── chain_distances_all.png            # All pairs overlaid
├── contact_frequency.png              # Bar chart of contact %
└── chain_distances_individual.png     # Grid of individual plots
```

### Interpreting Results

#### Terminal Output Example

```
A-B:
  Contact percentage: 45.2%
  Average distance: 0.42 nm
  Min distance: 0.28 nm

A-D (Protein-Peptide):
  Contact percentage: 89.7%
  Average distance: 0.31 nm
  Min distance: 0.23 nm
```

**Interpretation:**
- **High contact % (>70%)**: Stable interaction throughout simulation
- **Medium contact % (30-70%)**: Transient/dynamic interaction
- **Low contact % (<30%)**: Rare or no direct interaction

#### Plot 1: All Chain Distances

![Example](https://via.placeholder.com/800x400?text=All+Chain+Distances)

- **Stable interactions**: Flat lines below cutoff
- **Dynamic interactions**: Fluctuating around cutoff
- **Dissociation events**: Sharp increases above cutoff

#### Plot 2: Contact Frequency Bar Chart

![Example](https://via.placeholder.com/800x400?text=Contact+Frequency)

- Quickly identify which chain pairs interact most
- Compare relative interaction strengths
- Identify dominant binding interfaces

#### Plot 3: Individual Chain Pair Plots

![Example](https://via.placeholder.com/800x400?text=Individual+Plots)

- Green shading: Frames in contact
- Identify specific time windows of interaction
- Correlate with RMSD/other properties

---

## Advanced Usage

### Custom Chain Selections

Edit `analyze_chain_contacts.py` to analyze specific regions:

```python
chain_pairs = [
    ("chain A and resid 10 to 50", "chain B and resid 60 to 100", "A(10-50)-B(60-100)"),
    ("chain A and name CA", "chain D", "A_backbone-D"),
]
```

### Adjusting Contact Cutoff

```python
# In analyze_chain_contacts.py
CONTACT_CUTOFF = 0.5  # Change to 5.0 Angstrom
```

### Residue-Level Contact Analysis

For detailed per-residue contacts, use GROMACS mindist with `-or` option:

```bash
gmx mindist \
    -f md_processed.xtc \
    -s md.tpr \
    -od distance.xvg \
    -or residue_contacts.xvg \
    -d 0.35
```

### Integration with Other Analyses

Combine contact analysis with:

1. **RMSD**: Correlate contact changes with structural changes
2. **RMSF**: Identify flexible regions affecting contacts
3. **Hydrogen bonds**: Detailed interaction characterization
4. **Distance analysis**: Specific atom pair distances

---

## Troubleshooting

### Issue: VMD script fails to load

**Solution:** Check file paths in script:

```tcl
# Edit vmd_chain_contacts.tcl
set topology "md.tpr"  # Use absolute path if needed
set trajectory "md_processed.xtc"
```

### Issue: No contacts detected

**Possible causes:**
1. Contact cutoff too small → Increase cutoff
2. Chains don't interact → Verify chain IDs with `gmx check`
3. PBC artifacts → Use nojump-processed trajectory

### Issue: Chain IDs not recognized

**Solution:** Check chain assignments in your system:

```bash
# List chains in structure
gmx check -f md.gro
```

Or in VMD:
```tcl
# In VMD TkConsole
set sel [atomselect top "protein"]
$sel get chain
```

### Issue: Analysis script takes too long

**Solution:** Reduce trajectory size:

1. Use already downsampled trajectory (dt=100 ps)
2. Analyze specific time window:
   ```bash
   gmx trjconv -f md_processed.xtc -s md.tpr -o short.xtc -b 50000 -e 60000
   ```

---

## Example Workflow

### Complete Analysis Pipeline

```bash
# 1. Process trajectory with nojump (if not done)
python scripts/pbc_process.py -f md.xtc -s md.tpr -o processed/ --nojump

# 2. Visual exploration
cd processed/
vmd -e ../scripts/vmd_chain_contacts.tcl

# 3. Quantitative analysis
python scripts/analyze_chain_contacts.py

# 4. Review results
ls chain_contact_analysis/
```

### Combining with RMSD Analysis

```python
# After running both analyses
import matplotlib.pyplot as plt

# Load RMSD data
time_rmsd, rmsd = parse_rmsd("rmsd_backbone.xvg")

# Load contact data
time_contact, distance_AB = parse_mindist("chain_contact_analysis/mindist_A-B.xvg")

# Plot correlation
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

ax1.plot(time_rmsd, rmsd, label='RMSD')
ax1.set_ylabel('RMSD (nm)')
ax1.legend()

ax2.plot(time_contact, distance_AB, label='A-B distance', color='red')
ax2.axhline(0.35, linestyle='--', color='black', label='Contact cutoff')
ax2.set_ylabel('Distance (nm)')
ax2.set_xlabel('Time (ns)')
ax2.legend()

plt.tight_layout()
plt.savefig('rmsd_vs_contacts.png', dpi=300)
```

---

## Biological Insights

### What to Look For

1. **Stable binding interfaces**:
   - High contact % (>80%)
   - Low distance fluctuation
   - Example: Essential protein-protein interactions

2. **Transient interactions**:
   - Medium contact % (30-70%)
   - Periodic binding/unbinding
   - Example: Regulatory interactions, allosteric effects

3. **Dissociation events**:
   - Contact % drops over time
   - Increasing distances
   - Example: Ligand unbinding, complex dissociation

4. **Correlation with structural changes**:
   - Contact loss coincides with RMSD increase
   - Suggests conformational change affects binding

---

## Performance Notes

**Typical analysis time (100 ns, 5 chains, 9 pairs):**
- VMD visualization: Real-time playback
- Contact distance calculation: ~2-5 minutes per pair
- Plot generation: ~10 seconds
- **Total: ~20-45 minutes for complete analysis**

---

## References

- [VMD User Guide](https://www.ks.uiuc.edu/Research/vmd/current/ug/)
- [GROMACS mindist](http://manual.gromacs.org/current/onlinehelp/gmx-mindist.html)
- [Protein-protein interaction analysis](https://doi.org/10.1016/j.sbi.2011.01.007)

---

**AfterMD - Comprehensive MD Analysis Toolkit**
