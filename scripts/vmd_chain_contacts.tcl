# VMD script for visualizing dynamic chain-chain interactions
# Usage: vmd -e vmd_chain_contacts.tcl

# Configuration
set topology "md.tpr"
set trajectory "md_processed.xtc"

# Load trajectory
mol new $topology type tpr waitfor all
mol addfile $trajectory type xtc waitfor all

# Set display style
display projection orthographic
display depthcue off
axes location off

# Define chain representations with distinct colors
mol delrep 0 top

# Chain A - Red (protein 1)
mol representation NewCartoon
mol color ColorID 1
mol selection "protein and chain A"
mol material Opaque
mol addrep top

# Chain B - Blue (protein 2)
mol representation NewCartoon
mol color ColorID 0
mol selection "protein and chain B"
mol material Opaque
mol addrep top

# Chain C - Green (protein 3)
mol representation NewCartoon
mol color ColorID 7
mol selection "protein and chain C"
mol material Opaque
mol addrep top

# Chain D - Yellow (peptide/ligand 1)
mol representation NewCartoon
mol color ColorID 4
mol selection "protein and chain D"
mol material Opaque
mol addrep top

# Chain E - Orange (peptide/ligand 2)
mol representation NewCartoon
mol color ColorID 3
mol selection "protein and chain E"
mol material Opaque
mol addrep top

# Add contact visualization between chains
# Contact distance threshold: 3.5 Angstrom
set cutoff 3.5

# Chain A-B contacts (Red-Blue)
mol representation DynamicBonds $cutoff 0.3
mol color ColorID 11
mol selection "protein and ((chain A and within $cutoff of chain B) or (chain B and within $cutoff of chain A))"
mol material Transparent
mol addrep top

# Chain A-C contacts (Red-Green)
mol representation DynamicBonds $cutoff 0.3
mol color ColorID 12
mol selection "protein and ((chain A and within $cutoff of chain C) or (chain C and within $cutoff of chain A))"
mol material Transparent
mol addrep top

# Chain B-C contacts (Blue-Green)
mol representation DynamicBonds $cutoff 0.3
mol color ColorID 13
mol selection "protein and ((chain B and within $cutoff of chain C) or (chain C and within $cutoff of chain B))"
mol material Transparent
mol addrep top

# Chain A-D contacts (Protein-Peptide)
mol representation DynamicBonds $cutoff 0.3
mol color ColorID 10
mol selection "protein and ((chain A and within $cutoff of chain D) or (chain D and within $cutoff of chain A))"
mol material Transparent
mol addrep top

# Chain A-E contacts (Protein-Peptide)
mol representation DynamicBonds $cutoff 0.3
mol color ColorID 9
mol selection "protein and ((chain A and within $cutoff of chain E) or (chain E and within $cutoff of chain A))"
mol material Transparent
mol addrep top

# Optional: Add surface representation for better visualization
# Uncomment if needed
# mol representation Surf
# mol color Chain
# mol selection "protein"
# mol material Transparent
# mol addrep top

# Center view
mol top 0
display resetview

puts "========================================================================"
puts "VMD Chain Contact Visualization Loaded"
puts "========================================================================"
puts ""
puts "Chain color scheme:"
puts "  Chain A: Red"
puts "  Chain B: Blue"
puts "  Chain C: Green"
puts "  Chain D: Yellow"
puts "  Chain E: Orange"
puts ""
puts "Contact representations (distance < ${cutoff} A):"
puts "  A-B contacts: Pink"
puts "  A-C contacts: Cyan"
puts "  B-C contacts: Purple"
puts "  A-D contacts: Lime"
puts "  A-E contacts: Silver"
puts ""
puts "Use VMD timeline to play trajectory and observe dynamic contacts"
puts "========================================================================"
