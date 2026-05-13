#!/bin/bash
# Usage: bash compare_apbs.sh protein1.pdb protein2.pdb


PDB1=$(basename "$1" .pdb)
PDB2=$(basename "$2" .pdb)

# Generate PQR files 
pdb2pqr --ff=parse "$1" "${PDB1}.pqr"
pdb2pqr --ff=parse "$2" "${PDB2}.pqr"

# Create APBS input files 

cat > "${PDB1}.in" << EOF
read
  mol pqr ${PDB1}.pqr
end
elec
  mg-auto
  dime 129 129 129
  fglen 80 80 80
  cglen 100 100 100
  mol 1
  lpbe
  bcfl sdh
  pdie 2.0
  sdie 78.54
  srfm smol
  chgm spl2
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
  write pot dx ${PDB1}_pot
end
quit
EOF

cat > "${PDB2}.in" << EOF
read
  mol pqr ${PDB2}.pqr
end
elec
  mg-auto
  dime 129 129 129
  fglen 80 80 80
  cglen 100 100 100
  mol 1
  lpbe
  bcfl sdh
  pdie 2.0
  sdie 78.54
  srfm smol
  chgm spl2
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
  write pot dx ${PDB2}_pot
end
quit
EOF

# Run APBS
echo ">>> Running APBS calculations..."
apbs "${PDB1}.in" > "${PDB1}_apbs.log"
apbs "${PDB2}.in" > "${PDB2}_apbs.log"

