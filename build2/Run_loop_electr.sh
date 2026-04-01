#!/bin/bash

# ---- USER SETTINGS ----
NX=20
NY=10

EXEC=./NeutronSource
GDML_DIR=gdml_out          # where Python generated files
OUTPUT_DIR=loopoutputs

mkdir -p $OUTPUT_DIR

# ---- LOOP ----
for ((ix=0; ix<=NX; ix++))
do
  for ((iy=0; iy<=NY; iy++))
  do
    echo "Running x=${ix}, y=${iy}"

    GDML_FILE=${GDML_DIR}/Electronics_x$(printf "%02d" $ix)_y$(printf "%02d" $iy).gdml
    MACRO_FILE=${OUTPUT_DIR}/run_x${ix}_y${iy}.mac

    cat > $MACRO_FILE <<EOF
#
# Macro file for "NeutronSource.cc"
#
#/random/setSeeds 19577794 424238336

/control/verbose 2
/run/verbose 1
#
/testhadr/phys/thermalScattering true

/testhadr/det/readFile ${GDML_FILE}

# To save the energy deposition in the slilcion slabs for dose calculations
# Set-> SetSiliconSlabs 1 and saveSiliconData 1 and saveFluxData 0
# To save the particles emerging from the world volume
# Set-> SetSiliconSlabs 0 and saveSiliconData 0 and saveFluxData 1
#/testhadr/det/SetSiliconSlabs 1
#/stepping/saveSiliconData 1
#/stepping/saveFluxData 0
#
/run/initialize
#
# Set a very high time threshold to allow all decays to happen
#/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year
#
/gps/position 100 150 20 mm
/gps/pos/type Plane
/gps/pos/shape Rectangle
/gps/pos/halfx 15 mm
/gps/pos/halfy 15 mm
/gps/particle neutron
/gps/ene/mono 0.025 eV
/gps/direction 0 0 -1
#
/analysis/setFileName ${OUTPUT_DIR}/Electronics_x${ix}_y${iy}
/analysis/h1/set 4 100  0. 10.  MeV #gammas
/analysis/h1/set 6  60  0. 12.  MeV #neutrons

#
/tracking/verbose 0
#
#/run/printProgress 100000
/run/beamOn 100

EOF

    $EXEC $MACRO_FILE

  done
done