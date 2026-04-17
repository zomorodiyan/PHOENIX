#!/bin/bash
# PHOENIX Build Script — compile only
set -e

# Clean build artifacts only (preserve results)
find . -maxdepth 1 -name "*.o" -exec rm -f '{}' \;
find . -maxdepth 1 -name "*.mod" -exec rm -f '{}' \;
find . -maxdepth 1 -name "cluster_main*" -exec rm -f '{}' \;

echo "Compiling PHOENIX..."

gfortran -fopenmp -O3 -march=native -c \
    thermal_fluid_solver/mod_precision.f90 \
    thermal_fluid_solver/mod_const.f90 \
    thermal_fluid_solver/mod_cfd_utils.f90 \
    thermal_fluid_solver/mod_param.f90 \
    thermal_fluid_solver/mod_geom.f90 \
    thermal_fluid_solver/mod_field_data.f90 \
    thermal_fluid_solver/mod_coeff_data.f90 \
    thermal_fluid_solver/mod_sim_state.f90 \
    thermal_fluid_solver/mod_init.f90 \
    thermal_fluid_solver/mod_laser.f90 \
    thermal_fluid_solver/mod_dimen.f90 \
    thermal_fluid_solver/mod_adaptive_mesh.f90 \
    thermal_fluid_solver/mod_resid.f90 \
    species_solver/mod_species.f90 \
    thermal_fluid_solver/mod_prop.f90 \
    thermal_fluid_solver/mod_bound.f90 \
    thermal_fluid_solver/mod_discret.f90 \
    thermal_fluid_solver/mod_entot.f90 \
    thermal_fluid_solver/mod_predict.f90 \
    thermal_fluid_solver/mod_sour.f90 \
    thermal_fluid_solver/mod_flux.f90 \
    thermal_fluid_solver/mod_revise.f90 \
    thermal_fluid_solver/mod_solve.f90 \
    thermal_fluid_solver/mod_print.f90 \
    thermal_fluid_solver/mod_converge.f90 \
    thermal_fluid_solver/mod_toolpath.f90 \
    thermal_fluid_solver/mod_timing.f90 \
    thermal_fluid_solver/mod_defect.f90 \
    mechanical_solver/mod_mech_material.f90 \
    mechanical_solver/mod_mechanical.f90 \
    mechanical_solver/mod_mech_io.f90 \
    thermal_fluid_solver/main.f90

gfortran -fopenmp -O3 -march=native *.o -o cluster_main

echo "Build complete: cluster_main"
echo ""
echo "Usage: bash run.sh <case_name> [thermal_threads] [mech_threads] &"
echo "  Example: bash run.sh baseline 4 &           # serial mechanical"
echo "  Example: bash run.sh baseline 10 10 &        # parallel mechanical"
echo "  Stop all:  kill \$(pgrep -f cluster_main)"
echo "  Stop one:  ps aux | grep cluster_main  then  kill <PID>"
