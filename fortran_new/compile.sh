#!/bin/bash
# PHOENIX Build Script — compile only
set -e

# Clean build artifacts only (preserve results)
find . -maxdepth 1 -name "*.o" -exec rm -f '{}' \;
find . -maxdepth 1 -name "*.mod" -exec rm -f '{}' \;
find . -maxdepth 1 -name "cluster_main*" -exec rm -f '{}' \;

echo "Compiling PHOENIX..."

gfortran -fopenmp -O3 -march=native -c \
    mod_precision.f90 \
    mod_const.f90 \
    mod_cfd_utils.f90 \
    mod_param.f90 \
    mod_geom.f90 \
    mod_field_data.f90 \
    mod_coeff_data.f90 \
    mod_sim_state.f90 \
    mod_init.f90 \
    mod_laser.f90 \
    mod_dimen.f90 \
    mod_adaptive_mesh.f90 \
    mod_resid.f90 \
    species_solver/mod_species.f90 \
    mod_prop.f90 \
    mod_bound.f90 \
    mod_discret.f90 \
    mod_entot.f90 \
    mod_predict.f90 \
    mod_sour.f90 \
    mod_flux.f90 \
    mod_revise.f90 \
    mod_solve.f90 \
    mod_print.f90 \
    mod_converge.f90 \
    mod_toolpath.f90 \
    mod_timing.f90 \
    mod_defect.f90 \
    mechanical_solver/mod_mech_material.f90 \
    mechanical_solver/mod_mechanical.f90 \
    mechanical_solver/mod_mech_io.f90 \
    main.f90

gfortran -fopenmp -O3 -march=native *.o -o cluster_main

echo "Build complete: cluster_main"
echo ""
echo "Usage: bash run.sh <case_name> [thermal_threads] [mech_threads] &"
echo "  Example: bash run.sh baseline 4 &           # serial mechanical"
echo "  Example: bash run.sh baseline 10 10 &        # parallel mechanical"
echo "  Stop all:  kill \$(pgrep -f cluster_main)"
echo "  Stop one:  ps aux | grep cluster_main  then  kill <PID>"
