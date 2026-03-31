# PHOENIX

**PHOENIX: Process-resolved Hybrid Omniphysics Engine for Nonlinear In-situ X-evolution**

PHOENIX is a hybrid physics-based solver for simulating thermo-fluid processes in laser-based additive manufacturing.

## Key Features

- **Multi-physics**: Coupled heat transfer, fluid flow (Navier-Stokes), phase change (melting/solidification), and species transport
- **Marangoni convection**: Surface-tension-driven flow from thermal and solutal gradients
- **Dissimilar metal mixing**: Species transport with composition-dependent material properties (two-way coupling)
- **Adaptive mesh**: Movable structured mesh in X-Y that follows the laser/melt pool, providing fine resolution where needed
- **Defect prediction**: Post-simulation analysis of lack-of-fusion and keyhole porosity from max temperature field
- **Melt pool tracking**: Time-series logging of melt pool length, depth, width, volume, and D/W aspect ratio with auto-generated plots
- **OpenMP parallelization**: Shared-memory parallelism for all hot loops (TDMA, discretization, source terms)
- **VTK output**: Binary structured-grid VTK files for ParaView visualization

## Physics

The solver models:

- **Heat transfer**: Conduction with temperature-dependent properties, volumetric laser heat source (Gaussian distribution)
- **Fluid flow**: Incompressible Navier-Stokes with SIMPLE algorithm on staggered grid
- **Phase change**: Enthalpy-based method with Darcy-type resistance in mushy zone
- **Surface effects**: Thermal Marangoni (dg/dT) and solutal Marangoni (dg/dC) stress on free surface
- **Buoyancy**: Boussinesq approximation for natural convection
- **Species transport**: Convection-diffusion of concentration field with molecular diffusivity
- **Powder layer**: Distinct thermal properties for unconsolidated powder
- **Defect prediction**: Post-simulation detection of lack-of-fusion ($T_{max} < T_s$, incomplete melting) and keyhole porosity ($T_{max} > T_b$, excessive vaporization) from peak temperature history within the build layer

## Quick Start

```bash
cd fortran_new
bash compile.sh              # Build
bash run.sh mycase 4 &       # Run with 4 OpenMP threads
```

Results are written to `fortran_new/result/mycase/`.

## Project Structure

```
PHOENIX/
├── fortran_new/          # Active simulation code
│   ├── main.f90          # Entry point
│   ├── mod_*.f90         # Fortran modules (one per file)
│   ├── compile.sh        # Build script
│   ├── run.sh            # Run script
│   ├── clean.sh          # Clean build artifacts
│   ├── inputfile/        # Simulation parameters
│   └── ToolFiles/        # Toolpath files (.crs)
├── legacy/               # Reference code (read-only)
├── projects/             # Task tracking
└── docs/                 # This documentation
```
