# PHOENIX Project

PHOENIX: Process-resolved Hybrid Omniphysics Engine for Nonlinear In-situ X-evolution — Fortran + OpenMP.

## Principles

- **First principles**: Understand the physics before writing code. Every equation and boundary condition must have a clear physical basis.
- **Occam's razor**: The simplest correct solution wins. Do not add complexity unless it is proven necessary.
- **Consistency, simplicity, and maintainability over features**: Clean, readable code is more valuable than new functionality. Do not sacrifice code quality to ship faster. Refactoring to keep things clean is always justified.

## Project Structure

```
PHOENIX/
├── fortran_new/          # Active simulation code
│   ├── main.f90          # Entry point
│   ├── mod_*.f90         # Fortran modules (one per file)
│   ├── compile.sh        # Build script (gfortran + OpenMP)
│   ├── clean.sh          # Remove build artifacts and results
│   ├── run.sh            # Run script: bash run.sh <case_name> [omp_threads]
│   ├── inputfile/
│   │   └── input_param.txt   # All simulation parameters
│   ├── ToolFiles/
│   │   ├── B26.crs           # Active toolpath (hardcoded in mod_toolpath.f90)
│   │   └── toolpath_generator_rectangle.py  # Toolpath generator
│   └── result/               # Output directory (VTK, reports, etc.)
├── legacy/               # Old/reference code (read-only)
└── projects/             # Task tracking (one folder per project)
```

## Build & Run

```bash
cd fortran_new

# Compile (cleans first, then builds all modules + main)
bash compile.sh

# Run simulation
bash run.sh <case_name> [omp_threads] &
# Example: bash run.sh baseline 4 &

# Clean build artifacts
bash clean.sh
```

## Key Files

- **mod_toolpath.f90**: Reads `./ToolFiles/B26.crs` (hardcoded path)
- **input_param.txt**: All simulation parameters (geometry, material, numerics, laser, etc.)
- **mod_sim_state.f90**: Global state arrays and constants
- **mod_laser.f90**: Laser beam positioning from toolpath data

## Toolpath Generator

```bash
cd fortran_new/ToolFiles

# Generate a toolpath (all units in meters and seconds)
python3 toolpath_generator_rectangle.py \
  --start_x 0.0005 --start_y 0.0005 --start_z 0.0006975 \
  --size_x 0.003 --size_y 0.003 \
  --scan_axis x --bidirectional \
  --hatch_spacing 0.0001 --scan_speed 1.23 \
  --turnaround_time 0.0005 --rotation_angle 0 \
  --output B26.crs

# Run all test cases
python3 toolpath_generator_rectangle.py --test
```

Output: `.crs` file (time, x, y, z, laser on/off) + `.png` visualization.
- `.png` is only generated when using `toolpath_generator_rectangle.py` — manually created `.crs` files do not have a `.png`.

## Project Management

- Each project: `projects/YYYYMMDD_PROJECT_NAME/task.md`
- On completion: rename to `YYYYMMDD_PROJECT_NAME-COMPLETEDYYYYMMDD`
- `task.md` contains objective, numbered tasks, and notes

## Species Transport Design Decisions

- Concentration `C` is mass fraction of **primary material**: C=1 = primary (`&material_properties`), C=0 = secondary (`mod_species.f90` constants)
- `hlatnt` (latent heat) remains a **scalar** — assume both materials have similar values
- `dgdt` (dγ/dT, thermal Marangoni coefficient) remains a **scalar** — same assumption
- `dgdc` (dγ/dC, solutal Marangoni coefficient) is a **scalar constant** defined in `mod_species.f90`
- `beta` (thermal expansion coefficient) and `emiss` (emissivity) remain **scalars**, using primary material values
- These simplifications avoid modifying `source_enthalpy` (mod_sour.f90) and keep `bound_uv` changes minimal

## Code Conventions

- All code comments and documentation in English
- One module per `.f90` file, named `mod_<name>.f90`
- `fortran_new/` is the active development folder — all code changes go here
- `legacy/` is read-only reference code — **never modify** files in `legacy/`, including renaming, refactoring, or any text substitutions. Legacy files must remain exactly as they are.
- Do not modify Fortran source code unless the task explicitly requires it
- When adding features, minimize impact on existing modules
- Prefer editing existing files over creating new ones
- **All codebase modifications must be fully documented in the corresponding `task.md`**. Users must be able to track the complete modification progress by reading `task.md` alone. New changes must be written to `task.md` before or alongside code changes.
- **All code changes must be synced to `docs/`** (MkDocs documentation). When adding/modifying modules, input parameters, output fields, or solver behavior, update the corresponding documentation pages. Key files: `docs/architecture/modules.md` (module reference), `docs/architecture/solver.md` (call flow), `docs/input-reference.md` (input parameters), `docs/output.md` (output files), `docs/species/overview.md` (species transport).
