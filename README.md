# SHEFF Boundary Value Problem (`sheff-bvp`)

`sheff-bvp` provides the top-level finite element analysis drivers, boundary condition formulations, and time/load-step controllers for solving linear and nonlinear boundary value problems in the **Scalable Hierarchical Enrichment Finite-element Framework (SHEFF)**.

## Capabilities

- **Boundary Value Problem Drivers (`bvp/bvp/`, `bvp/Factory.h`)**:
  - `BoundaryValueProblem` base architecture and centralized state container (`BVPData`) coordinating the mesh, finite element layouts, DoF manager, assembler, global solution/residual vectors, and solver backend.
  - Concrete problem drivers for linear static elasticity (`LinearElasticBVP`), materially nonlinear implicit analyses (`MateriallyNonlinearBVP`), and progressive cohesive/XFEM fracture propagation (`FractureBVP`).
- **Boundary Conditions (`bvp/bc/`)**:
  - Essential and natural boundary conditions with time-dependent load scaling:
    - Prescribed nodal displacements/field values (`DirchletBC.h`).
    - Concentrated nodal forces/fluxes (`NeumannBC.h`).
    - Distributed surface tractions integrated over element faces (`TractionBC.h`).
    - Multi-point periodic boundary conditions for RVE/unit-cell homogenization (`PeriodicBC.h`).
- **Time-Stepping & State Serialization (`TimeDiscretization.h`, `Serialization.h`)**:
  - Manages multi-step loading schedules, adaptive time-step cutbacks on nonlinear or fracture events, and full binary/JSON checkpoint serialization of BVP states and boundary conditions.

## Project Structure

- `include/bvp/`: Public headers organized into `bvp/` (problem classes and `BVPData`), `bc/` (boundary conditions), `TimeDiscretization.h`, and `Factory.h`.
- `source/`: Compiled BVP step-solve loops, Newton/PETSc nonlinear callbacks, and boundary condition constraint/load builders.
- `unit_tests/`: Unit tests covering linear elastic, nonlinear, periodic BC, traction BC, and XFEM fracture BVPs.
- **CMake Package / Target**: `Bvp` / `Bvp::bvp` (static library).

## Dependencies

- **Build System**: [`sheff-build-tools`](https://github.com/afrl/sheff-build-tools) (`SHEFF_TOOLS`), CMake 3.17+, C++17/C++20 compiler.
- **First-Party SHEFF Modules**: [`sheff-assembler`](https://github.com/afrl/sheff-assembler) (`Assembler::assembler`), [`sheff-xfem`](https://github.com/afrl/sheff-xfem) (`Xfem::xfem`).

## Building

When building as part of [SPAWC](https://github.com/afrl/spawc), `sheff-bvp` is configured and built automatically via `make.py`.

To build and test standalone (with `SHEFF_TOOLS` configured and a sparse solver backend enabled):

```bash
python3 make.py --config --build --use-mkl --make-tests --run-tests
```

## Public Release Information

Distribution Statement A. Approved for public release: distribution is unlimited. Case #: AFRL-2026-4513
