# CDE — Chemical Discovery Engine

CDE (Chemical Discovery Engine) is a collection of Fortran routines for chemical reaction-path analysis and reaction discovery.

!!! warning "Legacy software"
    This is a legacy version of the reaction-discovery code. An updated Python version is currently in development.

## What can CDE do?

CDE supports several types of chemical-reaction calculations:

- **Double-ended mechanism searching**
- **Generation of initial approximate minimum-energy paths (MEPs)**
- **Linear interpolation and IDPP path generation**
- **Nudged elastic band (NEB) calculations**
- **Reaction-path optimisation**

Single-ended graph-driven sampling is currently unavailable in this legacy version.

## External programs

CDE interfaces with several external programs for energy and force calculations:

- ORCA
- Psi4
- LAMMPS
- DFTB+
- Molpro

## Documentation

<div class="grid cards" markdown>

- :material-rocket-launch: **Getting Started**

    Install CDE and run your first calculation.

    [:octicons-arrow-right-24: Getting started](getting-started.md)

- :material-book-open-variant: **User Guide**

    Input files, PES definitions, path optimisation and configuration.

    [:octicons-arrow-right-24: User guide](user-guide/index.md)

- :material-school: **Tutorials**

    Complete worked examples of CDE calculations.

    [:octicons-arrow-right-24: Tutorials](tutorials/index.md)

- :material-code-braces: **API Reference**

    Documentation generated from the CDE Fortran source.

    [:octicons-arrow-right-24: API reference](api/index.md)

</div>
