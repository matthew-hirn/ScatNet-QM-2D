Atomic densities are spherical
==============================

The **total** atomic density ``rho(r)`` of an atom has *spherical* symmetry. As a consequence, one only needs to compute that density on a radial grid, starting with the position of the nuclei at r = 0.

How the density were obtained
=============================

The total atomic densities were obtained with **Gaussian09**, using the hybrid *PBE0* functional (called PBE1PBE in Gaussian) using a very large ``cc-pV5Z`` basis set. The use of this basis set basically allows one to produce "exact" atomic densities (by "exact" we mean the "true" *PBE0* density, as if the basis set were complete).

Radial grid used
================

Given the sizes of the shell structures for our atomic radial densities (``4 pi r**2 * rho(r)``), we decided to use to following radial grid::

    1. the grid starts at r = 0 and goes to r = 10 Angstroms
    2. the grid "step" is dr = 0.001 Angstroms
    3. as a consequence, the grid has 10001 points

Important note on using Gaussian09 to compute densities
=======================================================

It is important to note that, because of the fundamental multi-reference character of the many-electron wavefunction of atoms with non-closed shells, the single Slater determinant used in the DFT KS wavefunction is incompatible with the spherical symmetry of the total density. The only exceptions in our set of atoms (H, C, N, O, S, Cl) were H and N (these have exactly half-filled shells which restores the spherical symmetry of the wavefunction).

As a consequence, we used Gaussian09 with all possible symmetries "on" when computing the ground state density. To be exact we used the following command::

    PBE1PBE / cc-pV5Z SCF=(VeryTight,Symm,FSymm,DSymm,IDSymm,IntRep) SCF=QC Density=SCF

Empirical observation showed that when generating densities in x, y and z directions, two of those densities were similar, while the last was qualitatively similar in appearance (same number of shell structures in the radial density profile) but quantitatively slightly different.

In order to generate a density profile that integrates to the proper number of electrons, we averaged the normalized density profiles from the 3 directions. The algorithm was:

    1. Compute the integrated charges in x, y and z directions from the x, y and z densities
    2. Normalize each density profile so as to recover the proper number of electrons (from the "true" value of the number of electrons and the values found in 1.)
    3. Compute the total density by averaging the normalized densities from 2.

Conservation of the number of electrons
=======================================

Given the definition of the total atomic density, we made a check that the total number of electron was conserved. More specifically, we found that ``4 pi r**2 * rho(r)`` integrated from r = 0 to r = 10 gave::

    Q = 1 for Hydrogen (H)
    Q = 6 for Carbon (C)
    Q = 7 for Nitrogen (N)
    Q = 8 for Oxygen (O)
    Q = 16 for Sulfur (S)
    Q = 17 for Chlorine (Cl)

Format of the mat files
=======================

Each mat file contains 3 records. These records are:

    1. an array containing the values of the density ``rho(r)``
    2. an array containing the radial grid used to compute the density in 1.
    3. a string giving the units for the radial grid (here *bohr*

Please do **not** confuse the density ``rho(r)`` from the *radial* density ``4 * pi * r**2 * rho(r)``. It is the latter one that needs to be integrated from 0 to infinity to recover the total number of electrons. The unit of the density is (# of electrons per bohr cubed).

Generation of 2D densities
==========================

In order to use the exact atomic densities in a purely 2D context (i.e. within the 2D molecular plane), we integrated the three dimensional density found above in the z-direction. This allowed us to generate purely 2D densities. These 2D densities -- because of the spherical symmetry of the three dimensional density -- have cylindrical symmetry. More specifically the 2D density is defined from the 3D density by the following integration::

    rho_2D(r) = 2. * Integral( rho_3D((r**2 + z**2)**(1/2)) , z from 0 to infinity)

As a consequence, the 2D density verifies the following identity, which expresses the conservation of the number of electrons in the 2D case::

    N = # of electrons = Integral( 2. * pi * r * rho_2D(r), r from 0 to infinity)

The 2D densities were produced for a square window of side length 5 Angstroms, with a step of 0.01 Angstrom. The integrate charge was consistent with the theoretical value for the respective atom to within 1e-2 units. The atom position is at the center of each 501x501 2D density array.

Generation of 2D densities that keep the spherical symmetry
===========================================================

Upon a discussion with Stephane, it turned out that the following 2D density might be more appropriate than the above defined 2D density.

We define that density by "condensating the charge mass" from the sphere of radius *r* to the 2D circle of radius *r*, by defining ``rho_2D(r)`` to be such that::

    dQ = 4. * pi * r**2 * rho_3D(r) * dr = 2. * pi * r * rho_2D(r) * dr

which leads to::

    rho_2D(r) = 2. * r * rho_3D(r)

The 2D densities were produced with the same parameters as above. They consists in 2D arrays of shape 501x501 corresponding to squares of 5 Angstroms of side length with a step of 0.01 Angstrom.
