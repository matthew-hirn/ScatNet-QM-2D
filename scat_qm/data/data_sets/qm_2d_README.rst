Version Log
===========

June 28th, 2014::

    Add electron density within the molecular plane. Each density file in MATLAB format consists of a dictionnary with 2 keys.
    The first keys is called ``shape`` and corresponds to the shape of the 2D density array. The second key is called ``density``
    and consists of the density array.
    The dataset now contains a new attribute called ``D`` (for density). This array contains the names of the density mat file in directory ``densities`` that correspond to the given molecule.
    Regenerated new splits and changed the order of molecules in the dataset to be consistent and reproducible.

June 8th, 2014::

    revision of the atomization energies. There was a mistake in the way the total energies of the atomic species (here only
    H, C, N, O, S, Cl) were computed (the multiplicities of C, N, O and S were wrong). As a consequence, the atomization energies
    were wrong as well. The dataset targets 'T' have been updated consequently. A plot showing the change in atomization energies
    has been added for visual inspection. Regenerated the splits atomization profile PNG image as well.

June 6th, 2014::

    first draft of the full dataset in MATLAB format

How was the dataset assembled?
==============================

The mother database from which ``QM2D`` has been created is ``GDB13``. The following steps will describe the methodology that
has been followed in assembling the dataset.

1. From the **SMILES** characters of ``GDB13`` for all molecules up to and including 9 "heavy" atoms, a rough first initial
   guess for the cartesian coordinates were obtained using the Open Source code [OpenBabel]_. More specifically, the command
   line used to extract an ``smi`` file called ``myfile.smi`` was::

   obabel -ismi myfile.smi -oxyz myfile.xyz --gen3d -c

2. From the XYZ formatted files representing the 319398 molecules found, an SVD operation was performed on each molecule's
   coordinate matrix in order to extract the smallest of the 3 singular values. A strictly 2D molecule would have a value of
   0 for its smallest singular value, indicating that all the atoms belong to a 2D vector space (here a plane). We chose to
   screen all the molecules and only keep the ones with a smallest singular values smaller than 0.3.

3. The geometry of each of the 6002 molecules found at step 2 were fully optimized. First at a semi-empirical level using
   ``PM6`` and then a ``PBE0/6-31+G*`` level of theory. The reason for using the quite reduced ``6-31+G*`` basis set was
   motivated by its excellent accuracy and computational efficiency as outline in [Ref1]_. 5980 out of the initial 6002
   could be fully relaxed. Gaussian09 was used with options ``Opt(Tight,CalcFC,maxstep=100)``.

4. Out of the optimized geometries, 4371 molecules were strictly 2D. The z axis was chosen consistently as the axis normal to
   the molecular plane. For these 4371 molecules, the Standard Deviation of the atoms z coordinates are less than 0.0005 angstroms.

5. The next step consisted in selecting only the **stable** molecules by performing a Frequency calculation and making sure
   that no negative frequencies were found in the set of frequencies for which the translational and rotational degrees of freedom
   were taken out. We ended up with a final number of 4357 molecules.

6. The last step included a single point energy calculation at the ``PBE0/aug-cc-pVTZ`` level of theory for the relaxed geometry.
   The sum of the atomic energies of all the atoms in the molecule (computed at the same level of theory) was subtracted from the
   total energy to produce the molecule's atomization energy.


Format of the dataset
=====================

We kept the same format as the original ``QM7`` dataset for generating the Python Pickle file and Matlab mat file.


How to create "homogeneous" splits for the cross validation?
============================================================

As was explained in one of QM7 original paper [Ref2]_, achieving cross validation splits for which the distribution of organic
functional groups is homogeneous is critical for improving the training of Machine Learning models, and especially their transferability.
In order to address this issue, we created our splits by first exhaustively describing all the functional groups (out of more than
200 distinct organic functional groups) contained in each molecule of the dataset. Once the raw statistics was obtained we chose to
distribute molecules with the rarest functional groups homogeneously over the K splits. Then came the turn of all the molecules from
the rest of the dataset with the second rarest functional group to be distributed homogeneously over the K splits. We continued the
process until no more molecule was available.

5 splits were created, guided by::

    1. a goal of homogeneous representation of the different functional groups in the dataset
    2. a distribution of atomization energies as equivalent as possible between the splits

Some Statistics about the dataset in terms of functional groups
===============================================================

The dataset contains 60 distinct functional groups. The functional groups were found by analyzing the structure of each molecule using the ``checkmol`` utility from **Norbert Haider**.

+-----------------------------------------------------------+------------------------------------+
|                                   Functional Group        |       # molecules with that group  |
+-----------------------------------------------------------+------------------------------------+
|                                        1,2-diphenol       |              13                    |
+-----------------------------------------------------------+------------------------------------+
|                                        acyl cyanide       |              16                    |
+-----------------------------------------------------------+------------------------------------+
|                                            aldehyde       |             560                    |
+-----------------------------------------------------------+------------------------------------+
|                                              alkene       |            1527                    |
+-----------------------------------------------------------+------------------------------------+
|                                      alkyl chloride       |               2                    |
+-----------------------------------------------------------+------------------------------------+
|                                              alkyne       |            1169                    |
+-----------------------------------------------------------+------------------------------------+
|                                              aminal       |              37                    |
+-----------------------------------------------------------+------------------------------------+
|                                   aromatic compound       |            3208                    |
+-----------------------------------------------------------+------------------------------------+
|                                       aryl chloride       |             221                    |
+-----------------------------------------------------------+------------------------------------+
|                                        azo compound       |               5                    |
+-----------------------------------------------------------+------------------------------------+
|                                       carbamic acid       |               4                    |
+-----------------------------------------------------------+------------------------------------+
|                      carbamic acid ester (urethane)       |               2                    |
+-----------------------------------------------------------+------------------------------------+
|                                        carbonitrile       |             717                    |
+-----------------------------------------------------------+------------------------------------+
|                                     carboxylic acid       |              46                    |
+-----------------------------------------------------------+------------------------------------+
|                             carboxylic acid amidine       |              32                    |
+-----------------------------------------------------------+------------------------------------+
|                               carboxylic acid ester       |              68                    |
+-----------------------------------------------------------+------------------------------------+
|                           carboxylic acid hydrazide       |               4                    |
+-----------------------------------------------------------+------------------------------------+
|                carboxylic acid imide, N-substituted       |               6                    |
+-----------------------------------------------------------+------------------------------------+
|              carboxylic acid imide, N-unsubstituted       |              17                    |
+-----------------------------------------------------------+------------------------------------+
|                                             enamine       |             224                    |
+-----------------------------------------------------------+------------------------------------+
|                                                enol       |              26                    |
+-----------------------------------------------------------+------------------------------------+
|                                          enol ether       |              46                    |
+-----------------------------------------------------------+------------------------------------+
|                                           guanidine       |              11                    |
+-----------------------------------------------------------+------------------------------------+
|                                  halogen derivative       |               2                    |
+-----------------------------------------------------------+------------------------------------+
|                                          hemiacetal       |               1                    |
+-----------------------------------------------------------+------------------------------------+
|                                          hemiaminal       |               2                    |
+-----------------------------------------------------------+------------------------------------+
|                               heterocyclic compound       |            3608                    |
+-----------------------------------------------------------+------------------------------------+
|                                hydrazine derivative       |             286                    |
+-----------------------------------------------------------+------------------------------------+
|                                           hydrazone       |              45                    |
+-----------------------------------------------------------+------------------------------------+
|                                     hydroxamic acid       |               4                    |
+-----------------------------------------------------------+------------------------------------+
|                                       hydroxylamine       |              65                    |
+-----------------------------------------------------------+------------------------------------+
|                                         imido ester       |              15                    |
+-----------------------------------------------------------+------------------------------------+
|                                               imine       |              91                    |
+-----------------------------------------------------------+------------------------------------+
|                                     imino(het)arene       |              17                    |
+-----------------------------------------------------------+------------------------------------+
|                                             isourea       |              15                    |
+-----------------------------------------------------------+------------------------------------+
|                         ketene acetal or derivative       |              90                    |
+-----------------------------------------------------------+------------------------------------+
|                                              ketone       |             289                    |
+-----------------------------------------------------------+------------------------------------+
|                                              lactam       |              74                    |
+-----------------------------------------------------------+------------------------------------+
|                                             lactone       |              21                    |
+-----------------------------------------------------------+------------------------------------+
|                                      nitro compound       |             132                    |
+-----------------------------------------------------------+------------------------------------+
|                                    nitroso compound       |               5                    |
+-----------------------------------------------------------+------------------------------------+
|                     orthocarboxylic acid derivative       |              34                    |
+-----------------------------------------------------------+------------------------------------+
|                                               oxime       |             125                    |
+-----------------------------------------------------------+------------------------------------+
|                                         oxime ether       |              16                    |
+-----------------------------------------------------------+------------------------------------+
|                                       oxo(het)arene       |             447                    |
+-----------------------------------------------------------+------------------------------------+
|                           phenol or hydroxyhetarene       |             456                    |
+-----------------------------------------------------------+------------------------------------+
|                                     primary alcohol       |              46                    |
+-----------------------------------------------------------+------------------------------------+
|                primary aliphatic amine (alkylamine)       |               9                    |
+-----------------------------------------------------------+------------------------------------+
|                                       primary amine       |             162                    |
+-----------------------------------------------------------+------------------------------------+
|                              primary aromatic amine       |             153                    |
+-----------------------------------------------------------+------------------------------------+
|                       primary carboxylic acid amide       |              53                    |
+-----------------------------------------------------------+------------------------------------+
|                                   secondary alcohol       |               2                    |
+-----------------------------------------------------------+------------------------------------+
|                     secondary carboxylic acid amide       |             130                    |
+-----------------------------------------------------------+------------------------------------+
|                      tertiary carboxylic acid amide       |              19                    |
+-----------------------------------------------------------+------------------------------------+
|                                        thioaldehyde       |              14                    |
+-----------------------------------------------------------+------------------------------------+
|                                   thiocarbamic acid       |               4                    |
+-----------------------------------------------------------+------------------------------------+
|                                          thioketone       |               3                    |
+-----------------------------------------------------------+------------------------------------+
|                                            thiourea       |               7                    |
+-----------------------------------------------------------+------------------------------------+
|                                    thioxo(het)arene       |              26                    |
+-----------------------------------------------------------+------------------------------------+
|                                                urea       |              24                    |
+-----------------------------------------------------------+------------------------------------+


.. [OpenBabel] This useful utility can be found at: http://openbabel.org/wiki/Main_Page
.. [Ref1] Riley KE et al. JCTC, 3(2), pp 407-433, 2007
.. [Ref2] Assessment and Validation of Machine Learning Methods for Predicting Molecular Atomization Energies, JCTC 9, 3404-3419 (2013)
