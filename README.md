# paper-Glerum-Sensitivity-surface-deformation-mantle-dynamics

This repository belongs to the paper

    Sensitivity of horizontal surface deformation to mantle dynamics:
    3D instantaneous dynamics modeling of the eastern Mediterranean
    by
    A. Glerum, W. Spakman, D. J. J. van Hinsbergen, C. Thieulot and C. Pranger.

and contains the input files, postprocessing scripts and source code to reproduce the reported results.

Contents
--------
``prm/``

This directory contains the input files to the ASPECT code used for the numerical modeling
(version 1.5.0-pre, commit e1d63fb4dcda9ec98e34 b6ba5451236e77bc59f9).
The file prm_to_model_names.txt links the prm file names to the names used in the paper.

``code/`` 

This directory contains custom ASPECT plugins that were used as shared libraries. The CMakeLists.txt
file can be used to compile them as shared libraries to be called from the prm files. Note that the
prm files now contain paths to the shared libraries that will not exist on your system and need to be
replaced.

``python_visu_scripts/``

This directory contains the ParaView python scripts that were used to process the data and produce
most of the figures in the paper. 

``python_stat_scripts/``

This directory contains the python v2.6.6 scripts that were used to statistically analyse the 
reference model result.

``ASPECT_install_configuration.txt``

This file details the ASPECT and deal.II installation configurations as used for the paper.
