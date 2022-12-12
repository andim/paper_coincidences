# Measures of epitope binding degeneracy from T cell receptor repertoires

This repository contains the source code associated with the manuscript

Mayer, Callan: [Measures of epitope binding degeneracy from T cell receptor repertoires](https://doi.org/10.1101/2022.07.25.501373), bioRxiv preprint, 2022

It allows reproduction of the statistical analyses and numerical results reported in the manuscript.

For a number of figures final assembly and cosmetic changes were done in Inkscape as a postprocessing step. In these cases the figures will not be reproduced precisely. To help reuse the final edited figures are provided in png/svg format.

## Installation requirements

The software is written in Python, and was run on Python version 3.6. The code relies on [Pyrepseq](https://github.com/andim/pyrepseq) (version 1.0), a python package for the analysis of immune repertoire sequencing data that we have released to accompany this paper. Other packages used include:

- [numpy](http://github.com/numpy/numpy/)
- [scipy](https://github.com/scipy/scipy)
- [pandas](http://github.com/pydata/pandas)
- [matplotlib](http://github.com/matplotlib/matplotlib)
- [seaborn](http://github.com/mwaskom/seaborn)
- [networkx](https://github.com/networkx/networkx)

All can be installed using:

`pip install -r requirements.txt`

## Running the code

Data download and preprocessing is handled by the [Snakemake](https://snakemake.github.io/) workflow manager, with the help of scripts located within the `scripts` directory. Data visualization is done within [Jupyter](https://jupyter.org/) notebooks provided for each of the figures.

## Contact

If you run into any difficulties running the code, please contact us at `andimscience@gmail.com`.

## License

The source code is freely available under an MIT license. The plots are licensed under a Creative Commons attributions license (CC-BY).
