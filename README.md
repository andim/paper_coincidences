# Measures of epitope binding degeneracy from T cell receptor repertoires

This repository contains the source code associated with the manuscript

Mayer, Callan: [Measures of epitope binding degeneracy from T cell receptor repertoires](https://doi.org/10.1101/2022.07.25.501373), bioRxiv preprint, 2022

It allows reproduction of the statistical analyses and numerical results reported in the manuscript.

For a number of figures final assembly and cosmetic changes were done in Inkscape as a postprocessing step. In these cases the figures will not be reproduced precisely. To help reuse the final edited figures are provided in png/svg format.

## Installation requirements

The code makes heavy use of [Pyrepseq](https://github.com/andim/pyrepseq), a python package for the analysis of immune repertoire sequencing data that we have released to accompany this paper.

`pip install pyrepseq==1.0`


Python 3.6+.

A number of standard scientific python packages are needed for the statistical analyses and visualizations. An easy way to install many of these is to install a Python distribution such as [Anaconda](https://www.continuum.io/downloads).

- [numpy](http://github.com/numpy/numpy/)
- [scipy](https://github.com/scipy/scipy)
- [pandas](http://github.com/pydata/pandas)
- [matplotlib](http://github.com/matplotlib/matplotlib)
- [seaborn](http://github.com/mwaskom/seaborn)
- [networkx](https://github.com/networkx/networkx)



## Contact

If you run into any difficulties running the code, please contact us at `andimscience@gmail.com`.

## License

The source code is freely available under an MIT license. The plots are licensed under a Creative Commons attributions license (CC-BY).
