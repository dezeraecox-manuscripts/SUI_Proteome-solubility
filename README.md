# XS_Proteome-solubility
This repo contains the code and example data associated with the proteome solubility project, led by Dr Xiaojing Sui. This manuscript has been submitted for publication.

Manuscript preprint: https://www.biorxiv.org/content/biorxiv/early/2019/07/26/692103.full.pdf

Final version: _upon acceptance_

## Prerequisites

This analysis assumes a standard installation of Python 3 (=> 3.6), as well as [FiJi](https://fiji.sc/) (v 2.0.0) <sup>[1](footnote_1)</sup> and [CellProfiler](https://cellprofiler.org/) (v 3.1.9) <sup>[2](footnote_2)</sup>. For specific package requirements, see the `environment.yml` file.

## Test data

[Example raw images](raw_data) and [partially processed results](imageJ_results) are provided here to test the included analysis scripts. For the complete dataset used in the manuscript including ROIs manually defined at cell boundaries, please download the `.zip` file from [DOI](). These files can then be processed using the workflow below from `pixel_calculator.py` onwards. Outputs will be deposited in the 

## Workflow

Initial preprocessing of the raw `LIFF` files was completed in ImageJ, using the Bioformats importer (available via the drag-and-drop interface). Stacked `TIFF` files were then exported for further processing as described below.

1. Eif4a3 

| Script                        | Language/Interpreter | Description                                                                   |
|-------------------------------|----------------------|-------------------------------------------------------------------------------|
| pixel_collector.py            | FiJi                 | Define ROIs, collect per-pixel information                                    |
| pixel_calculator.py           | Python               | Filter pixels to define nuclei, cell ROIs and calculate initial summary stats |
| plot_EIF4A3_histogram.py      | Python               | Generate kernel density estimate for nuclear and cytoplasmic pixels           |
| pixel_intensity_proportion.py | Python               | Calculate proportion of pixels in nucleus above threshold                     |
| plot_eIF4A3_single_cells.py   | Python               | Plot individual cells with thresholded pixels overlayed                       |


2. Fus


| Script                        | Language/Interpreter | Description                                                                     |
|-------------------------------|----------------------|---------------------------------------------------------------------------------|
| pixel_collector.py            | FiJi                 | Define ROIs, collect per-pixel information                                      |
| pixel_calculator.py           | Python               | Filter pixels to define nuclei, cell ROIs and calculate initial summary stats   |
| pixel_ROI_mask.py             | Python               | Generate binary masks using pixel coordinates for cell, nucleus, cytoplasm ROIs |
| spot_counting.cpipe | CellProfiler         | Detect spots, filter using binary masks                                         |
| spot_cleanup.py               | Python               | Filter CellProfiler output, calculate per cell/treatment means                  |
| spot_to_ROI.py                | Python               | Read individual spot ROIs to regenerate binary mask for spots per cell          |
| plot_FUS_single_cells.py      | Python               | Plot individual cells with spots overlayed                                      |


## References

<a name="footnote_1">1.</a> Schindelin, J.; Arganda-Carreras, I. & Frise, E. et al. (2012), "Fiji: an open-source platform for biological-image analysis", Nature methods 9(7): 676-682, PMID 22743772, doi:10.1038/nmeth.2019 (on Google Scholar).

<a name="footnote_2">2.</a> McQuin C, Goodman A, Chernyshev V, Kamentsky L, Cimini BA, Karhohs KW, Doan M, Ding L, Rafelski SM, Thirstrup D, Wiegraebe W, Singh S, Becker T, Caicedo JC, Carpenter AE (2018). CellProfiler 3.0: Next-generation image processing for biology. PLoS Biol. 16(7):e2005970 / doi. PMID: 29969450
