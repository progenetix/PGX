## PGX - Genome visualisation from the Progenetix project

This repository contains code and resources for visualizing different types of genome data, focussed on copy number variants (CNV). The software is mostly used for the [progenetix.org](http://progenetix.org) visualizations; however, users can just download the package & use with their own data (though documentatuion is ... sparse).

Current implementations (in `bin`):

#### Genome arrays

* example implementations
  - [arrayplotter.pl](bin/arrayplotter.pl)
  
#### Copy number histograms & (clustered) samples

* example implementations
  - [segfilePlotter.pl](bin/segfilePlotter.pl)
  - for progenetix and TCGA style segmentation & probe files
  

### Example use
  
* download `PGX` to whichever location (`___my_path___`)
* run as in example below, replacing the segment and label files with your own

```
cd ___my_path___/PGX/bin/
perl segfilePlotter.pl -outdir ./out -f ./data/testfile_segments.tab  -sf ./data/testfile_labels.tab -min_group_no 1
```
