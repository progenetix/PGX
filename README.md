## PGX - Genome visualisation from the Progenetix project

This repository contains code and resources for visualizing different types of genome data, focussed on copy number variants (CNV). The software is mostly used for the [progenetix.org](http://progenetix.org) visualizations; however, users can just download the package & use with their own data (though documentatuion is ... sparse).

Current implementations (in `bin`):

### Genome arrays

* example implementations
  - [arrayplotter.pl](bin/arrayplotter.pl)
  
### Copy number histograms & (clustered) samples

* example implementations
  - [CNVsegfilePlotter.pl](bin/CNVsegfilePlotter.pl)
  - for progenetix and TCGA style segmentation & probe files

#### Example
  
* download `PGX` to whichever location (`___my_path___`)
* run as in example below, replacing the segment and label files with your own

```
cd ___my_path___/PGX/bin/
perl CNVsegfilePlotter.pl -outdir ./out -f ./data/testfile_segments.tab  -sf ./data/testfile_labels.tab -min_group_no 1
```

### Sample profiles from Progenetix-style database


#### Example

* installation as above
* for more about Progenetix see
  - <http://docs.progenetix.org>
  - <http://github.com/progenetix/schemas>
  - <http://github.com/progenetix/bycon>

```
cd ___my_path___/PGX/bin/
perl CNVdbPlotter.pl -query icdom-81 -randno 500
```
