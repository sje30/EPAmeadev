# EPAmeadev

This is the repository to accompany our 2016 paper freely available at
http://dx.doi.org/10.1177/1087057116640520


    Characterization of Early Cortical Neural Network Development in
    Multiwell Microelectrode Array Plates

    Ellese Cotterill, Diana Hall, Kathleen Wallace,
    William R Mundy, Stephen J Eglen, Timothy J Shafer
	
	Journal of Biomolecular Screeing (2016)

You are free to use any of the data or resources in this repository.
We do request however that if you use this material, you cite the
above paper in any work that you publish.


	
## Data files

The data files for this project are stored in
[allH5Files](allH5Files/).  They are stored in the HDF5 format as
outlined in our
[2014 paper](http://dx.doi.org/10.1186/2047-217X-3-3).  Our
[analysis package, sjemea](http://github.com/sje30/sjemea) can read in these files as follows:

```{r}
library(sjemea)
s = h5.read.spikes("allH5Files/CO_20140212_MW1007-51_DIV12_00_001.h5")
fourplot(s)
```
These are the first multi-well datasets that we have published, so
some of the data in the HDF5 fields has not previously been
documented.  The new fields include:

1. `/dose` dosage of any drug (not used here) to each well.
2. `/treatment` any manipulation performed to a well (all "control" in
   this study).

## TODO

1. Upload rest of analysis code
2. Push to zenodo


