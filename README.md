
Python scripts and accessory files released to accompany the manuscript entitled “Disruption of ultraconserved elements by copy number variation or chromosome rearrangement is associated with neurodevelopmental phenotypes” by Ruth B. McCole, Wren Saylor, Claire Redin, Chamith Y. Fonseka, Harrison Brand, Jelena Erceg, Michael E. Talkowski, and C.-ting Wu.

Correspondence to the corresponding author:
twu@genetics.med.harvard.edu. 
Department of Genetics,
77 Avenue Louis Pasteur,
NRB 264, Boston,
MA 02115, USA. 
Phone: (617) 432-4431. 
Fax (617) 432-7663.

or the first author:
rmccole@genetics.med.harvard.edu

Distributed under the following license:

Copyright 2017 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

___________________

Accessory text files have suffix .txt or .bed

Python script requirements:
Python version 2.7

Modules required:

argparse

collections

logging

math

matlab.engine

matplotlib

numpy

os

pandas

pybedtools

random

scipy

seaborn

sys
___________________

1. Enrichment analysis (see Methods)

Scripts:
collapsecoordinates.py
coordinateoverlaps.py
enrichment_driver.py
enrichment.py

Accessory files:
hg18_exons_awk1based.txt
hg18_genic_awk1based.bed
hg18_intergenic_awk1based.bed
hg18_introns_awk1based.bed
hg18_nonN_awk1based.txt

Instuctions:

collapsecoordinates.py and coordinateoverlaps.py prepare files containing regions of interest such as CNVs for enrichment analysis. coordinateoverlaps.py requires accessory file hg18_nonN_awk1based.txt. 

Enrichment analysis is performed with enrichment_driver.py. All accessory files are required, in addition to files listing the coordinates of the UCEs of interest.

___________________

2. UCE-to-breakpoint distances analysis (see Methods)

Scripts:
distances1.py
distances2.py
distances3.py
AndersonDarling.py

Instructions:

distances1.py filters a set of breakpoints to retain only one breakpoint from a cluster of a specified proximity on the chromosome. 

distances2.py calculates distances between a set of UCEs, a set of breakpoints, and optionally, between sets of random regions and breakpoints.

distances3.py bins and calculates frequency counts that can be plotted as a histogram of the distances data output from distances2.py.

AndersonDarling.py performs Anderson-Darling test on distances produced by distances2.py.

___________________

3. Partial correlation (see Methods)

Scripts:
partial1.py
partial2.py

Accessory file:
hg18_nonN_0based.txt

Instructions:

partial1.py calculates densities of provided features in preparation for partial correlation analysis. Accessory file hg18_nonN_0based.txt is required.

partial2.py carries out partial correlation analysis.

________________________

General usage:

Each script is equipped with a help message, accessed from the command line by using the argument -h.

Example (command prompt represented by $):

$ python clustermodule.py -h
