# BPPS-tools
---

### Downstream analysis of mcBPPS/omcBPPS sequence comparison results.

___pttrn2pdb.py:___ maps residue patterns identified from mcBPPS/omcBPPS onto PDB structures and outputs a PyMol script (.pml).

___annt_parti.py:___ analyzes sequence partitions output from mcBPPS/omcBPPS and annotates each partition in an hpt file with: 1) the # of sequence hits to families in a MAPGAPS profile and 2) the lowest common taxonomic group represented by sequences in the partition.

___mcbpps2annt_tsv.py:___ outputs residue pattern data from a .lpr file as a .tsv file 

Requirements
---

### pttrn2pdb.py
* __Biopython__ 
* __mapgaps__ (http://www.igs.umaryland.edu/labs/neuwald/software/mapgaps/)
* __MAPGAPS sequence profile (Default: ePKf)__ (http://www.igs.umaryland.edu/labs/neuwald/software/mapgaps/)

### annt_parti.py
* __Biopython__ 
* __ete3__ 
* __tweakcma__ (http://www.igs.umaryland.edu/labs/neuwald/software/mapgaps/)
* __mapgaps__ (http://www.igs.umaryland.edu/labs/neuwald/software/mapgaps/)
* __MAPGAPS sequence profile (Default: ePKf)__ (http://www.igs.umaryland.edu/labs/neuwald/software/mapgaps/)

These tools were designed for downstream analysis of mcBPPS and omcBPPS progrmas, which were developed by Dr. Andrew F. Neuwald. Relevant software can be found here:

MAPGAPS, new version of BPPS: http://www.igs.umaryland.edu/labs/neuwald/software/.

mcBPPS: http://www.chain.umaryland.edu/mcbpps/index.html.

omcBPPS: http://www.chain.umaryland.edu/omcbpps/.



