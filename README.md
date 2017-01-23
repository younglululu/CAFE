# **CAFE**
**aCcelerated Alignment-FrEe sequence analysis**

===================

Thank you for downloading this tool for sequence distance/dissimialrity measures using state-of-art Alignment-Free methods. This software provides the well-optimized programs to compute overall 27 measures including (1) Conventional measures based on k-mer counts, (2) Newly developed measures based on background adjusted k-mer counts, and (3) Measures based on presence/absence of k-mers. The detailed definitions can be found in the paper or in the later section. 

CAFE works with sequence data, both long genomic sequences and shotgun sequence reads from NGS technologies, and subsequently generates pairwise dissimilarities among the sequences as output. CAFE provides four types of visualized downstream analysis, including heatmap, two dimensional projection using principal coordinate analysis (PCoA), network display and sequence clustering into a dendrogram by using the neighbour-joining algorithm. All the analysis can be performed by simply clicking through well-designed graphical user interface (GUI) or invoking a stand-alone command line executable program on three common operating systems (Linux, Mac, and Windows)


----------

Compatibility
=============

The stand-alone command line executable program is written in c++, which has been fully tested under gcc version 4.8.1. It works under Windows, Linux or Mac environment. Precompiled executables are provided for these platforms, and users have the option to compile the source code for their specific platform if desired (see the Installation section below). The graphical user interface (GUI) is written in Python Tkinter library, the built-in GUI library for all standard Python distributions.

----------

Dependencies
============

Fundamental dependencies
------------------------

> Python v2.7.*
> gcc

These items are prerequisites for the installation of CAFE as described below. The former is essential for the graphical user interface and the later is essential for the stand-alone command line executable program.

Python package
------------------------

> numpy
> scipy 
> biopython
> scikit-learn
> matplotlib
> networkx

These are the python packages that need to be installed in order to run the graphical user interface. If the user is only running the stand-alone command line executable program, please just ignore this step.

Additional software (optional)
------------------------

The software relies on JELLYFISH (http://www.cs.cmu.edu/~ckingsf/software/jellyfish/) to count kmer efficiently. With built-in kmer counting program, CAFE can still proceed to work **WITHOUT** the installation of JELLYFISH in a slower manner.

----------

Installation
============

Installation of the stand-alone executable program
------------------------
####<i class="icon-cog"> Compile program directly using g++

```sh
$ make
```

Installation of python dependencies using Anaconda 
------------------------

Anaconda is a tool to isolate your python installation, which allows you to have multiple parallel installations using different versions of different packages, and gives you a very convenient and fast way to install the most common scientific python packages. 

To install Anaconda on your computer, perform the following steps:

 1. Download Anaconda for your platform from [here](https://www.continuum.io/downloads)
 2. After the file is completely downloaded, install Anaconda:
	- Windows users can double-click on the installer and follow the on-screen instruction
	- Mac users can double-click the .pkg file and follow the instructions displayed on screen
	- Linux users can run the following command:
```sh
$ bash <downloaded_file>
```
More detailed installation instructions can be found [here](https://docs.continuum.io/anaconda/install.html)

After installing Anaconda, including the corresponding dependent python package varies on different systems, and described in this README is only how to proceed with a linux (ubuntu) distribution. 

Create a new environment:
```sh
$ conda create -n cafe_env python=2.7.6
```

After choosing to proceed, run the suggested command:
```sh
$ source activate cafe_env 
```

Then install the concoct dependencies into this environment:
```sh
$ conda install numpy scipy biopython scikit-learn matplotlib networkx
```

Usage
=====

Usage of Graphical User Interface
------------------------

The main window has the layout shown in the following figure, containing five parts in terms of functionality:

 1. The red area contains four button involving add/remove sequence data files. The sequence data files can be either long genomic sequences or shotgun sequence reads from NGS technologies, with the file extension '.fasta', '.fa' or '.fna'. The first button indicates adding one single sequence data file. The second button indicates adding all sequence data files in the specified directory. The third button indicates removing the selected sequence data files shown in the green panel below. And the fourth button indicates removing all sequence data files shown in the green panel below.
 2. The yellow area contains the configuration might involved in various distance measures. For example, the selection of 27 distance measures, the selection of kmer length, the selection of possible markov order describing the sequence generating model, the threshold to cutoff the kmer occurrence, and whether to consider the reverse complementary of each kmer, which is a common choice in the scenario of dealing with shotgun sequence reads from NGS technologies
 3. The green area contains  all added sequence data files as mentioned.
 4. The blue area keeps tracking of the running information when calculating the distance measures. 
 5. The purple area are the key to visualize the pairwise distance result. Specifically,  CAFE provides four types of visualized downstream analysis, including heatmap, two dimensional projection using principal coordinate analysis (PCoA), network display and sequence clustering into a dendrogram by using the neighbour-joining algorithm. Each analysis is shown in the respective tabbed window.

<p align="center">
  <img src="https://raw.githubusercontent.com/younglululu/CAFE/master/example/visualization/gui_snapshot0.jpg"/>
</p>


Here we go through a test example by step-by-step guidance.  You can find a folder named "example" in the package, which contains 30 virus genomes. There are three subfolders, "data", "hash" and "visualization". The folder "data" contains the corresponding genome sequence files. The folder "hash" contains the corresponding binary kmer count files. The folder "visualization" contains the visualization-related files including web pages and associated javascript/css files.

We first add the virus genome sequence files by specifying the directory.

<p align="center">
  <img src="https://raw.githubusercontent.com/younglululu/CAFE/master/example/visualization/gui_snapshot1.jpg"/>
</p>

We nest specify the directory of pre-computed binary kmer count files to save the time. The setting panel contains the configuration of Jellyfish executable file path, the saved hash directory and the quantile of edges to display in the network analysis.

<p align="center">
  <img src="https://raw.githubusercontent.com/younglululu/CAFE/master/example/visualization/gui_snapshot2.jpg"/>
</p>

After necessary information has been set, here we use D2star distance measures, we simply click the 'Run' button. From the following four snapshots, we can see the running information has been displayed, and the pairwise distance measures have been further processed into multiple way. Here is the dendrogram by using the neighbour-joining algorithm.

<p align="center">
  <img src="https://raw.githubusercontent.com/younglululu/CAFE/master/example/visualization/gui_snapshot3.jpg"/>
</p>

Here is the two dimensional projection using principal coordinate analysis (PCoA).

<p align="center">
  <img src="https://raw.githubusercontent.com/younglululu/CAFE/master/example/visualization/gui_snapshot4.jpg"/>
</p>

Here is the heatmap, very straightforward.

<p align="center">
  <img src="https://raw.githubusercontent.com/younglululu/CAFE/master/example/visualization/gui_snapshot5.jpg"/>
</p>

Here is the network analysis with respect to the 10% quantile of the edges with smallest distance as weight.

<p align="center">
  <img src="https://raw.githubusercontent.com/younglululu/CAFE/master/example/visualization/gui_snapshot6.jpg"/>
</p>



Usage of  Stand-alone Executable Program
------------------------

> **Command:  ** ./cafe [options]* -D  < dist > -I < fa_files > -K  < intK >

> - Main arguments:
	- **-D** < dist >  Comma-separated list of distance measurements,  **E.g.** -D D2star,Ma,CVtree. The options include: 
		 - Conventional measures based on kmer counts :			 
				 1. Ch: Chebyshev distance
				 2. Canberra: Canberra distance
				 3. Chisq: Chi-Square distance
				 4. Cosine: Cosine distance
				 5. Co-phylo: Co-phylog distance
				 6. D2: D2 distance
				 7. Eu: Euclidean distance
				 8. FFP: Feature frequency profiles (FFP)
				 9. JS: Jensen-Shannon divergence
				 10. Ma: Manhattan distance
				 11. Pearson: Pearson distance
		 - Newly developed measures based on background adjusted kmer counts: 
				 1. CVtree: CVtree distance
				 2. D2shepp: D2shepp distance
				 3. D2star: D2star distance
		 - Measures based on presence/absence of kmers:
				 1. Anderberg: Anderberg distance
				 2. Antidice: anti-Dice distance
				 3. Dice: Dice distance
				 4. Gower: Gower distance
				 5. Hamman: Hamman distance
				 6. Hamming: Hamming distance
				 7. Jaccard: Jaccard distance
				 8. Kulczynski: Kulczynski distance
				 9. Matching: Matching distance
				 10. Ochiai: Ochiai distance
				 11. Phi: Pearson Phi distance
				 12. Russel: Russel-Rao distance
				 13. Sneath: Sneath-Sokal distance
				 14. Tanimoto: Rogers-Tanimoto distance
				 15. Yule: Yule distance
	- **-I** < fa_files > Comma-separated list of sequence fasta files, e.g. -I speciesA.fa,speciesB.fa,speciesC.fa. Pairwise similarity is calculated based upon the sequences specified with this option.
	- **-K** < intK > Kmer Length.

> - Options:
	- **-J** < jfexe_path > Use jellyfish to accelerate kmer counting. <jfexe_path> denotes the file path of jellyfish executable file, e.g. jellyfish-2.2.4/bin/./jellyfish
	- **-L** < lower > Only consider k-mer with occurrence >= <lower>. The default value is 1.
	- **-M** < order > Markov Order involved in D2star, D2shepp and JS. There are two possible options. The first option is one single value indicating that all the sequences use the same order. The second option is comma-separated list of orders. Notice that the length of the list should match the number of fasta files. The order value could be non-negative integer but less than Kmer length or \"-1\" with the special intention to automatically infer the suitable order (not suitable for JS). The default Markov Order is 0 (i.i.d. model).
	- **-R** Consider Reverse Complement in kmer counting.
	- **-S** < dir > Save/Load calculated k-mer count binary files to the folder < dir >. Each input fasta file corresponds to particular model.
	- **-O** < path > Output results to file at < path >.
	- **-T** < type > The output type as the input to downstream analysis, including: plain, [phylip](http://evolution.genetics.washington.edu/phylip.html) (as hierarchical clustering), [cytoscape](www.cytoscape.org/) (as network analysis) and mds (Multidimensional Scaling as 2D plotting). E.g. -T mds. The default type is plain.

> - Examples:
	- ./cafe -M 0 -O output_path -S model_dir -T plain -I speciesA.fa,speciesB.fa -J /panfs/cmb-panasas2/ylu465/jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star,Ma,Hao
	- ./cafe -M 0 -S model_dir -I speciesA.fa,speciesB.fa -J /panfs/cmb-panasas2/ylu465/jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star,Ma,Hao
	- ./cafe -M 0 -L 2 -I speciesA.fa,speciesB.fa -J /panfs/cmb-panasas2/ylu465/jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star,Ma,Hao -R


A test example
==============

You can find a folder named "example" in the package, which contains 30 virus genomes. There are three subfolders, "data", "hash" and "visualization". The folder "data" contains 30 genome fasta files. The folder "hash" contains the corresponding binary kmer count files. The folder "visualization" contains the visualization-related files including web pages and associated javascript/css files.



----------


Contacts and bug reports
========================

Please send bug reports, comments, or questions to 

Yang Lu: [ylu465@usc.edu](mailto:ylu465@usc.edu)

Prof. Fengzhu Sun: [fsun@usc.edu](mailto:fsun@usc.edu)


----------

Copyright and License Information
=================================


Copyright (C) 2016 University of Southern California, Yang Lu

Authors: Yang Lu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.

Last update: 29-Dec-2016
