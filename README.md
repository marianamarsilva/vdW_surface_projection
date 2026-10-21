# **Molecular van der Waals Footprint Calculator**



The script calculates the projected van der Waals (vdW) footprint of a molecule using a Boolean union approach. 

The molecular structure is projected onto a user defined plane. The individual atomic vdW spheres are represented as circles using Bondi radii. Overlapping regions are removed using a Boolean union operation.



## Requirements



R version 4.0 and above

Package sf and ggplot2 (in the R console run install.packages("sf"), install.packages("ggplot2") and then library("sf") and library ("ggplot2"))





## Input 



A coordinate file in .txt format (this could be for instance obtained from the .xyz exported from Avogrado).



Example:

"19

&#x09;Energy:      27.1381945

N         -0.67102        0.10721       -0.28641

N          1.44670       -0.20605        0.04064

C         -3.61150       -2.77591       -0.69820

H          0.59738        0.89733       -1.52670

..."



## Procedure



Select the atoms defining the projection plane (e.g., check the number attribute to each atom in the Avogrado projection)

In the R console:

&#x09;Define the working directory (e.g. setwd())					

&#x09;Define the plane atoms (plane\_atoms) 

&#x09;Define the specimen name (specimen\_name)

E.g. 



> plane\_atoms <- c(1,5,7,3)

> specimen\_name <- c("1-ethyl-3-methylcyclohexane")



Run the scrip and selecting the correct .txt file 



## Output



The script generates the following information:

&#x09;- Number of atoms and element counts

&#x09;- Normal vector of the plane

&#x09;- Projected xx and yy radius per atom (vdW radii are assigned according to Bondi\[1])

&#x09;- Boolean union area (in Angstrom squared)

&#x09;- Removed overlapping area

&#x09;- Removed percentage

&#x09;- Converted footprint in square meter

&#x09;- A .png file with the vdW projection of the molecule, with atoms coloured by element, and the area (specimen\_name\_vdW\_projection.png)



## Assumptions and limitations



The molecular geometries must be optimised.

For calculations, Bondi vdW raddi are used \[1], some elements may not be covered.

Area footprints are based on 2D projections.

Plane selection influences the calculated footprint.







## References

\[1] Bondi, A. van der Waals Volumes and Radii. Journal of Physical Chemistry 1964, 68(3), 441-451. DOI: 10.1021/j100785a001

&#x09;- 

