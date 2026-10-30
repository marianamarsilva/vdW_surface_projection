# **Molecular van der Waals Footprint Calculator**



This R script calculates the projected van der Waals (vdW) footprint of a molecule using a Boolean union approach. 

The molecular structure is projected onto a user defined plane. 
Individual atoms are represented as circles with a radii corresponding to their vdW radii (Bondi radii\[1]).
Overlapping regions are removed using a Boolean union operation to calculate the corrected milecular footprint area.



##### Requirements



R version 4.0 and above

Package sf and ggplot2.





##### Input 



The script requires a coordinate file in .txt format.
Such files can be generated from an optimised molecular structure, for example by exporting an .xyz file from Avogrado.



Example:

"19

&#x09;Energy:      27.1381945

N         -0.67102        0.10721       -0.28641

N          1.44670       -0.20605        0.04064

C         -3.61150       -2.77591       -0.69820

H          0.59738        0.89733       -1.52670

..."



##### Procedure



Identify the atoms that define the projection plane.

In the R console:


&#x09;Define the plane atoms (plane\_atoms), where numbers corresopnf to the atom indices in the coordinate file.

&#x09;Define the specimen name (specimen\_name)

Important: At least 3 atoms are required to define the reference plane.

E.g. 


> plane\_atoms <- c(plane_atoms <- c(1,2,10,11,13))

specimen_name <- c("1-octyl-3-methylimidazolium")



Run the scrip and select the desired coordinate .txt file when prompted.

The script will automatically:

&#x09;- Import the molecular coordinates;
&#x09;- Fit a least-squares plane through the selected atoms;
&#x09;- Project all atoms onto the plane;
&#x09;- Assign vdW radii;
&#x09;- Calculate overlaped Boolean union area;
&#x09;- Generate a graphical representation of the projection footprint.


##### Output


The script generates the following information

&#x09;- Number of atoms and elemental composition

&#x09;- Plane normal vector

&#x09;- Projected atomic coordinates

&#x09;- Assigned vdW radii (Bondi radii\[1])

&#x09;- Boolean union area (in Angstrom squared)

&#x09;- Removed overlapping area

&#x09;- Removed percentage

&#x09;- A .png file with the vdW projection of the molecule, with atoms coloured by element, and the area (specimen\_name\_vdW\_projection.png)



##### 

##### Assumptions and limitations



The molecular geometries must be optimised.

For calculations, Bondi vdW raddi are used \[1], some elements may not be covered.

Area footprints are based on 2D projections.

Plane selection influences the calculated footprint.




##### 

##### References

\[1] Bondi, A. van der Waals Volumes and Radii. Journal of Physical Chemistry 1964, 68(3), 441-451. DOI: 10.1021/j100785a001

&#x09;- 

