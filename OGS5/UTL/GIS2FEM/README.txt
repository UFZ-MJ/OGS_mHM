

This is a brief instruction of how to convert the triangular wise recharge data 
into the nodal source terms of a three dimensional finite element model for OGS.

Assume the file name of the  the 3D mesh file name is 'model170.msh', 
and the file containing the triangular wise recharge data is called 
'gwR_WesternDS_coordinates1_1977_2010_withoutIDs.txt'.  
 The first line of file 'gwR_WesternDS_coordinates1_1977_2010_withoutIDs.txt' 
must contains three strings for coordinates, then pairs of two strings of times, e.g. 'Okt 77'. 
The second line must be the total number of rows of the following data. 
The data process takes the following procedure.
     

Firstly, the face elements of the top surface have to be subtracted from the three
dimensional mesh and written to a file. This task is performed by using command
 'GIS2FEM' with option 6 as  
   GIS2FEM  6  G:\Code\temp\Agnes\JAMS-Model\model170_redzMatGr\model170
where 'G:\Code\temp\Agnes\JAMS-Model\model170_redzMatGr\' is the path to the directory 
where file 'model170.msh' exists. After executing this command, a file called 
'model170_sfc.msh' is generated.
 
Secondly, a new input file has to be prepared. The file has four lines as:

  Ratio: 1.e-3
  Unit: month
  model170_sfc
  gwR_WesternDS_coordinates1_1977_2010_withoutIDs.txt

The first line has one keyword 'Ratio' and one value. The value can be used to 
scale the recharge data. For instance, the unit of the recharge data is mm/d 
but the unit of m/d is required in the simulation by OGS, then we set 1.e-3 
for 'Ratio'. The second line have a keyword 'Unit' and it is for time unit.
The third line specifies the file name (without extension) of the newly 
generated file that contains the face elements of the top surface. 
While the forth line gives the name of the file that contains the triangular 
wise recharge data. Assuming this new file is saved as 'monthly.txt', the data 
conversion can then be performed as the last step.

Last, with the two files prepared in the previous two steps, the recharge data 
is then converted into the OGS input files by using the command with option '7' as 
   GIS2FEM  7  G:\Code\temp\Agnes\JAMS-Model\model170_redzMatGr\monthly.txt
Program GIS2FEM first reads the top surface elements, which could be triangle or 
quadrilateral, and the  triangular wise recharge data. Then it assign the proper 
recharge value to the top surface elements by checking the coordinates of the centroid 
of each top surface element. For each surface element, if its centroid is in a triangle 
of the grid of  the recharge data, the value of the triangle is assigned to the surface 
element.  Once all top surface elements have been processed, the elements that have been 
assigned with the recharge values are involved in the face integration calculation, 
in which the recharge data is converted into nodal source terms. The calculated results 
are written to different asci files for the different measurement times given in file 
 gwR_WesternDS_coordinates1_1977_2010_withoutIDs.txt.  Besides the files for nodal 
source terms,  one more control file  is generated as well after the command is done. 
The additional file is named automatically by adding string '_monthly.ifl' to the name of 
the recharge data file. In this example, it is
   gwR_WesternDS_coordinates1_1977_2010_withoutIDs.txt_monthly.ifl
The file starts with a keyord 'DIRECT' in the first line, and follows with the names of the generated 
files and their corresponding times, e.g.
DIRECT
14 Dez78.asc
15 Jan79.asc
17 Mrz79.asc
25 Nov79.asc
26 Dez79.asc
27 Jan80.asc
28 Feb80.asc
29 Mrz80.asc
38 Dez80.asc


OGS simulation:
For the simulation, one only needs to fill '$DIS_TYPE' in .st input file with key word 'PRECIPITATION' 
and the name of the control file name that generated  the last step of the above data conversion. 
Here is an example:
 
 #SOURCE_TERM
 $PCS_TYPE
  GROUNDWATER_FLOW
 $PRIMARY_VARIABLE
  HEAD
 $DIS_TYPE
   PRECIPITATION gwR_WesternDS_coordinates1_1977_2010_withoutIDs.txt_monthly.ifl
 #STOP


   
   
 
 
 
 
   

     
      

 