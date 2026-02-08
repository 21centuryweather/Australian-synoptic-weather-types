# Australian Synoptic Weather Types

This repository houses the Australian Synoptic Weather Types (SWTs), an index of the daily circulation across the Australian continent. The code herein can be used to recreate these types from ERA5 as well as assign other gridded 850hPa wind data to a synoptic weather type.  

30 synoptic weather types are defined which are categorised into 8 groups. 5 of the 8 SWTs identify the position of large anticyclones within the domain reflecting the different phasing of extratropical wave packets over the Australian region. We refer to these SWTs as the western high (WH), central high (CH), eastern high (EH), Tasman high (TH) and flanking high (FH) regimes. 2 of the 8 groups of SWTs are comprised of closed cyclonic circulations at the surface. The most frequent of the two closed cyclonic circulation SWTs is the cut-off low (COL) regime, representing a deep cut-off with a surface expression. The second is characterized by a strong trough over Western Australia, which extends from the tropics over the western part of the continent, promoting a quasi-closed circulation south of the southern coastline, called west coast trough (WCT) SWTs. The final set of SWTs is the active Australian monsoon (AM) SWT, which represents 60\% of the total summer climatology. It is with tropical moist margin lies over the northern part of the Australian continent and well-defined monsoon trough over northern Australia with westerly flow over or north of the maritime continent.

![Synoptics_SWT_ERA5-2](https://github.com/user-attachments/assets/d436d2b1-0046-4b26-9ae2-88c9f713d7e7)

More infromation on the development and interpretation of the SWTs can be found here: 
Barnes, Michael A., Franciscus Liqui Lung, Christian Jakob, Andrew Gunn, and Michael J. Reeder. 2025. “Australian Synoptic Weather Types.” Journal of Geophysical Research: Atmospheres 130 (24): e2025JD043873. https://doi.org/10.1029/2025JD043873.

# Using the SWTs

SWT climatology: An updated climatology of the SWTs can be found in SWT_climatology in either csv or netcdf format.
SWT reassignment: Any gridded 850hPa wind field can be assigned to a SWT. 

For more information contact:
Michael Barnes (michael.barnes@monash.edu) or Christian Jakob (christian.jakob@monash.edu)
