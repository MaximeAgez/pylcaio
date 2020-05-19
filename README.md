# pyLCAIO
An object class that can structure, manipulate, and facilitate the hybridization of lifecycle assessment (LCA) and environmentally  extended input-output (EEIO) matrices.

* Only works with Python
* Read in, combine, organize, manipulate and concatenate LCA foreground and background matrices
* Combine LCA system with EEIO matrices
* Automate hybridization and correction for double-counting with two available methods (STAM and binary)
* Can currently hybridize ecoinvent3.5 with exiobase3
* Includes extrapolated additional environmental extensions for EXIOBASE
* Includes matching of ecoinvent and EXIOBASE to Impact World+
* Can accept capitals-endogenized version of EXIOBASE
* Includes regionalized characterization matrices for use with Impact World+

The software is still under development but has already 2 operational versions.

# System requirements
Memory: 8GB RAM

# Dependencies
* Python 3
* Pandas
* Numpy
* Scipy
* pymrio
* ecospold2matrix
* pickle

# Related publications
* Majeau-Bettez, G., Agez, M., Wood, R., Södersten, C., Margni, M., Strømman, A. H., & Samson, R. (2017). Streamlined Hybridization software: merging Ecoinvent and Exiobase. In Biennial Conference of the International Society for Industrial Ecology.
* Agez, M., Majeau-Bettez, G., Margni, M., Strømman, A. H., & Samson, R. (2019). Lifting the veil on the correction of double counting incidents in hybrid Life Cycle Assessment. Journal of Industrial Ecology, 1–17. https://doi.org/https://doi.org/10.1111/jiec.12945
* Agez, M., Wood, R., Margni, M., Strømman, A. H., Samson, R., & Majeau-Bettez, G. (2019). Hybridization of complete LCA and MRIO databases for a comprehensive product system coverage. Journal of Industrial Ecology, 1–17. https://doi.org/10.1111/jiec.12979




