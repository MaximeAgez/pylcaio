# pyLCAIO
An object class to hybridize lifecycle assessment (LCA) and environmentally extended input-output (EEIO) databases.


* Create your own LCA-IO hybrid database (e.g. combining ecoinvent and exiobase data)
* Automates hybridization and correction for double-counting with two available methods (STAM and binary)
* Default parameters only allow the hybridization of ecoinvent3.5 with EXIOBASE
* Can accept capitals-endogenized version of EXIOBASE
* Includes extrapolated additional environmental extensions for EXIOBASE (from USEEIO)
* Includes matching of ecoinvent and EXIOBASE environmental flows to Impact World+
* Includes regionalized characterization matrices for use with Impact World+
* Can be exported to brightway2

If you are just interested in the default hybrid database (if you do not want to or cannot run the code) you can find it here: https://zenodo.org/record/3890379

This software is still under development.

# System requirements
Under 8GM of RAM you will most likely run into a MemorryError, making it impossible to generate a database

# Dependencies
* Python 3
* Pandas
* Numpy
* Scipy
* pymrio
* ecospold2matrix
* pickle
* brightway2
* bw2agg

# Related publications
* Majeau-Bettez, G., Agez, M., Wood, R., Södersten, C., Margni, M., Strømman, A. H., & Samson, R. (2017). Streamlined Hybridization software: merging Ecoinvent and Exiobase. In Biennial Conference of the International Society for Industrial Ecology.
* Agez, M., Majeau-Bettez, G., Margni, M., Strømman, A. H., & Samson, R. (2019). Lifting the veil on the correction of double counting incidents in hybrid Life Cycle Assessment. Journal of Industrial Ecology, 1–17. https://doi.org/https://doi.org/10.1111/jiec.12945
* Agez, M., Wood, R., Margni, M., Strømman, A. H., Samson, R., & Majeau-Bettez, G. (2019). Hybridization of complete LCA and MRIO databases for a comprehensive product system coverage. Journal of Industrial Ecology, 1–17. https://doi.org/10.1111/jiec.12979




