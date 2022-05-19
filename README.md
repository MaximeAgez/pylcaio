# pyLCAIO
An object class to hybridize lifecycle assessment (LCA) and environmentally extended input-output (EEIO) databases.


* Create your own LCA-IO hybrid database (e.g., combining ecoinvent and exiobase data)
* Automates hybridization and correction for double-counting with two available methods (STAM and binary)
* Default parameters only allow the hybridization of ecoinvent 3.5, 3.6, 3.7, 3.7.1 and 3.8 with EXIOBASE3.7+ (v3.7 and higher)
* The resulting hybrid-ecoinvent database can be exported to brightway2 and the GUI activity-browser
* Includes matching of ecoinvent and EXIOBASE environmental flows to Impact World+

Specific additional features, _**only available**_ while hybridization ecoinvent3.5 with exiobase
* Can accept capitals-endogenized version of EXIOBASE
* Includes extrapolated additional environmental extensions for EXIOBASE (from USEEIO)
* Includes _**regionalized**_ characterization matrices for use with Impact World+

This library will be regularly updated to provide support for newer versions of ecoinvent.

# System requirements
Under 12GM of RAM you will most likely run into a MemoryError, making it impossible to generate a database

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
* Agez, M., Majeau-Bettez, G., Margni, M., Strømman, A. H., & Samson, R. (2019). Lifting the veil on the correction of double counting incidents in hybrid Life Cycle Assessment. Journal of Industrial Ecology, 24(3), 517–533. https://doi.org/https://doi.org/10.1111/jiec.12945
* Agez, M., Wood, R., Margni, M., Strømman, A. H., Samson, R., & Majeau-Bettez, G. (2020). Hybridization of complete LCA and MRIO databases for a comprehensive product system coverage. Journal of Industrial Ecology, 24(4), 774–790. https://doi.org/10.1111/jiec.12979
* Agez, M., Muller, E., Patouillard, L., Södersten, C. J. H., Arvesen, A., Margni, M., Samson, R., & Majeau-Bettez, G. (2021). Correcting remaining truncations in hybrid LCA database compilation. Journal of Industrial Ecology. https://doi.org/10.1111/jiec.13132