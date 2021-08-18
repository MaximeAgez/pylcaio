#### This folder contains all information and support module needed to compile a newer version of hybrid-ecoinvent.

Follow the notebook "How to update hybrid-ecoinvent" it will guide you through the needed steps.

###### More information on the modifications needed when compiling a newer version:

- need an eco-x-_exio-y- folder in the src/Data/ folder for the newer version. Note that the code reads the version of 
ecoinvent declared by the variable entered when creating the object of the class DatabaseLoader so make sure that numbers
of version and format match.

- within the previously mentioned folder the same txt and xlsx files are needed as for previous versions. Some of these 
files might not to be completed, while others always have to be modified:
    - **_countries.txt_** regroups the countries of exiobase. If exiobase did not add new countries, no modification needed
    - **_countries_per_region.txt_** regroups the regions used by ecoinvent (e.g., RER) and translates them in terms of 
    countries covered by exiobase (with exiobase countries), e.g., RNA is translated to CA and US.
    - **_Filter.xlsx_** this file dictates which processes of ecoinvent will be hybridized or not. As unique identifiers
    of ecoinvent processes change with each new version of ecoinvent, this file **must** be updated with every new version.
    See the appropriate notebook to see how to update this file.
    - **_geography_replacements_** contains the countries (and not regions since they are in countries_per_region.txt)
    that are defined by ecoinvent, but that are not explicitly covered by exiobase, e.g., there are processes for Thailand
    in ecoinvent, but Thailand is not explicitly in exiobase as it is regrouped with other asian countries in the exiobase
    region WA. This file must thus be updated if new ecoinvent versions include new processes with new countries not 
    explictly covered by exiobase. See the appropriate notebook to see how to update this file.
    - **_Product_concordances.xlsx_** matches ecoinvent processes to exiobase sectors to allow the hybridization. The 
    productIds and activityNameIds of ecoinvent processes, used to match processes to sectors, do not change with newer
    versions of ecoinvent. Nevertheless, the metadat associated to these identifiers can change. For example, from 
    ecoinvent3.4 to ecoinvent3.5, "green bell pepper" became "bell pepper", but the identifier for this process 
    remain unchanged. While the metadata itself is not used for the matching, it is recommended to update the metadata 
    using the appropriate notebook. 
    In addition, with newer version of ecoinvent come new processes and very likely new products which need to be matched
    to exiobase sectors. This matching cannot realistically be done by a machine and will thus have to be done manually.
    Support functions exist though to identify which are the new products covered by ecoinvent. In summary, this file 
    **must** be updated with each newer version of ecoinvent. See the appropriate notebook to see how to update this file.
    - **_STAM_categories_** regroups exiobase sectors per bigger categories and is required for correction of double
    counting. Unless exiobase introduces new sectors or modifies the names of some sectors this file does not need to 
    be updated.
    - **_STAM_function_categories_** identifies particular exiobase sectors and is required for correction of double
    counting. Unless exiobase introduces new sectors or modifies the names of some sectors this file does not need to 
    be updated.

- Need to update the inflation value of euros if the year selected with exiobase is not currently supported. Currently
supported years: 1995 to 2021