# Methodology of pyLCAIO

This document describes in details the methodology that is followed by pyLCAIO to hybridize ecoinvent with exiobase.

The methodology is divided into multiple sections:
1. [General methodology](#general-methodology)
2. [Which processes to hybridize](#which-processes-to-hybridize)
3. [Special case for electricity and heat](#special-case-for-electricity-and-heat)
4. [Get price data](#get-price-data)
5. [Mapping our way through](#mapping-our-way-through)
6. [Uncorrected upstream cutoff](#uncorrected-upstream-cutoff)
7. [Determining covered inputs](#determining-covered-inputs)
8. [Correcting for double counting](#correcting-for-double-counting)
9. [Writing into brightway2](#writing-into-brightway2)

Electricity vs heat, how to know which is internal or to-be-hybridized

### General methodology
For each process to hybridize in ecoinvent, we will identify its corresponding exiobase sector, then add **_all_** inputs
from the IO sector to the LCA description, using the price provided of the process provided by ecoinvent to convert from
€ to physical units. Doing so necessarily double counts inputs (i.e., inputs that are covered by both the LCA and the IO
part). To correct for this, we determine which inputs are covered by the LCA and then remove these inputs from the IO 
description, to only add from IO what is missing from LCA.

### Which processes to hybridize
Not all processes within the ecoinvent database must be hybridized. Doing so would trigger another type of double counting
where the added IO inputs (a.k.a upstream cutoffs) are added multiple times for the same transformation/production step.<br>
First, market processes must not be hybridized. These processes do not represent a transformation/fabrication of a product,
but simply represent a market distribution of the technologies to produce a product. Furthermore, their price is identical
to the price of the production processes. Hybridizing those would thus add the upstream cutoffs twice.


Second, empty or aggregated (S) processes must not be hybridized neither. Indeed, if the process is empty in ecoinvent, 
it reflects the cutoff aspect of the ecoinvent database where secondary products come free of environmental burdens
(think about thr recycled content cut-off process). Hybridizing such a process would simply add environmental burdens
which is inconsistent with how ecoinvent operates. For aggregated processes, we only have access to life cycle 
elementary flows which thus makes it impossible to know which inputs are already covered by LCA and which ones should
be missing, according to IO. It is thus impossible to hybridize them.


Lastly, processes which we refer to as "internal processes" or "activities" in ecoinvent, should also not be hybridized.
These processes are "fake" processes that ecoinvent created to dismember production processes into many small steps 
which in reality are all performed within a single process in a single factory. They also often do not represent the
real output of a sector. Think for example of the process "drying of feed grain". This is a process that describes the
drying step of feed grains. In reality though, this drying is executed by the same company that produces the feed
grain. In the corresponding IO sector, "Cereal grains nec", yes the drying step must be accounted for within the sector
but it is clearly not the main "output" of the sector, meaning that the average "use" of the sector will not be driven by
the drying step but most likely by the production of the grain itself. Hybridizing the drying of the grain with the 
"Grain" sector overall would not make sense. Moreover, more often than not, there are no specific prices for these 
"internal processes". Their hybridization is thus not performed. Their identification was performed mostly manually. In 
the end, which processes are "to hybridize" or "not to hybridize" can be seen in the filters.xlsx files, in the 
Data/ecoinvent/ei3.x/mappings/ folders.

### Special case for electricity and heat
For electricity and heat, not all processes are to be hybridized neither. Only processes which actually produced
energy that is then distributed to users through a grid must be hybridized. Energy processes that are re-used internally
within the factory must not be hybridized. But how do we determine which energy processes are being sold to other
consumers and which ones are only consumed internally? To do so we rely on the capacity of the generating apparatus.
If the latter exceeds 10MV, it is considered to be for the purpose of actually selling energy. If it is lower, then it
is most likely a small operation that only focuses on re-using energy within a facility. There are three exceptions:
- if the description f the ecoinvent process specifically states that the energy produced is sent to the grid 
(for electricity) or a district heat network (for heat), then this process is hybridized
- heat and power from wood is hybridized as heat nd electricity are the only two outputs of the process, so the operation
is running solely to produce (and thus sell) at the very least electricity. We assume that heat is also sold instead of
being used internally to produce electricity.
- photovoltaic installation of 570kWp are hybridized. That is because their use within ecoinvent is to form a large
PV farm which goal is obviously to sell electricity.

### Get price data
To connect information of exiobase (in €) to the processes of ecoinvent (in various physical units), we require prices
for each of the processes that were identified for hybridization. We simply took the prices provided by the ecoinvent 
database itself for each of their processes. These prices can only be found directly in the .ecospold files themselves. 
They are not visible on any LCA software or directly visible on the ecoinvent website. These prices are old. They are
all expressed in euros of the year 2005. We thus converted those to the reference year of the exiobase version used,
through a simple inflation rate.

These prices are also very coarse. For instance, all electricity production processes have the same price of 0.0977€/kWh,
disregarding technologies of production and countries, which makes little sense.

There is probably better price data out there. For instance, price data could be extracted from the UN COMTRADE or BACI 
databases, even though it would only provide a price paid at the border, which might be different from the price of 
commodities sold domestically for example. One issue to handle though (amongst many), is that there is not a 1-for-1 
concordance with commodities of the trade database (using the HS6 classification) and ecoinvent processes. So not all
processes could have prices extracted this way.

### Mapping our way through
A mapping between ecoinvent processes and exiobase sectors is required. This mapping was produced manually and is 
available in the filters.xlsx files within the Data/ecoinvent/ei3.x/mappings/ folders.

### Uncorrected upstream cutoff
Once we have collected all necessary information,w e can finally start the automated hybridization. We simply loop 
through the processes needing hybridization and for each, we identify its corresponding sector and geography, through 
the mapping files, which allows us to identify the specific column of exiobase to copy paste. In the case of processes
from ecoinvent from broad regions (e.g., GLO, RER, RoW, IAI Area, Asia, without China and GCC, etc.) we determined which
countries do these regions correspond to, and recreated an aggregate of the corresponding countries of exiobase, using 
the distribution of these regions. So if we have the production of 1-butanol in RER. We look at all European countries
of exiobase and more specifically to their production of "Chemicals nec". Then, if the distribution of their output
is as such Germany 35%, Belgium 25%, France 20%, Italy 10%, Spain 10%, then we will use these output shares to recreate 
a "Chemicals nec" sector for the RER region.

Note that this is also done for the dynamic region of RoW (Rest of the World) which covered countries changes depending 
on what specific country the ecoinvent database covers for a particular commodity. So for 1-butanol, ecoinvent
could cover one specific country A and RoW would then be all countries in the world except A. But for wheat, ecoinvent
could cover countries A and B and thus RoW would be all countries in the world except A and B.

So, by simply either selecting the corresponding sector in exiobase, or the aggregate of sectors (for broad regions) and
multiplying these inputs of these sectors by the price of the process, we obtain the uncorrected upstream cutoffs, which
are all added to the ecoinvent process.

### Determining covered inputs
After adding all the uncorrected upstream cutoffs, there are necessarily double counted inputs. To determine which inputs
are double counted or not, we first need to determine which inputs are covered by LCA, in terms of IO sectors. This step
is trickier than we can think. The ecoinvent structure is in our way. Let me explain.

Take the example of an agricultural process in ecoinvent, such as the "wheat grain production". What inputs do you 
directly have in this process? Ok we have some pesticides and fertilizers, some irrigation, some chemicals. No trace of
a tractor though. Why is that? That's because the tractor is actually not a direct input of this process. It's coming 
along with the process "tillage, ploughing" or through "sowing". If we were to simply look at the directly covered 
inputs, the code would thus not see any tractor, while there are tractors used in the "Wheat" sector of exiobase, and it
would thus add tractor inputs to this process. This would count the tractor twice, as the tractor would be a direct input,
and would also be added through "tillage" or "sowing". So how to remediate to this kind of issue, where characteristic
inputs are not directly available within the process, but are essentially "passed" through other technosphere inputs?

To do that, we can use a mathematical matrix trick that will automatically add these characteristic inputs directly to 
the process itself. So the trick is that these processes that "pass" characteristic inputs are the "internal processes"
that we mentioned in the previous section. In essence, we thus need to add all inputs of the "internal processes" to 
their final product, that is, the next process identified as needing hybridization. This also allows adding inputs of 
"market for" processes to their final products, i.e., transportation to each product.

SO what is the mathemetical matrix trick to do this? We do 2 copies of the technosphere matrix of ecoinvent. In the
1st copy, we will force all processes to-hybridize to zero, call this matrix A_not_hyb. A_not_hyb thus only contains
the inputs of processes that are ot to hybridize, such as markets or internal processes. In the 2nd matrix, we do the
opposite, thus creating A_hyb, which only contains the inputs of processes-to-hybridize. Then, we invert A_not_hyb,
thus connecting life cycle stages of non-hybridized inputs together. Multiplying (dot product) A_hyb with the inverse 
of A_not_hyb then add all inputs of non-hybridized inputs to the different hybridized inputs using them. That way, the 
tractor input from "tillage" is now directly contained in the description of "wheat grain production" inside A_hyb.

Finally, we can then simply look at the inputs covered by LCA for each process within A_hyb, and this tells us which
inputs not to add from IO, because they are already covered by LCA.

### Correcting for double counting
After determining the covered inputs, we need to actually correct for the double counting. In its simplest form, this
corresponds to only deleting upstream cutoffs corresponding to whatever is already covered by LCA. This is what we 
typically refer to as the "binary" method and this is one of the two arguments possible for the .correct_double_counting()
function.

The second option is the STAM (Similar Technological Attribute Method). This second approach has the same starting point 
as the binary, i.e., simply deleting all IO inputs already covered by LCA, but it also includes heuristics to refine
the correction even more. One of the issue with the simple binary approach is that sometimes missing inputs in LCA 
descriptions are not missing data but rather a specificity of a technology. Take the electricity car operation. No diesel
or gasoline inputs in there. But that's not because it was forgotten. Simply because it's an electric car. The binary
approach would see such a car and add inputs of fossil fuels to its description because the average "Motor vehicles" in
exiobase has fossil fuel inputs. the STAM is further detailed in this article https://doi.org/https://doi.org/10.1111/jiec.12945.

### Writing into brightway2
Once all the data is available in the brightway2 dictionaries format, it's time to import all the relevant databases
into the SQL database itself. While it could be a simple bw2.Database().write() this would import the full version of
exiobase and of the hybridized version of ecoinvent. Here, pyLCAIO chooses to impose a cutoff to limit the amount of 
data written in the database. This kind of goes against the principle of hybrid LCA, where we want to avoid such cutoffs,
but we do it as a practical choice. The full version of exiobase and hybrid ecoinvent would includesomething like 8000
to 10000 inputs per process of ecoinvent/exiobase. You can imagine that calculations are longer as a result, especially
if you are performing a "FT contributions" on the activity-browser. As a result, we implemented by default two cutoffs,
which can be changed by the user.

First, we have a "culling_threshold" argument for the .import_exiobase_in_bw2() function. By default, it was set to
1e-7€. Meaning that any input of exiobase smaller than 1e-7€ is not imported in the brightway2 database. Of course, 
1e-7€ is absolutely nothing, but the aggregation of all these small values actually still impact the result by about 
0.5% on the climate change indicator when we compared results. If users are not comfortable with this variation, they
can lower this culling_threshold, or even set it to zero to not have any cutoff at all.

In a second time, we have a "cutoff" argument for the .import_hybridized_database() function. This cutoff works a little
bit differently than the culling_threshold of exiobase, simply because we have to work in relative terms here. We cannot 
cull anything lower than 1e-7€ in a hybridized process of ecoinvent if the price of the process overall is 1e-6€. 
Instead, the cutoff is defined as a % of the price. By default, it is set to 0.00001, meaning that upstream cutoffs that
are smaller than 0.0001% of the price of the process are not imported in the database. Obviously, care was taken to
not remove LCA inputs that would not cross that threshold, which could maybe happen with extrmely small inputs such as
infrastructure (4e-16 unit factory ot something).
