This folder contains the R function corresponding to the NDEXPO application found at www.expostats.ca under the "other tools" tab.

NDEXPO can be references to using the following : 

Lavoué, J (2013) dealing with non-detects in industrial hygiene datasets. Exposure (official magazine of the British occupational hygiene society), December 2013, 13-16.

Note added March 31st 2022 : on the ordering of the original data in the results

This wa prompted by a user who found out that the "original data" section of the results were not the same as the input data. Here is the clarification below.

Very early in the process, the data are ordered in order to try to save the initial order (the procedure plays a lot with reordering). This first ordering is very basic and does not take censoring into account. With increasing sample size, it becomes likely that a tie will happen. The issue occurs when a detected and a non detected values are numerical ties (e.g. 0.15 and <0.15) : NDexpo doesn't know well how to order these.

This does not mean the function is mistaken : in the final data set there will be a 0.15 value that was untouched, and one that was transformed according to ROS. The issue is that the function is not able to reorder them exactly as in the initial sample.

Note added July 2026 : on the values reported as "original data" when few values are detected

When fewer than 3 values are detected (or, for samples of fewer than 5 observations),
ROS is not attempted and each non-detect is replaced by half its reporting limit.

Until July 2026 this substitution was written into the working vector before the results
were assembled, so the substituted values ended up in $x as well as in $xfin. The
"original data" column therefore showed half the reporting limit rather than the limit as
entered - most visibly as a minimum reported at half its true value in datasets with a
high proportion of non-detects.

$x now holds the values as entered and only $xfin carries the substituted ones, which is
what the output description above has always specified. $xfin, $yfin, $pp, $alpha and
$beta are unchanged; anything computed from $x (minimum, maximum, proportion above a
limit) will differ from previous versions, and correctly so.

Further questions : jerome.lavoue@umontreal.ca
