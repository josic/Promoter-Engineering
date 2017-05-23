# Promoter-Engineering
R Code for the data analysis in "Tuning the dynamic range of bacterial promoters regulated by ligand-inducible transcription factors" by Chen, et al

The two excel files contain the data that is imported by the two R files. 

The inducible promoter data is oragnized according to two main "Conditions": Ara or Las.
The third and fourth column in the spreasdsheet determine -35 and -10 condition (in that order), 
and the fith column gives the measured promoter strength.


There are four conditions for the hybrid promoter data, as explained in more detail in the 
manuscript and supplementary material. The four conditions are ++, +−, −+, and −− with
the first sign denoting inducer absent/present, and the second sign 
denoting repressed state (ie repressed = -, unrepressed = +).  In the spreasheet

1 -  is the -/- state in the notes, repressed and no activator
2 -  is the +/- state in the notes, repressed with activator
3 -  is the -/+ state in the notes, uninduced
4 -  is the +/+ state in the notes, activator, but no repression

The R files contain the code that was used to fit the model to the data by minimizing
the objective function specified in the Supplentary material. 
