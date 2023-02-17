# MDSubSampler
Molecular Dynamics Sub Sampler

MDSubSampler is a Python library and toolkit for a posteriori subsampling of multiple trajectory data for further analysis. This toolkit implements uniform, random, stratified sampling, bootstrapping and targeted sampling to preserve the original distribution of relevant geometrical or energetic properties.

Input:  

Molecular dynamics trajectory 
<br />Reference structure [optional] 
<br />Range of subsample sizes (or percentages) 
<br />Dissimilarity threshold [optional] 
       

Output: 

.dat file(s) with calculated property for all sample sizes input
<br />.dat file with calculated property for full trajectory (user input)
<br />.xtc file(s) with sample trajectory for all sample sizes
<br />.npy file(s) with sample trajectory for all sample sizes 
<br />.json file report with important statistics from the analysis
<br />.txt log file with essential analysis steps
        
