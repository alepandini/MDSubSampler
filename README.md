# MDSubSampler
Molecular Dynamics Sub Sampler

MDSubSampler is a Python library and toolkit for a posteriori subsampling of multiple trajectory data for further analysis. This toolkit implements uniform, random, stratified sampling, bootstrapping and targeted sampling to preserve the original distribution of relevant geometrical or energetic properties.

Input:  Molecular dynamics trajectory 
        Reference structure [optional] 
        Range of subsample sizes (or percentages) 
        Dissimilarity threshold [optional] 
       

Output: .dat file(s) with calculated property for all sample sizes input
        .dat file with calculated property for full trajectory (user input)
        .xtc file(s) with sample trajectory for all sample sizes
        .npy file(s) with sample trajectory for all sample sizes 
        .json file report with important statistics from the analysis
        .txt log file with essential analysis steps
        
