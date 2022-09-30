# def run_subsampler(p_data, property_class, sampler_class):
#     """
#     Method that uses the user input for property calculation, sampling method
#     and size of sample and returns a subsample trajectory along with a log
#     file with diagnostics for the particular property."""
#     log.logging.info("The MDSubsampler is running:")
#     property.calculate_property()
#     property_sample.calculate_property()
#     print(f"Calculating {property_class.display_name}")
#     print(f"Applying {property_class.display_name}")
#     property.calculate_property()


# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "11%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='BootstrappingSampler' --n-iterations=50 --size=11000 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "0.1%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1991 --size=100 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "0.5%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1992 --size=500 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "1%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1993 --size=1000 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "5%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1994 --size=5000 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "10%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1995 --size=10000 --dissimilarity='BhattaCoefficient'

# Results

# 0.05: 0.4847413787731183
# 0.1:0.1892646854226224
# 0.5:0.03408998524296803
# 1: 0.008814553862456592
# 5:0.003142312466437116
# 10:0.001085055552200252