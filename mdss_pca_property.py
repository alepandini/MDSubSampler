from mdss_property import ProteinProperty
import MDAnalysis.analysis.pca as pca
import numpy as np


class TrjPCA(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to perform PCA on the protein
    trajectory and calculates a pca_space vector with number of columns equal
    the number of PCs and the number of rows equal the number of frames in the
    trajectory.

    Returns
    ----------
    pca_vector: A vector that will be used to sample the trajectory of the protein
    """

    display_name = "TrjPCA"

    def _run_pca(self):
        self.pca_model = pca.PCA(self.protein_data.trajectory_data, self.atom_selection)
        self.pca_model.run()

    def calculate_property(self, var_threshold = 0.8, pc_filter = False):
        self._run_pca()
        self.n_pcs = self.pca_model.n_components
        self.ed_npcs = np.where(self.pca_model.results.cumulated_variance > var_threshold)[0][0]
        atomgroup = self.protein_data.trajectory_data.select_atoms(self.atom_selection)
        if pc_filter:
            selected_pcs = self.ed_npcs
        else:
            selected_pcs = self.n_pcs
        self.property_vector = self.pca_model.transform(atomgroup, n_components=selected_pcs)

