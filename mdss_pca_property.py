from mdss_property import ProteinProperty
import MDAnalysis.analysis.pca as pca
import numpy as np


class TrjPCA(ProteinProperty):

    display_name = "TrjPCA"

    def __init__(self, protein_data, atom_selection='name CA'):
        super().__init__(protein_data, atom_selection)
        self.pca_model = None
        self.n_pcs = None
        self.ed_n_pcs = None
        self.property_matrix = None

    def _run_pca(self):
        self.pca_model = pca.PCA(self.protein_data.trajectory_data, self.atom_selection)
        self.pca_model.run()

    def select_pc(self, pc_index=1):
        if not self.n_pcs:
            print("Warning: values not available. Please calculate property values.")
        else:
            if pc_index > self.property_matrix.shape[1]:
                print("Warning: PC index larger than number of PCs.")
            else:
                self.property_vector = self.property_matrix[:,(pc_index-1)]
                self.pc_index = pc_index

    def calculate_property(self, var_threshold=0.8, pc_filter=False):
        self._run_pca()
        self.n_pcs = self.pca_model.n_components
        self.ed_n_pcs = np.where(self.pca_model.results.cumulated_variance > var_threshold)[0][0]
        if pc_filter:
            selected_pcs = self.ed_n_pcs
        else:
            selected_pcs = self.n_pcs
        self.property_matrix = self.pca_model.results.p_components[:,:selected_pcs]
        self.select_pc()


class TrjPCAProj(TrjPCA):
    """
    A Subclass of ProteinProperty class used to perform PCA on the protein
    trajectory and calculates a pca_space vector with number of columns equal
    the number of PCs and the number of rows equal the number of frames in the
    trajectory.

    Returns
    ----------
    pca_vector: A vector that will be used to sample the trajectory of the protein
    """

    display_name = "TrjPCAProj"

    def calculate_property(self, var_threshold=0.8, pc_filter=False):
        self._run_pca()
        self.n_pcs = self.pca_model.n_components
        self.ed_n_pcs = np.where(self.pca_model.results.cumulated_variance > var_threshold)[0][0]
        atomgroup = self.protein_data.trajectory_data.select_atoms(self.atom_selection)
        if pc_filter:
            selected_pcs = self.ed_n_pcs
        else:
            selected_pcs = self.n_pcs
        self.property_matrix = self.pca_model.transform(atomgroup, n_components=selected_pcs)
        self.select_pc()

