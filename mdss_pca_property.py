from mdss_property import ProteinProperty
import MDAnalysis.analysis.pca as pca


class PCA(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to perform PCA on the protein
    trajectory and calculates a pca_space vector with number of columns equal
    the number of PCs and the number of rows equal the number of frames in the
    trajectory.

    Returns
    ----------
    pca_vector: A vector that will be used to sample the trajectory of the protein
    """

    display_name = "PCA"

    def calculate_property(self):

        prot_data = self.protein_data.trajectory_data
        protein_selection_pca = pca.PCA(prot_data, self.atom_selection)
        protein_selection_pca.run()
        n_pcs = np.where(protein_selection_pca.results.cumulated_variance > 0.95)[0][0]
        atomgroup = self.protein_data.trajectory_data.select_atoms(self.atom_selection)
        pca_vector = protein_selection_pca.transform(atomgroup, n_components=n_pcs)
        self.property_vector = pca_vector
        return pca_vector

