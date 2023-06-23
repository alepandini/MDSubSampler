"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Namir Oues
    
    This file is part of the MDSubSampler software 
    (https://github.com/alepandini/MDSubSampler).
    Copyright (c) 2023 Namir Oues and Alessandro Pandini.

    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
from mdss.property import ProteinProperty
import MDAnalysis.analysis.pca as pca
import numpy as np
from mdss.log_setup import log


class TrjPCAProj(ProteinProperty):
    """
    Subclass of ProteinProperty class representing PCA property. Calculates pca_space
    vector while setting number of columns as number of PCs and number of rows as number
    of frames in the trajectory.

    Attributes
    ----------
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    atom_selection : str, optional
        Atom selection for property calculation. Default selection is CA atoms.
    """

    display_name = "TrjPCAProj"

    def __init__(self, protein_data, atom_selection="name CA"):
        """
        Initialize the TrjPCAProj object.

        Parameters
        ----------
        protein_data : ProteinData
            An instance of the ProteinData class representing the protein data.
        atom_selection : str, optional
            Atom selection for property calculation. Default selection is CA atoms.
        """
        super().__init__(protein_data, atom_selection)
        self.pca_model = None
        self.n_pcs = None
        self.ed_n_pcs = None
        self.property_matrix = None

    def _run_pca(self):
        """
        Perform Principal Component Analysis for a given protein trajectory and selection of atoms.
        Retrieve principal components.
        """
        self.pca_model = pca.PCA(self.protein_data.trajectory_data, self.atom_selection)
        self.pca_model.run()

    def select_pc(self, pc_index=1):
        """
        Select principal components.

        Parameters
        ----------
        pc_index : int, optional
            Principal component index. Default option is index 1.

        Notes
        -----
        If property values are not available, a warning message is printed.
        If the specified PC index is larger than the number of PCs, a warning message is printed.
        The selected principal component is assigned to the `property_vector` attribute.
        The PC index is stored in the `pc_index` attribute.
        """
        if not self.n_pcs:
            print("Warning: values not available. Please calculate property values.")
            log.warning(
                "{:12s} Values not available. Please calculate property values".format(
                    "STEPS"
                )
            )
        else:
            if pc_index > self.property_matrix.shape[1]:
                print("Warning: PC index larger than number of PCs.")
                log.warning(
                    "{:12s} PC index larger than number of PCs.".format("STEPS")
                )
            else:
                self.property_vector = self.property_matrix[:, (pc_index - 1)]
                self.pc_index = pc_index

    def write_property_vector(self, outfilepath):
        if outfilepath is None:
            return

        with open(outfilepath, "w") as f:
            for i, values in zip(self.protein_data.frame_indices, self.property_matrix):
                values_display = " ".join(str(v) for v in values)
                f.write("{} {}\n".format(i, values_display))

    def calculate_property(self, var_threshold=0.8, pc_filter=True):
        """
        Perform Principal Component Analysis with all trajectory frames.

        Parameters
        ----------
        var_threshold : float, optional
            Threshold to compare cumulated variance to pick the right number of principal components.
            Default is 0.8.
        pc_filter : bool, optional
            Flag to indicate whether to use the filtered number of principal components.
            If True, `self.ed_n_pcs` will be used; otherwise, `self.n_pcs` will be used.
        """

        self._run_pca()
        self.n_pcs = self.pca_model.n_components
        self.ed_n_pcs = np.where(
            self.pca_model.results.cumulated_variance > var_threshold
        )[0][0]
        atomgroup = self.protein_data.trajectory_data.select_atoms(self.atom_selection)
        if pc_filter:
            selected_pcs = self.ed_n_pcs
        else:
            selected_pcs = self.n_pcs
        self.property_matrix = self.pca_model.transform(
            atomgroup, n_components=selected_pcs
        )
        self.select_pc()
        return self.property_matrix
