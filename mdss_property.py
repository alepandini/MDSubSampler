import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis import distances


class ProteinProperty:
    """
    A class used to calculate a property for the protein based on the frame selection
    from the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    vector:
        Simple vector used to calculate statistics
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = None

    def __init__(self, protein_data, vector, frame_list, atom_selection="name CA"):
        self.protein_data = protein_data
        self.atom_selection = atom_selection
        self.property_vector = []

    def _add_reference_to_protein_data(self):
        """
        Method that links the ProteinProperty and ProteinData classes
        """
        self.protein_data._add_property_dummy(self, self.property_name)

    def discretize_vector(self, min_value=None, max_value=None):
        """
        Method that discretises a vector used for Bhatta and KL distances

        Returns
        ----------------------------
        return the discretised vector
        """
        if min_value is None:
            min_value = self.min_value
        if max_value is None:
            max_value = self.max_value
        bin_size = (max_value - min_value) / 100.0
        bin_vector = np.arange(min_value, max_value, bin_size)
        counts, bins = np.histogram(self.property_vector, bins=bin_vector)
        self.property_vector_discretized = dict(
            zip(bins, (counts / len(self.property_vector)))
        )
        return self.property_vector_discretized

    def _property_statistics(self):
        """
        Method that calculates the minimum, maximum and average values of a vector
        """
        self.min_value = np.min(self.property_vector)
        self.max_value = np.max(self.property_vector)
        self.avg_value = np.average(self.property_vector)

    def set_reference_coordinates(self):
        """
        Method that sets us on the first frame and extracts a copy of the coordinates
        of the first frame only for a selection of atoms
        """
        self.protein_data.trajectory_data.trajectory[0]
        self.ref_coordinates = self.protein_data.trajectory_data.select_atoms(
            self.atom_selection
        ).positions.copy()

    def write_property_vector(self, outfilename):
        """
        Method that saves the vector with the calculations of a specific property for
        a protein in a file
        """
        with open(outfilename, "w") as f:
            for i, value in enumerate(self.property_vector):
                f.write("{} {}\n".format(i, value))

    def write_property_discretised_vector(self, outfilename):
        """
        Method that saves the discretised vector with the calculations of a specific
        property for a protein in a file
        """
        discr_vector = self.discretize_vector()
        with open(outfilename, "w") as f:
            for (key, value) in discr_vector.items():
                f.write("{} {}\n".format(key, value))


class RMSDProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the RMSD value for each frame in the
    protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = "RMSD"

    def __init__(self, protein_data, frame_list, atom_selection="name CA"):
        super().__init__(protein_data, frame_list, atom_selection)
        self.set_reference_coordinates()

        for frame in frame_list:
            """
            Go through the trajectory and for each frame I compare with my reference frame
            """
            self.protein_data.trajectory_data.trajectory[frame]
            self.property_vector.append(
                rms.rmsd(
                    self.protein_data.trajectory_data.select_atoms(
                        atom_selection
                    ).positions,
                    self.ref_coordinates,
                )
            )

        self._property_statistics()
        self.discretize_vector()


class DistanceProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the distance value between two
    atoms or two group of atoms for each frame in the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    def __init__(
        self,
        protein_data,
        frame_list,
        atom_selection_1="name CA",
        atom_selection_2="name CA",
    ):

        selection_of_atoms_1 = protein_data.trajectory_data.select_atoms(
            atom_selection_1
        )

        selection_of_atoms_2 = protein_data.trajectory_data.select_atoms(
            atom_selection_2
        )
        dist = distances.distance_array(
            selection_of_atoms_1.positions, selection_of_atoms_2.positions
        )

        print(type(dist))


class RMSFProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the RMSF value for the particles
    in the atomgroup in each frame in the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = "RMSF"

    def __init__(self, protein_data, frame_list, atom_selection="name CA"):
        super().__init__(protein_data, frame_list, atom_selection)
        self.set_reference_coordinates()

        u = self.protein_data.trajectory_data
        protein = u.select_atoms("protein")
        # probably changes the u so maybe I need to make a copy of it
        prealigner = align.AlignTraj(u, u, select=atom_selection, in_memory=True).run()
        reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
        reference = mda.Merge(protein).load_new(
            reference_coordinates[:, None, :], order="afc"
        )
        aligner = align.AlignTraj(
            u, reference, select="protein and name CA", in_memory=True
        ).run()
        ca = protein.select_atoms(atom_selection)
        print(ca)
        rmsfer = rms.RMSF(ca, verbose=True).run()
        self.property_vector.append(rmsfer.results.rmsf)

        self._property_statistics()
        self.discretize_vector()


class RadiusOfGyrationProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the Radius of Gyration value for each frame
    in the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = "rog"

    def __init__(self, protein_data, frame_list, atom_selection="name CA"):

        super().__init__(protein_data, frame_list, atom_selection)

        self.time = []
        for frame in frame_list:
            protein_data.trajectory_data.trajectory[frame]

            self.time.append(protein_data.trajectory_data.trajectory.time)
            self.property_vector.append(
                protein_data.trajectory_data.select_atoms(
                    atom_selection
                ).radius_of_gyration()
            )

        self._property_statistics()
        self.discretize_vector()