# class AsterError in libaster


class AsterError(Exception):
    pass

    # Method resolution order:
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object

    # Data descriptors defined here:


# class ConvergenceError in libaster


class ConvergenceError(AsterError):
    pass

    # Method resolution order:
    #     ConvergenceError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object


# class IntegrationError in libaster


class IntegrationError(AsterError):
    pass

    # Method resolution order:
    #     IntegrationError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object


# class SolverError in libaster


class SolverError(AsterError):
    pass

    # Method resolution order:
    #     SolverError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object


# class ContactError in libaster


class ContactError(AsterError):
    pass

    # Method resolution order:
    #     ContactError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object


# class TimeLimitError in libaster


class TimeLimitError(AsterError):
    pass

    # Method resolution order:
    #     TimeLimitError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object


# built-in function raiseAsterError in libaster


def raiseAsterError(idmess="VIDE_1", valk=[], vali=[], valr=[]):
    pass


# class PythonBool in libaster


class PythonBool:
    """Enumeration that represents an extended boolean."""

    # Method resolution order:
    #     PythonBool
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    FALSE = 0

    NONE = -1

    TRUE = 1


# class DataStructure in libaster


class DataStructure:
    pass

    # Method resolution order:
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def addDependency(self, ds):
        """Add a dependency to a *DataStructure*.

        Arguments:
            ds (*DataStructure*): Parent *DataStructure* to depend on.
        """

    def build(self):
        """Update the *DataStructure* attributes from the *Jeveux* objects.
        *Only use internally after calling fortran subroutines*.

        Returns:
            bool: *True* if all went ok, *False* otherwise.
        """

    def debugPrint(self, unit=6, synchro=True):
        """Print the raw content of a *DataStructure* on the selected file.

        Args:
            unit (int): File number (default: 6, means stdout).
            synchro (bool): To synchronize prints between processors (default: True).
        """

    def getDependencies(self):
        """Return the explicit dependencies.

        Returns:
            list[*DataStructure*]: List of parents (dependencies) *DataStructure*.
        """

    def getName(self):
        """Return the internal (*Jeveux*) name of the *DataStructure*.

        Returns:
            str: Internal/*Jeveux* name.
        """

    def getTitle(self):
        """Return the tile of the *DataStructure* .

        Returns:
            str: Title of the *DataStructure*.
        """

    def getType(self):
        """Return the name of the *DataStructure* type.

        Returns:
            str: Name of the *DataStructure* type.
        """

    def id(self):
        """Return the identity of the object.

        Returns:
            int: Identifier (address as int).
        """

    def removeDependency(self, ds):
        """Remove a dependency to a *DataStructure*.

        Arguments:
            ds (*DataStructure*): Parent *DataStructure* to be removed from
                dependencies.
        """

    def resetDependencies(self):
        """Clear the list of explicit dependencies."""

    def setTitle(self, title):
        """Set the tile of the *DataStructure* .

        Arguments:
            title [str]: Title of the *DataStructure*.
        """

    # ----------------------------------------------------------------------
    # Data descriptors defined here:

    @property
    def userName(self):
        """str: Name of the user variable that holds this object."""


# class DSWithCppPickling in libaster


class DSWithCppPickling(DataStructure):
    pass

    # Method resolution order:
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""


# built-in function debugJeveuxContent in libaster


def debugJeveuxContent(arg0):
    pass


# built-in function debugJeveuxExists in libaster


def debugJeveuxExists(arg0):
    pass


# built-in function use_count in libaster


def use_count(*args, **kwargs):
    """Overloaded function.

    1. use_count(arg0: Mesh) -> int

    2. use_count(arg0: Model) -> int

    3. use_count(arg0: DOFNumbering) -> int

    4. use_count(arg0: ElementaryMatrix<double, (PhysicalQuantityEnum)4>) -> int

    5. use_count(arg0: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>) -> int

    6. use_count(arg0: ElementaryMatrix<double, (PhysicalQuantityEnum)6>) -> int

    7. use_count(arg0: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>) -> int

    8. use_count(arg0: AssemblyMatrix<double, (PhysicalQuantityEnum)4>) -> int

    9. use_count(arg0: AssemblyMatrix<std::complex<double>, (PhysicalQuantityEnum)4>) -> int

    10. use_count(arg0: AssemblyMatrix<double, (PhysicalQuantityEnum)6>) -> int

    11. use_count(arg0: AssemblyMatrix<std::complex<double>, (PhysicalQuantityEnum)6>) -> int

    12. use_count(arg0: AssemblyMatrix<double, (PhysicalQuantityEnum)5>) -> int

    13. use_count(arg0: AssemblyMatrix<std::complex<double>, (PhysicalQuantityEnum)5>) -> int
    """


# class PhysicalQuantityManager in libaster


class PhysicalQuantityManager:
    pass

    # Method resolution order:
    #     PhysicalQuantityManager
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    # ----------------------------------------------------------------------
    # Static methods defined here:


# class Node in libaster


class Node:
    pass

    # Method resolution order:
    #     Node
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def getId(self):
        """Return the Id of the node.

        Returns:
            int: local id of the node.
        """

    def getValues(self):
        """Return coordinates as (x,y,z.)

        Returns:
            list[float]: (x,y,z).
        """

    def x(self):
        """Return coordinate x.

        Returns:
            float: x.
        """

    def y(self):
        """Return coordinate y.

        Returns:
            float: y.
        """

    def z(self):
        """Return coordinate z.

        Returns:
            float: z.
        """


# class EntityType in libaster


class EntityType:
    """Enumeration for entity type."""

    # Method resolution order:
    #     EntityType
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    AllMeshEntitiesType = 2

    CellType = 3

    GroupOfCellsType = 1

    GroupOfNodesType = 0

    NoType = 5

    NodeType = 4


# class MeshEntity in libaster


class MeshEntity:
    pass

    # Method resolution order:
    #     MeshEntity
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def getNames(self):
        pass

    def getType(self):
        pass


# class AllMeshEntities in libaster


class AllMeshEntities(MeshEntity):
    pass

    # Method resolution order:
    #     AllMeshEntities
    #     MeshEntity
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""


# class BaseMesh in libaster


class BaseMesh(DataStructure):
    pass

    # Method resolution order:
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def build(self):
        """Build list of Tables based on the mesh

        Returns:
            bool: true if building is ok
        """

    def check(self, tolerance):
        """Check some properties of the mesh.

        Arguments:
            tolerance (float): Tolerance used to detect flat cells.
        """

    def getCellType(self, index):
        """Return the type of the given cell

        Arguments:
            index (int) : index of the cell (0-based)

        Returns:
            int : the cell type
        """

    def getCellTypeName(self, index):
        """Return the type name of the given cell

        Arguments:
            index (int) : index of the cell (0-based)

        Returns:
            str : name of the cell type (stripped)
        """

    def getConnectivity(self):
        """Return the connectivity of the mesh as Python lists.

        Returns:
            list[list[int]]: List of, for each cell, a list of the nodes indexes.
        """

    def getCoordinates(self):
        """Return the coordinates of the mesh.

        Returns:
            MeshCoordinatesField: Field of the coordinates.
        """

    def getDimension(self):
        """Return the dimension of the mesh.

        Returns:
            int: 2 or 3
        """

    def getLocalToGlobalCellIds(self):
        """Returns local to global IDs mapping for cells

        Returns:
            list[int]: local to global IDs mapping.
        """

    def getLocalToGlobalNodeIds(self):
        """Returns local to global node Ids mapping

        Returns:
            list[int]: local to global IDs mapping.
        """

    def getMedCellsTypes(self):
        """Return the Med type of each cell.

        Returns:
            list[int]: List of Med types.
        """

    def getMedConnectivity(self):
        """Return the connectivity of the mesh as Python lists following the Med IDs.

        Returns:
            list[list[int]]: List of, for each cell, a list of the nodes indexes.
        """

    def getMinMaxEdgeSizes(self, arg0):
        """Get minimum and maximum length of edges in group of cells

        Returns:
            tuple(real): values of min and max edges
        """

    def getNumberOfCells(self):
        """Return the number of cells of the mesh.

        Returns:
            int: Number of cells.
        """

    def getNumberOfNodes(self):
        """Return the number of nodes of the mesh.

        Returns:
            int: Number of nodes.
        """

    def getOriginalToRestrictedCellsIds(self):
        """If the mesh is created as restriction of an initial mesh,
        It returns a dict between the cell id of the initial mesh and the current cell id.

        Returns:
            dict[int]: a dict between the cell id of the initial mesh and the current cell id.
        """

    def getOriginalToRestrictedNodesIds(self):
        """If the mesh is created as a restriction of an initial mesh,
        It returns a dict betweenn the node id of the initial mesh and the current node id.

        Returns:
            dict[int]: a dict betweenn the node id of the initial mesh and the current node id.
        """

    def getRestrictedToOriginalCellsIds(self):
        """If the mesh is created as restriction of an initial mesh,
        It returns for each cells, the cell id of the initial mesh.

        Returns:
            list[int]: for each cells, the cell id of the initial mesh.
        """

    def getRestrictedToOriginalNodesIds(self):
        """If the mesh is created as a restriction of an initial mesh,
        It returns for each nodes, the node id of the initial mesh.

        Returns:
            list[int]: for each nodes, the node id of the initial mesh.
        """

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """

    def hasCellsOfType(self, type):
        """Return True if mesh contains at least one cell of given type

        Arguments:
            type (str) : cell type

        Returns:
            bool : *True* if mesh contains at least one cell of given type, else *False*
        """

    def isConnection(self):
        """Function to know if a mesh is a ConnectionMesh"""

    def isIncomplete(self):
        """Tell if the mesh is complete on parallel instances.

        Returns:
            bool: *False* for a centralized or parallel mesh, *True* for an incomplete mesh.
        """

    def isParallel(self):
        """Tell if the mesh is distributed on parallel instances.

        Returns:
            bool: *False* for a centralized mesh, *True* for a parallel mesh.
        """

    def printMedFile(self, fileName, local=True, version=[0, 0, 0], mode="a"):
        """Print the mesh in the MED format

        Arguments:
            filename (Path|str): Name of the file
            local (bool=True) : print local values only (relevant for a ParallelMesh only)
            version (list): list of size 3 ([major, minor, release])
            mode (str): 'a' or 'w' ('w' to overwrite file)

        Returns:
            Bool: True if of
        """

    def show(self, verbosity=1):
        """Show mesh informations.

        Arguments:
            verbosity (int): Verbosity level (default: 1)
        """

    def updateInternalState(self):
        """Update the internal state of the datastructure.

        Returns:
            bool: *True* in case of success, *False* otherwise.
        """


# class Mesh in libaster


class Mesh(BaseMesh):
    pass

    # Method resolution order:
    #     Mesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Mesh) -> None

        2. __init__(self: libaster.Mesh, arg0: str) -> None
        """

    def addCellLabels(self, cell_labels):
        """Add cell labels.

        Arguments:
            cell_labels (list) : Cell labels.
        """

    def addNodeLabels(self, node_labels):
        """Add node labels.

        Arguments:
            node_labels (list) : Node labels.
        """

    def convertToBiQuadratic(self, info=1):
        """Convert the mesh to a bi-quadratic one.
        For cells that have no bi-quadratic version, the quadratic version is used.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            Mesh: the bi-quadratic mesh.
        """

    def convertToCubic(self, info=1):
        """Convert the mesh to a cubic one.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            Mesh: the cubic mesh.
        """

    def convertToLinear(self, info=1):
        """Convert the mesh to a linear one.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            Mesh: the linearized mesh.
        """

    def convertToQuadratic(self, info=1):
        """Convert the mesh to a quadratic one.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            Mesh: the quadratic mesh.
        """

    def fix(
        self,
        remove_orphan=True,
        positive_measure=True,
        outward_normal=True,
        double_nodes=True,
        double_cells=True,
        tole=1e-07,
        info=1,
    ):
        """Fix potential problems.

        Arguments:
            remove_orphan (bool) : remove orphelan nodes.
            positive_measure (bool) : reorder nodes to have a positive measure of cells.
            outward_normal (bool) : reorder nodes to have an outward normal for boundary faces.
            double_nodes (bool) : merge double nodes with almost same coordinates.
            double_cells (bool) : merge double cells with same nodes.
            tole (float) : tolerance for double nodes
            info (int) : verbosity mode (0 or 1 or 2).

        Returns:
            Mesh: fixed mesh
        """

    def getCells(self, *args, **kwargs):
        """Overloaded function.

        1. getCells(self: libaster.Mesh, group_name: str) -> list[int]


        Return the list of the indexes of the cells that belong to a group of cells.

        Arguments:
            group_name (str): Name of the local group.

        Returns:
            list[int]: Indexes of the cells of the local group.


        2. getCells(self: libaster.Mesh, groups_name: list[str] = []) -> list[int]


        Return the list of the indexes of the cells that belong to the groups of cells.

        Arguments:
            groups_name (str): Name of the local groups.

        Returns:
            list[int]: Indexes of the cells of the local groups.
        """

    def getGroupsOfCells(self, local=False):
        """Return the list of the existing groups of cells.

        Returns:
            list[str]: List of groups names (stripped).
        """

    def getGroupsOfNodes(self, local=False):
        """Return the list of the existing groups of nodes.

        Arguments:
            local (bool): not used (for compatibilty with ParallelMesh)

        Returns:
            list[str]: List of groups names (stripped).
        """

    def getInnerNodes(self):
        """Return the list of the indexes of the nodes in the mesh

        Returns:
            list[int]: Indexes of the nodes.
        """

    def getOctreeMesh(self, nb_max_pt=1, nb_max_level=20):
        """Get the octree mesh.

        Arguments:
            nb_max_pt (int) : maximum number of points for the last level.
            nb_max_level (int) : maximum number of level.

        Returns:
            Mesh: octree mesh.
        """

    def hasGroupOfCells(self, group_name, local=False):
        """The group exists in the mesh

        Arguments:
            group_name (str): Name of the group.
            local (bool): not used (for compatibilty with ParallelMesh)

        Returns:
            bool: *True* if exists, *False* otherwise.
        """

    def hasGroupOfNodes(self, group_name, local=False):
        """The group exists in the mesh

        Arguments:
            group_name (str): Name of the group.
            local (bool): not used (for compatibilty with ParallelMesh)

        Returns:
            bool: *True* if exists, *False* otherwise.
        """

    def isQuadratic(self, local=False):
        """Tells if the mesh contains quadratic cells.

        Arguments:
            local (bool): not used (for compatibilty with ParallelMesh)

        Returns:
            bool: *True* if the mesh contains quadratic cells, *False* otherwise.
        """

    def printMedFile(self, fileName, local=True, version=[0, 0, 0], mode="a"):
        """Print the mesh in the MED format

        Arguments:
            filename (Path|str): Name of the file
            local (bool=True) : print local values only (relevant for a ParallelMesh only)
            version (list): list of size 3 ([major, minor, release])
            mode (str): 'a' or 'w' ('w' to overwrite file)

        Returns:
            Bool: True if of
        """

    def readAsterFile(self, filename):
        """Read a mesh file from ASTER format.

        Arguments:
            filename (Path|str): Path to the file to be read.

        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """

    def readGibiFile(self, filename):
        """Read a mesh file from GIBI format.

        Arguments:
            filename (Path|str): Path to the file to be read.

        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """

    def readGmshFile(self, filename):
        """Read a mesh file from GMSH format.

        Arguments:
            filename (Path|str): Path to the file to be read.

        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """

    def setGroupOfCells(self, group_name, cell_ids):
        """Set new group of cells in the mesh

        Arguments:
            group_name (str): Name of the new group.
            cell_ids (list[int]) : cell ids which are in the group
        """

    def setGroupOfNodes(self, group_name, node_ids, localNumbering=False):
        """Set new group of nodes in the mesh

        Arguments:
            group_name (str): Name of the new group.
            node_ids (list[int]) : node ids which are in the group
            localNumbering=false (bool): not used (for compatibilty with ParallelMesh)
        """


# built-in function getMedCouplingConversionData in libaster


def getMedCouplingConversionData(mesh):
    """Return three dictionnaries containing data to create an equivalent MedCoupling unstructured mesh.

    MedCoupling needs a mesh splitted by dimension for what concerns cells and groups of cells.
    The group of nodes all belongs to an unique level so there is no need to split them.

     - The first dictionnary (cells) contains for each dimension (the keys) :
       1. The connectivity
       2. The connectivity index
     - The second dictionnary (groups_c) contains for each dimension (the keys) a dictionnary
       which keys are the groups names at the items their cells.
     - The third dictionnary (groups_n) contains for each group of nodes (the keys)
       the nodes composing the group.

    Arguments:
        mesh (BaseMeshPtr): The aster mesh to be processed.

    Returns:
        tuple (cells, groups_c, groups_n) : The data to create the equivalent MedCoupling mesh.
    """


# class DiscreteComputation in libaster


class DiscreteComputation:
    pass

    # Method resolution order:
    #     DiscreteComputation
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, arg0):
        pass

    def __setstate__(self, arg0):
        pass

    def getAcousticDirichletBC(self, time_curr=0.0):
        """Return the imposed acoustic vector used to remove imposed DDL
        *for internal use - prefer to use getDirichletBC*

        Arguments:
              time_curr (float): Current time (default 0.0)

        Returns:
              FieldOnNodesComplex: imposed accoustic vector
        """

    def getAcousticImposedDualBC(self, assembly=True):
        """Return the acoustic imposed nodal BC elementary vector

        Arguments:
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorPressureComplex: imposed dual vector
        """

    def getAcousticNeumannForces(self, assembly=True):
        """Return the elementary acoustic Neumann forces vector

        Arguments:
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorPressureComplex: elementary Neumann forces vector
        """

    def getAcousticVolumetricForces(self, assembly=True):
        """Return the elementary acoustic volumetric forces vector

        Arguments:
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorPressureComplex: elementary volumetric forces vector
        """

    def getCompressibilityMatrix(self, groupOfCells=[]):
        """Return the elementary matrices for compressibility acoustic matrix.
        Option MASS_ACOU.

        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """

    def getContactForces(
        self, geom, displ_prev, displ_step, time_prev, time_step, data, coef_cont, coef_frot
    ):
        """Compute contact and friction forces

        Arguments:
            geom (MeshCoordinatesField): coordinates of mesh used to compute normal
            displ_prev (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            data (FieldOnCellsReal): contact data
            coef_cont (FieldOnNodesReal) : contact coefficient
            coef_frot (FieldOnNodesReal) : friction coefficient

        Returns:
            FieldOnNodesReal: contact and friction forces
        """

    def getContactMatrix(
        self, geom, displ_prev, displ_step, time_prev, time_step, data, coef_cont, coef_frot
    ):
        """Compute contact matrix

        Arguments:
            geom (MeshCoordinatesField): coordinates of mesh used to compute normal
            displ_prev (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            data (FieldOnCellsReal): contact data
            coef_cont (FieldOnNodesReal) : contact coefficient
            coef_frot (FieldOnNodesReal) : friction coefficient

        Returns:
            ElementaryMatrixDisplacementReal: contact and friction elementary matrix
        """

    def getDualElasticStiffnessMatrix(self):
        """Return elementary matrices for dual mechanical BC

        Returns:
            ElementaryMatrix: elementary matrices
        """

    def getDualForces(self, U):
        """Return the imposed displacement assembled vector

        Arguments:
              U (FieldOnNodes): current displacement vector

        Returns:
              FieldOnNodes: dual reaction vector (B^T*lambda)
        """

    def getDualLinearConductivityMatrix(self):
        """Return elementary matrices for dual thermal BC

        Returns:
            ElementaryMatrix: elementary matrices
        """

    def getDualLinearMobilityMatrix(self):
        """Return elementary matrices for dual acoustic BC

        Returns:
            ElementaryMatrix: elementary matrices
        """

    def getDualPrimal(self, U, scaling=1.0):
        """Return the Dirichlet load vector

        Arguments:
              U (FieldOnNodes): current displacement vector

        Returns:
              FieldOnNodes: Dirichlet load vector
        """

    def getElasticStiffnessMatrix(
        self, time_curr=0.0, fourierMode=-1, varc_curr=None, groupOfCells=[], with_dual=True
    ):
        """Return the elementary matrices for elastic Stiffness matrix.
        Option RIGI_MECA.

        Arguments:
              time_curr (float): Current time for external state variable
                evaluation (default: 0.0)
              fourierMode (int): Fourier mode (default: -1)
              varc_curr (FieldOnCellsReal): external state variables at current time
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it empty, the full model is used
              with_dual (bool): compute dual terms or not (default: True)
        Returns:
              ElementaryMatrix: elementary elastic Stiffness matrix
        """

    def getExternalStateVariablesForces(
        self,
        time_curr,
        varc_curr,
        varc_prev=None,
        vari_curr=None,
        stress_prev=None,
        mode=0,
        assembly=True,
        mask=None,
    ):
        """Compute load from external state variables

        Arguments:
              time_curr (float): Current time
              varc_curr (FieldOnCellsReal): external state variables at current time
              varc_prev (FieldOnCells): external state variables at begin of current time
              vari_curr (FieldOnCellsReal): internal state variables at current time
              stress_prev (FieldOnCellsReal): stress at begin of current time
              mode (int): fourier mode
              assembly (bool) : assemble or not
              mask (FieldOnCellsLongPtr): mask to assemble

        Returns:
              FieldOnNodes: load from external state variables
        """

    def getFluidStructureMassMatrix(self, varc_curr=None, groupOfCells=[]):
        """Return the elementary matrices for fluid-structure mass matrix.
        Option MASS_FLUI_STRUC.

        Arguments:
              varc_curr (FieldOnCellsReal): external state variables at current time
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it empty, the full model is used
        Returns:
              ElementaryMatrixReal: elementary fluid-structure mass matrix
        """

    def getFluidStructureStiffnessMatrix(self, fourierMode=-1, varc_curr=None, groupOfCells=[]):
        """Return the elementary matrices for fluid-structure stiffness matrix.
        Option RIGI_FLUI_STRUC.

        Arguments:
              fourierMode (int): Fourier mode (default: -1)
              varc_curr (FieldOnCells): internal state variables at current time
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it empty, the full model is used
        Returns:
              ElementaryMatrixReal: elementary fluid-structure Stiffness matrix
        """

    def getGeometricStiffnessMatrix(
        self, sief_elga, strx_elga=None, displ=None, modeFourier=-1, groupOfCells=[]
    ):
        """Return the elementary matrices for geometric Stiffness matrix.
        Option RIGI_MECA_HYST.

        Arguments:
            sief_elga (FieldOnCellsReal) : stress at Gauss points
            strx_elga (FieldOnCellsReal) : stress at Gauss points for structural element
            displ (FieldOnNodesReal) : displacement field
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixComplex: elementary geometric rigidity matrix
        """

    def getGyroscopicDampingMatrix(self, groupOfCells=[]):
        """Return the elementary matrices for gyroscopic damping matrix.
        Option MECA_GYRO.

        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixReal: elementary gyroscopic damping matrix
        """

    def getGyroscopicStiffnessMatrix(self, groupOfCells=[]):
        """Return the elementary matrices for gyroscopic Stiffness matrix.
        Option RIGI_GYRO.

        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixReal: elementary gyroscopic rigidity matrix
        """

    def getHystereticStiffnessMatrix(self, stiffnessMatrix, varc_curr=None, groupOfCells=[]):
        """Return the elementary matrices for viscoelastic Stiffness matrix.
        Option RIGI_MECA_HYST.

        Arguments:
            stiffnessMatrix : elementary stiffness matrix
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixComplex: elementary viscoelastic rigidity matrix
        """

    def getImpedanceBoundaryMatrix(self, groupOfCells=[], onde_flui=1):
        """Return the elementary matrices for impedance (mechanical) matrix.
        Option IMPE_MECA.

        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
            onde_flui (int): integer to indicate if we have an outgoing or incoming wave
        Returns:
            ElementaryMatrixReal: impedance mechanical matrix
        """

    def getImpedanceMatrix(self, arg0):
        """Return the elementary matrices for impedance (acoustic) damping matrix.
        Option AMOR_ACOU.

        Returns:
            ElementaryMatrixReal: elementary damping matrix
        """

    def getImpedanceWaveMatrix(self, groupOfCells=[]):
        """Return the elementary matrices for impedance (mechanical) matrix
        from an harmonic wave.
        Option ONDE_FLUI.

        Returns:
            ElementaryMatrixReal: impedance wave matrix
        """

    def getIncrementalDirichletBC(self, time_curr, disp):
        """Return the incremental imposed displacement vector used to remove imposed DDL
        for incremental resolution.

        incr_disp = getDirichletBC(time_curr) - disp, with 0.0 for DDL not imposed

        Arguments:
              time_curr (float): Current time
              disp (FieldOnNodes): displacement field at current time

        Returns:
              FieldOnNodes: incremental imposed displacement vector
        """

    def getInternalMechanicalForces(
        self,
        displ_prev,
        displ_step,
        stress,
        internVar,
        internVarIter,
        time_prev,
        time_step,
        varc_prev=None,
        varc_curr=None,
        groupOfCells=[],
    ):
        """Compute internal forces (integration of behaviour)

        Arguments:
            displ_prev (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): field of internal state variables at begin of current time
            internVarIter (FieldOnCells): field of internal state variables at begin of
                                          current newton iteration
            time_prev (float): time at begin of the step
            time_step (float): delta time between begin and end of the step
            varc_prev (FieldOnCells): external state variables at begin of current time
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.

        Returns:
            tuple (tuple): return code error (FieldOnCells),
            error code flag (int),
            internal state variables VARI_ELGA (FieldOnCells),
            Cauchy stress SIEF_ELGA (FieldOnCells),
            field of internal forces (FieldOnNodesReal),
        """

    def getInternalThermalForces(self, temp_prev, temp_step, varc_curr=None, groupOfCells=[]):
        """Compute internal thermal forces (integration of behaviour)
        Option RAPH_THER.

        Arguments:
            temp_prev (FieldOnNodes): thermal field at begin of current time
            temp_step (FieldOnNodes): field of increment of temperature
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used

        Returns:
            tuple (tuple):
            error code flag (int),
            fluxes FLUX_ELGA (FieldOnCellsReal),
            internal forces (FieldOnNodesReal),
        """

    def getLinearCapacityMatrix(self, time_curr, varc_curr=None, groupOfCells=[]):
        """Return the elementary matrices for linear Capacity matrix in thermal computation.
        Option MASS_THER.

        Arguments:
            time_curr (float): current time to evaluate rho_cp
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """

    def getLinearConductivityMatrix(
        self, time_curr, fourierMode=0, varc_curr=None, groupOfCells=[], with_dual=True
    ):
        """Return the elementary matices for linear thermal matrix.
        Option RIGI_THER.

        Arguments:
              time_curr (float): Current time
              fourierMode (int): Fourier mode (default: -1)
              varc_curr (FieldOnCellsReal): external state variables at current time
              groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
              with_dual (bool): compute dual terms or not (default: True)
        Returns:
              ElementaryMatrix: elementary linear thermal matrices
        """

    def getLinearMobilityMatrix(self, groupOfCells=[], with_dual=True):
        """Return the elementary matices for linear mobility acoustic matrix
        Option RIGI_ACOU.

        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
            with_dual (bool): compute dual terms or not (default: True)

        Returns:
            ElementaryMatrix: elementary linear acoustic matrices
        """

    def getMechanicalCouplingForces(self, displ_prev, displ_step, time_prev, time_step):
        """Compute coupling for LIAISON_MASSIF.

        Arguments:
            displ_prev (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step

        Returns:
            FieldOnNodesReal: coupling forces
        """

    def getMechanicalCouplingMatrix(self, displ_prev, displ_step, time_prev, time_step):
        """Compute coupling for LIAISON_MASSIF.

        Arguments:
            displ_prev (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step

        Returns:
            ElementaryMatrixDisplacementReal: coupling elementary matrix.
        """

    def getMechanicalDampingMatrix(
        self,
        getMechanicalMassMatrix=None,
        stiffnessMatrix=None,
        varc_curr=None,
        groupOfCells=[],
        flui_int=1,
        onde_flui=1,
    ):
        """Return the elementary matrices for damping matrix.
        Option AMOR_MECA.

        Arguments:
            getMechanicalMassMatrix : elementary mass matrix
            stiffnessMatrix : elementary stiffness matrix
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
            flui_int (int): integer to activate damping impedance fluid matrix
            onde_flui (int): integer to indicate if we have an outgoing or incoming wave
        Returns:
            ElementaryMatrixReal: elementary damping matrix
        """

    def getMechanicalDirichletBC(self, time_curr=0.0):
        """Return the imposed displacement vector used to remove imposed DDL
        *for internal use - prefer to use getDirichletBC*

        Arguments:
              time_curr (float): Current time (default 0.0)

        Returns:
              FieldOnNodesReal: imposed displacement vector
        """

    def getMechanicalForces(
        self, time_curr=0.0, time_step=0.0, theta=1.0, modeFourier=0, varc_curr=None
    ):
        """Return the total mechanical Neumann forces vector

        Arguments:
              time_curr (float): Current time
              time_step (float): Time increment
              theta (float): Theta parameter for time-integration
              modeFourier (int) : fourier mode
              varc_curr (FieldOnCellsReal): external state variables at current time

        Returns:
              FieldOnNodesReal: forces vector
        """

    def getMechanicalImposedDualBC(self, time_curr=0.0, assembly=True):
        """Return the mechanical imposed nodal BC elementary vector

        Arguments:
              time_curr (float): Current time (default: 0.0)
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorDisplacementReal: imposed dual vector
        """

    def getMechanicalLinearCouplingMatrix(self, varc_curr=None):
        """Compute coupling for LIAISON_MASSIF.

        Arguments:
              varc_curr (FieldOnCellsReal): external state variables for Nitsche method.

        Returns:
              ElementaryMatrixDisplacementReal: coupling elementary matrix.
        """

    def getMechanicalMassMatrix(self, diagonal, varc_curr=None, groupOfCells=[]):
        """Return the elementary matrices for mechanical mass matrix
        Option MASS_MECA.

        Arguments:
            diagonal (bool) : True for diagonal mass matrix else False.
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """

    def getMechanicalNeumannForces(
        self, time_curr=0.0, time_step=0.0, theta=1.0, mode=0, varc_curr=None, assembly=True
    ):
        """Return the elementary mechanical Neumann forces vector

        Arguments:
              time_curr (float): Current time
              time_step (float): Time increment
              theta (float): Theta parameter for time-integration
              mode (int) : fourier mode
              varc_curr (FieldOnCellsReal): external state variables at current time
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorDisplacementReal: elementary Neumann forces vector
        """

    def getMechanicalNodalForces(
        self,
        stress,
        disp=None,
        modeFourier=0,
        varc_curr=None,
        behaviourMap=None,
        groupOfCells=[],
        assembly=True,
    ):
        """Return the elementary mechanical nodal forces vector

        Arguments:
              stress (FieldOnCells): field of stresses
              disp (FieldOnNodes): displacement field (required for large strains hypothesis)
              modeFourier (int) : fourier mode
              varc_curr (FieldOnCellsReal): external state variables
              behaviourMap (FieldOnCellsReal): map for non-linear behaviour
              groupOfCells (list[str]): compute vector on given groups of cells.
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorDisplacementReal: elementary Neumann forces vector
        """

    def getMechanicalReactionForces(
        self,
        disp,
        stress,
        time_prev=0.0,
        time_curr=0.0,
        theta=1.0,
        modeFourier=0,
        varc_curr=None,
        behaviourMap=None,
    ):
        """Return the reaction forces

        Arguments:
              stress (FieldOnCells): field of stresses
              disp (FieldOnNodes): displacement field (required for large strains hypothesis)
              time_prev (float): time at begin of the step
              time_curr (float): time at end of the step
              theta (float): Theta parameter for time-integration
              modeFourier (int) : fourier mode
              varc_curr (FieldOnCellsReal): external state variables at current time
              behaviourMap (FieldOnCellsReal): map for non-linear behaviour

        Returns:
              FieldOnNodesReal: forces vector
        """

    def getMechanicalVolumetricForces(
        self, time_curr=0.0, time_step=0.0, theta=1.0, mode=0, varc_curr=None, assembly=True
    ):
        """Return the elementary mechanical Volumetric forces vector

        Arguments:
              time_curr (float): Current time
              time_step (float): Time increment
              theta (float): Theta parameter for time-integration
              mode (int) : fourier mode
              varc_curr (FieldOnCellsReal): external state variables at current time
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorDisplacementReal: elementary Volumetric forces vector
        """

    def getNonLinearCapacityForces(self, temp_prev, temp_step, varc_curr=None, groupOfCells=[]):
        """Compute internal thermal forces (integration of behaviour)
        Option MASS_THER_RESI.

        Arguments:
            temp_prev (FieldOnNodes): thermal field at begin of current time
            temp_step (FieldOnNodes): field of increment of temperature
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """

    def getPhysicalProblem(self):
        """Get physical probelm

        Returns:
              PhysicalProblem: physical problem
        """

    def getPredictionTangentStiffnessMatrix(
        self,
        displ_prev,
        displ_step,
        stress,
        internVar,
        time_prev,
        time_step,
        varc_prev=None,
        varc_curr=None,
        groupOfCells=[],
    ):
        """Compute jacobian matrix for Newton algorithm, Euler prediction

        Arguments:
            displ_prev (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): internal state variables at begin of current time
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            varc_prev (FieldOnCellsReal): external state variables at begin of current time
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.

        Returns:
            tuple (tuple): return code error (FieldOnCellsLong),
            error code flag (int),
            elementary tangent matrix (ElementaryMatrixDisplacementReal)
        """

    def getResidualReference(self, arg0):
        """Return the residual reference (for RESI_REFE_RELA)

        Arguments:
              vale_by_name : dict :
                             keys are component names
                             values are the given reference value corresponding to component name

        Returns:
              FieldOnNodesReal: residual reference forces vector
        """

    def getRotationalStiffnessMatrix(self, groupOfCells=[]):
        """Return the elementary matrices for rotational Stiffness matrix.
        Option RIGI_ROTA.

        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixReal: elementary rotational rigidity matrix
        """

    def getTangentCapacityMatrix(self, temp_prev, temp_step, varc_curr=None, groupOfCells=[]):
        """Return the elementary matrices for nonlinear Capacity matrix in thermal computation.
        Option MASS_THER_TANG.

        Arguments:
            temp_prev (FieldOnNodes): thermal field at begin of current time
            temp_step (FieldOnNodes): field of increment of temperature
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """

    def getTangentConductivityMatrix(
        self, temp_prev, temp_step, varc_curr=None, groupOfCells=[], with_dual=True
    ):
        """Return the elementary matrices for tangent conductivity.
        Option MASS_THER_TANG.

        Arguments:
            temp_prev (FieldOnNodes): thermal field at begin of current time
            temp_step (FieldOnNodes): field of increment of temperature
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
            with_dual (bool): compute dual terms or not (default: True)
        Returns:
            ElementaryMatrix: elementary mass matrix
        """

    def getTangentStiffnessMatrix(
        self,
        displ_prev,
        displ_step,
        stress,
        internVar,
        internVarIter,
        time_prev,
        time_step,
        varc_prev=None,
        varc_curr=None,
        groupOfCells=[],
    ):
        """Compute jacobian matrix for Newton algorithm

        Arguments:
            displ_prev (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): internal state variables at begin of current time
            internVarIter (FieldOnCells): field of internal state variables
                                          at begin of current newton iteration
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            varc_prev (FieldOnCellsReal): external state variables at begin of current time
            varc_curr (FieldOnCellsReal): external state variables at current time
            groupOfCells (list[str]): compute matrices on given groups of cells.

        Returns:
            tuple (tuple): return code error (FieldOnCellsLong),
            error code flag (int),
            elementary tangent matrix (ElementaryMatrixDisplacementReal)
        """

    def getThermalDirichletBC(self, time_curr=0.0):
        """Return the imposed thermal vector used to remove imposed DDL
        *for internal use - prefer to use getDirichletBC*

        Arguments:
              time_curr (float): Current time (default 0.0)

        Returns:
              FieldOnNodesReal: imposed thermal vector
        """

    def getThermalExchangeForces(self, temp_curr, time_curr=0.0, assembly=True):
        """Return the elementary thermal Exchange forces vector

        Arguments:
              temp_curr (FieldOnNodesReal): thermal field at current time
              time_curr (float): Current time
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorThermalReal: elementary Exchange forces vector
        """

    def getThermalExchangeMatrix(self, time_curr):
        """Return the elementary matices for exhange thermal matrix.

        Arguments:
            time_curr (float): Current time
        Returns:
            ElementaryMatrix: elementary exchange thermal matrices
        """

    def getThermalImposedDualBC(self, time_curr=0.0, assembly=True):
        """Return the thermal imposed nodal BC elementary vector

        Arguments:
              time_curr (float): Current time (default: 0.0)
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorThermalReal: imposed dual vector
        """

    def getThermalNeumannForces(self, time_curr=0.0, assembly=True):
        """Return the elementary thermal Neumann forces vector

        Arguments:
              time_curr (float): Current time (default: 0.0)
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorThermalReal: elementary Neumann forces vector
        """

    def getThermalNonLinearNeumannForces(self, temp_curr, time_curr, assembly=True):
        """Return the elementary field for nonlinear neuamnn forces.
        Option CHAR_THER_FLUTNL, CHAR_THER_RAYO_F, CHAR_THER_RAYO_R.

        Arguments:
            temp_curr (FieldOnNodesReal): thermal field at end of current time
            time_curr (float): Current time
            assembly (bool): assemble or not the field
        Returns:
            ElementaryVector: elementary field
        """

    def getThermalNonLinearVolumetricForces(self, temp_curr, time_curr, assembly=True):
        """Return the elementary field for nonlinear volumetric forces.
        Option CHAR_THER_SOURNL.

        Arguments:
            temp_curr (FieldOnNodesReal): thermal field at end of current time
            time_curr (float): Current time
            assembly (bool): assemble or not the field
        Returns:
            ElementaryVector: elementary field
        """

    def getThermalTangentNonLinearNeumannMatrix(self, temp_curr, time_curr, varc_curr=None):
        """Return the elementary matrices for tangent nonlinear neumann forces.
        Option MTAN_THER_FLUXNL, MTAN_THER_RAYO_R, MTAN_THER_RAYO_F.

        Arguments:
            temp_curr (FieldOnNodesReal): thermal field at end of current time
            time_curr (float): Current time
            varc_curr (FieldOnCellsReal): external state variables at current time

        Returns:
            ElementaryMatrix: elementary matrix
        """

    def getThermalTangentNonLinearVolumetricMatrix(self, temp_curr, time_curr):
        """Return the elementary matrices for tangent nonlinear volumetric forces.
        Option MTAN_THER_SOURNL.

        Arguments:
            temp_curr (FieldOnNodesReal): thermal field at end of current time
            time_curr (float): Current time

        Returns:
            ElementaryMatrix: elementary matrix
        """

    def getThermalVolumetricForces(self, time_curr=0.0, varc_curr=None, assembly=True):
        """Return the elementary thermal Volumetric forces vector

        Arguments:
              time_curr (float): Current time
              varc_curr (FieldOnCellsReal): external state variables at current time
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              ElementaryVectorThermalReal: elementary Volumetric forces vector
        """

    def getTransientThermalForces(self, *args, **kwargs):
        """Overloaded function.

        1. getTransientThermalForces(self: libaster.DiscreteComputation, time_curr: float, time_step: float, theta: float, varc_curr: FieldOnNodes<double> = None, previousPrimalField: FieldOnCells<double> = None) -> FieldOnNodes<double>


                    Compute Transient Thermal Load

                    Arguments:
                          time_curr (float): Current time
                          time_step (float): Time increment
                          theta (float): Theta parameter for integration
                          varc_curr (FieldOnCellsReal): external state variables at current time
                          previousPrimalField (fieldOnNodesReal): solution field at previous time

                    Returns:
                          FieldOnNodes: load


        2. getTransientThermalForces(self: libaster.DiscreteComputation, time_curr: float, time_step: float, theta: float, previousPrimalField: FieldOnNodes<double>, varc_curr: FieldOnCells<double> = None) -> FieldOnNodes<double>


                    Compute Transient Thermal forces due to time scheme
                    Option CHAR_THER_EVOL

                    Arguments:
                          time_curr (float): Current time
                          time_step (float): Time increment
                          theta (float): Theta parameter for integration
                          previousPrimalField (fieldOnNodesReal): solution field at previous time
                          varc_curr (FieldOnCellsReal): external state variables at current time

                    Returns:
                          FieldOnNodes: load
        """

    def getTransientThermalLoadForces(self, time_curr, temp_prev=None, assembly=True):
        """Compute Transient Thermal Load given by EVOL_CHAR.
        Option CHAR_THER.

        Arguments:
              time_curr (float): Current time
              temp_prev (FieldOnNodesReal): solution field at previous time
              assembly (bool) : if True return assembled vector (default: True)

        Returns:
              FieldOnNodes: load
        """


# class EquationNumbering in libaster


class EquationNumbering(DataStructure):
    pass

    # Method resolution order:
    #     EquationNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.EquationNumbering) -> None

        2. __init__(self: libaster.EquationNumbering, arg0: str) -> None
        """

    def getComponents(self):
        """Get list of components

        Returns:
            list[str]: list of components
        """

    def getComponentsIdToName(self):
        """Get map between id and name of components

        Returns:
            dict[int]: map between id and name
        """

    def getComponentsNameToId(self):
        """Get map between id and name of components

        Returns:
            dict[int]: map between id and name
        """

    def getDOFFromNodeAndComponent(self, local=True):
        """Return the dict of dofs with the pair (node id, component's name) as keys

        Arguments:
            local (bool) = True: if True use local dof index else use global index in HPC

        Returns:
            dict[int, str] : dofs id for each node id and component's name
        """

    def getDOFFromNodeAndComponentId(self, local=True):
        """Return the dict of dofs with the pair (node id, name id) as keys

        Arguments:
            local (bool) = True: if True use local DOF index else use global index in HPC

        Returns:
            dict[int, str] : dofs id for each node id and component id
        """

    def getDOFs(self, sameRank=False, list_cmp=[], list_grpno=[]):
        """Return list of DOFs

        Arguments:
            sameRank = False: Use only owned nodes / False: Use all nodes
            list_cmp = []: Use all cmp / keep only cmp given
            list_grpno = []: Use all nodes / keep only nodes given

        Returns:
            list[int]: list of dofs.
        """

    def getDOFsWithDescription(self, *args, **kwargs):
        """Overloaded function.

        1. getDOFsWithDescription(self: libaster.EquationNumbering, cmps: list[str] = [], groupNames: list[str] = [], local: bool = True, same_rank: int = <PythonBool.NONE: -1>) -> tuple[tuple[list[int], list[str]], list[int]]


                    Get the dofs associated to the given component restricted to the given group.

                    Arguments:
                        cmps (list[str]): components to extract.
                        groupNames (list[str]): group names to filter.
                        local (bool): if True use local dof index else use global index in HPC.

                    Returns:
                        pair[list[int], list[str]]: list of nodes and list of components.
                        list[int]: list of dofs.


        2. getDOFsWithDescription(self: libaster.EquationNumbering, cmps: list[str] = [], nodes: list[int] = [], local: bool = True, same_rank: int = <PythonBool.NONE: -1>) -> tuple[tuple[list[int], list[str]], list[int]]


                    Get the dofs associated to the given component restricted to the given nodes.

                    Arguments:
                        cmps (list[str]): components to extract.
                        nodes (list[int]): list of nodes to filter.
                        local (bool): if True use local dof index else use global index in HPC.

                    Returns:
                        pair[list[int], list[str]]: list of nodes and list of components.
                        list[int]: list of dofs.
        """

    def getMesh(self):
        pass

    def getModel(self):
        pass

    def getNoGhostDOFs(self, local=False):
        """Returns the indexes of the DOFs owned locally (aka not ghost).

        Returns:
        int: indexes of the DOFs owned locally.
        """

    def getNodeAndComponentFromDOF(self, *args, **kwargs):
        """Overloaded function.

        1. getNodeAndComponentFromDOF(self: libaster.EquationNumbering, local: bool = True) -> list[tuple[int, str]]


                    Return the list of node id and name of component for each dofs

                    Arguments:
                        local (bool) = True: if True use local node index else use global index in HPC

                    Returns:
                        list[tuple[int, str]] : node id and name of component for each dofs


        2. getNodeAndComponentFromDOF(self: libaster.EquationNumbering, dof: int, local: bool = True) -> tuple[int, str]


                    Return the node id and name of component for given DOF

                    Arguments:
                        dof (int): DOF index
                        local (bool) = True: if True use local node index else use global index in HPC

                    Returns:
                        tuple[int, str] : node id and name of component
        """

    def getNodeAndComponentIdFromDOF(self, *args, **kwargs):
        """Overloaded function.

        1. getNodeAndComponentIdFromDOF(self: libaster.EquationNumbering, local: bool = True) -> list[tuple[int, int]]


                    Return the list of node id and component id for each dofs

                    Arguments:
                        local (bool) = True: if True use local node index else use global index in HPC

                    Returns:
                        list[tuple[int, int]] : node id and component if for each dofs


        2. getNodeAndComponentIdFromDOF(self: libaster.EquationNumbering, dof: int, local: bool = True) -> tuple[int, int]


                    Return the node id and component id for given DOF

                    Arguments:
                        dof (int): DOF index
                        local (bool) = True: if True use local node index else use global index in HPC

                    Returns:
                        tuple[int, int] : node id and component if for each dofs
        """

    def getNumberOfDOFs(self, local=False):
        """Returns the number of DOFs.

        Arguments:
            local (bool): not used.

        Returns:
            int: number of DOFs.
        """

    def getPhysicalQuantity(self):
        """Returns the name of the physical quantity that is numbered.

        Returns:
            str: physical quantity name.
        """

    def isParallel(self):
        """The numbering is distributed across MPI processes for High Performance Computing.

        Returns:
            bool: *True* if used, *False* otherwise.
        """

    def setMesh(self, arg0):
        pass

    def setModel(self, arg0):
        pass

    def useLagrangeDOF(self):
        """Lagrange multipliers are used for BC or MPC.

        Returns:
            bool: *True* if used, *False* otherwise.
        """

    def useSingleLagrangeDOF(self):
        """Single Lagrange multipliers are used for BC or MPC.

        Returns:
            bool: *True* if used, *False* otherwise.
        """


# class MatrixStorage in libaster


class MatrixStorage(DataStructure):
    pass

    # Method resolution order:
    #     MatrixStorage
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0):
        pass


# class MorseStorage in libaster


class MorseStorage(MatrixStorage):
    pass

    # Method resolution order:
    #     MorseStorage
    #     MatrixStorage
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getDiagonalPositions(self):
        pass

    def getRows(self):
        pass


# class BaseDOFNumbering in libaster


class BaseDOFNumbering(DataStructure):
    pass

    # Method resolution order:
    #     BaseDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def computeNumbering(self, *args, **kwargs):
        """Overloaded function.

        1. computeNumbering(self: libaster.BaseDOFNumbering, model: Model, listOfLoads: ListOfLoads, verbose: bool = True) -> bool

        2. computeNumbering(self: libaster.BaseDOFNumbering, matrix: list[Union[ElementaryMatrix<double, (PhysicalQuantityEnum)4>, ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>, ElementaryMatrix<double, (PhysicalQuantityEnum)6>, ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>]], verbose: bool = True) -> bool
        """

    def computeRenumbering(self, model, listOfLoads, defiCont, vContElem, verbose=True):
        pass

    def getEquationNumbering(self):
        """Returns the global equation numbering object

        Returns:
            EquationNumbering: global equation numbering.
        """

    def getFiniteElementDescriptors(self):
        """Returns the objects defining the finite elements.

        Returns:
            list[FiniteElementDescriptor]: List of finite elements descriptions.
        """

    def getMesh(self):
        """Return the mesh

        Returns:
            MeshPtr: a pointer to the mesh
        """

    def getModel(self):
        pass

    def getMorseStorage(self):
        pass

    def getPhysicalQuantity(self):
        """Returns the name of the physical quantity that is numbered.

        Returns:
            str: physical quantity name.
        """

    def isParallel(self):
        """The numbering is distributed across MPI processes for High Performance Computing.

        Returns:
            bool: *True* if used, *False* otherwise.
        """

    def setFiniteElementDescriptors(self, descr):
        """Returns the object defining the finite elements.

        Arguments:
            descr (list[FiniteElementDescriptor]): List of finite elements descriptions.
        """

    def setModel(self, arg0):
        pass


# class DOFNumbering in libaster


class DOFNumbering(BaseDOFNumbering):
    pass

    # Method resolution order:
    #     DOFNumbering
    #     BaseDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.DOFNumbering) -> None

        2. __init__(self: libaster.DOFNumbering, arg0: str) -> None

        3. __init__(self: libaster.DOFNumbering, arg0: str, arg1: libaster.EquationNumbering, arg2: Model) -> None
        """

    def getComponentFromDOF(self, dof, local=False):
        """Returns the component name associated to a dof index.

        - If the dof is associated to a physical DOF, the name of the component is returned.

        - If the dof is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, the name of the component which is constrained by the multiplier is
          returned, precedeed by 'LAGR:', e.g. 'LAGR:DX'.

        - If the dof is associated to a Lagrange multiplier DOF for a multipoint-constraint
          (MPC) implying several DOF, 'LAGR:MPC' is returned (since no component can be
          identified).

        Arguments:
            dof (int): Index of the dof.
            local (bool, optional): not used (default: false).

        Returns:
            str: component name.
        """

    def getComponentFromNode(self, node, local=False):
        """Returns the components name associated to a node index.

        Arguments:
            node (int): Index of the node.
            local (bool, optional): not used (default: false).

        Returns:
            str: component names.
        """

    def getComponents(self):
        """Returns all the component names assigned in the numbering.

        Returns:
            str: component names.
        """

    def getDOFFromNodeAndComponent(self, node, cmp, local=False):
        """Returns the DOF index associated to a node and component.

        Arguments:
            node (int): Index of the node.
            cmp (str): name of the component
            local (bool, optional): not used (default: false).

        Returns:
            int: index of the dof.
        """

    def getDictOfLagrangeDOFs(self, local=False):
        """Returns the Rows Associated to the first and second Lagrange Multipliers Dof

        Arguments:
            local (bool, optional): not used (default: false).

        Returns:
            [dict]: {1 : indexes of the first Lagrange multipliers dof,
                     2 : indexes of the second Lagrange multipliers dof }
        """

    def getLagrangeDOFs(self, local=False):
        """Returns the indexes of the Lagrange multipliers dof.

        Arguments:
            local (bool, optional): not used (default: false).

        Returns:
            [int]: indexes of the Lagrange multipliers dof.
        """

    def getNoGhostDOFs(self, local=False):
        """Returns the indexes of the DOFs owned locally (aka not ghost).

        Returns:
          int: indexes of the DOFs owned locally.
        """

    def getNodeAndComponentFromDOF(self, *args, **kwargs):
        """Overloaded function.

        1. getNodeAndComponentFromDOF(self: libaster.DOFNumbering, local: bool = False) -> list[tuple[int, str]]


        Return the list of node id and name of component for each dofs

        Arguments:
            local (bool, optional): not used (default: false).
        Returns:
            list[tuple[int, str]] : node id and name of component for each dofs


        2. getNodeAndComponentFromDOF(self: libaster.DOFNumbering, dof: int, local: bool = False) -> tuple[int, str]


        Return the node id and name of component for given DOF

        Arguments:
            dof (int): DOF index
            local (bool, optional): not used (default: false).
        Returns:
            tuple[int, str] : node id and name of component
        """

    def getNodeFromDOF(self, dof, local=False):
        """Returns the node index associated to a dof index.

        Arguments:
            dof (int): Index of the dof.
            local (bool, optional): not used (default: false).

        Returns:
            int: index of the node.
        """

    def getNumberOfDOFs(self, local=False):
        """Returns the number of DOFs.

        Arguments:
            local (bool, optional): not used (default: false).

        Returns:
            int: number of DOFs.
        """

    def getPhysicalDOFs(self, local=False):
        """Returns the indexes of the physical dof.

        Arguments:
            local (bool, optional): not used (default: false).

        Returns:
            [int]: indexes of the physical dof.
        """

    def isPhysicalDOF(self, dof, local=False):
        """If the dof is associated to a physical DOF, return True

        If the dof is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, return False

        Arguments:
            dof (int): Index of the dof.
            local (bool, optional): not used (default: false).

        Returns:
            bool: True if the DOF is a physical DOF else False.
        """

    def useLagrangeDOF(self):
        """Lagrange multipliers are used for BC or MPC.

        Returns:
            bool: *True* if used, *False* otherwise.
        """

    def useSingleLagrangeDOF(self):
        """Single Lagrange multipliers are used for BC or MPC.

        Returns:
            bool: *True* if used, *False* otherwise.
        """


# class ElementaryCharacteristics in libaster


class ElementaryCharacteristics(DataStructure):
    pass

    # Method resolution order:
    #     ElementaryCharacteristics
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryCharacteristics, arg0: Model) -> None

        2. __init__(self: libaster.ElementaryCharacteristics, arg0: str, arg1: Model) -> None
        """

    def containsFieldOnCells(self):
        """Return True if ElementaryCharacteristics contains FieldOnCells"""

    def getMesh(self):
        pass

    def getModel(self):
        pass


# class FiniteElementDescriptor in libaster


class FiniteElementDescriptor(DataStructure):
    pass

    # Method resolution order:
    #     FiniteElementDescriptor
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FiniteElementDescriptor, arg0: libaster.BaseMesh) -> None

        2. __init__(self: libaster.FiniteElementDescriptor, arg0: str, arg1: libaster.BaseMesh) -> None

        3. __init__(self: libaster.FiniteElementDescriptor, arg0: libaster.FiniteElementDescriptor, arg1: list[str]) -> None

        4. __init__(self: libaster.FiniteElementDescriptor, model: Model, groupOfCells: list[str]) -> None
        """

    def getListOfGroupsOfElements(self):
        pass

    def getMesh(self):
        pass

    def getNumberOfCells(self):
        pass

    def getPhysics(self):
        pass

    def getVirtualCellsDescriptor(self):
        pass

    def restrict(self, *args, **kwargs):
        """Overloaded function.

        1. restrict(self: libaster.FiniteElementDescriptor, arg0: list[int]) -> libaster.FiniteElementDescriptor

        2. restrict(self: libaster.FiniteElementDescriptor, arg0: list[str]) -> libaster.FiniteElementDescriptor
        """

    def transferDofDescriptorFrom(self, arg0):
        pass


# class FiberGeometry in libaster


class FiberGeometry(DataStructure):
    pass

    # Method resolution order:
    #     FiberGeometry
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FiberGeometry) -> None

        2. __init__(self: libaster.FiberGeometry, arg0: str) -> None
        """


# class DataField in libaster


class DataField(DSWithCppPickling):
    pass

    # Method resolution order:
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.DataField) -> None

        2. __init__(self: libaster.DataField) -> None

        3. __init__(self: libaster.DataField, arg0: str) -> None

        4. __init__(self: libaster.DataField, arg0: str, arg1: str) -> None
        """

    def getFieldType(self):
        """Get field type between "ELEM", "ELGA", "ELNO", "NOEU", "CART".

        Returns:
            str: field type
        """


# class FieldOnCellsReal in libaster


class FieldOnCellsReal(DataField):
    pass

    # Method resolution order:
    #     FieldOnCellsReal
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __iadd__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnCellsReal) -> None

        2. __init__(self: libaster.FieldOnCellsReal, arg0: Model) -> None

        3. __init__(self: libaster.FieldOnCellsReal, arg0: Model, arg1: str, arg2: str) -> None

        4. __init__(self: libaster.FieldOnCellsReal, arg0: libaster.FiniteElementDescriptor, arg1: str, arg2: str) -> None

        5. __init__(self: libaster.FieldOnCellsReal, arg0: libaster.FieldOnCellsReal) -> None

        6. __init__(self: libaster.FieldOnCellsReal, arg0: tuple) -> None

        7. __init__(self: libaster.FieldOnCellsReal, model: Model, loc: str, quantity: str, behaviour: BehaviourProperty, elem_char: libaster.ElementaryCharacteristics) -> None

        8. __init__(self: libaster.FieldOnCellsReal, model: Model, loc: str, quantity: str, behaviour: BehaviourProperty) -> None

        9. __init__(self: libaster.FieldOnCellsReal, model: Model, loc: str, quantity: str, elem_char: libaster.ElementaryCharacteristics) -> None
        """

    def __isub__(self, arg0):
        pass

    def __len__(self):
        pass

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def asLocalization(self, loc):
        """Return a new field interpolated at the given localozation.

        Arguments:
            loc [str]: name of localization (ELEM, ELNO or ELGA)

        Returns:
            FieldOnCellsReal: new field with new localization.
        """

    def build(self, feds=[]):
        pass

    def checkInternalStateVariables(self, prevBehaviour, currBehaviour, newFEDesc):
        """Check consistency of internal states variables with behaviour.
        If you give previous behaviour, check is more precise (name of beahviour for instance)

        Arguments:
            prevBehaviour (ConstantFieldOnCellsChar16): previous behaviour
            currBehaviour (ConstantFieldOnCellsChar16): current behaviour
            newFEDesc (FiniteElementDescriptorPtr): new finite element descriptor
        """

    def compareShape(self, fieldModel, projectOnLigrel, paraName):
        """Compare structure of field with another one and project on new model if require

        Arguments:
            fieldModel (FieldOnCellsRealPtr): field as model
            projectOnLigrel (bool) : project field on new model (from model field)
            paraName (string) : name of parameter to complete the new values in field

        Returns:
            iret (integer) : error code
        """

    def copy(self):
        """Return a duplicated FieldOnCellsReal as a copy

        Returns:
            FieldOnCellsReal
        """

    def dot(self, other):
        """Return the dot product of two fields

        Arguments:
            field (FieldOnCells): other field

        Returns:
            float: dot product
        """

    def exists(self):
        """The field exists or not ?

        Returns:
            Bool
        """

    def getComponents(self):
        """Get list of components

        Returns:
            list[str]: list of components
        """

    def getDescription(self, *args, **kwargs):
        """Overloaded function.

        1. getDescription(self: libaster.FieldOnCellsReal) -> libaster.FiniteElementDescriptor


                    Return the descriptor associated with the FieldOnCellsReal object

                    Returns:
                        FiniteElementDescriptor: FiniteElementDescriptor Object


        2. getDescription(self: libaster.FieldOnCellsReal) -> libaster.FiniteElementDescriptor
        """

    def getLocalization(self):
        """Get localization between ELEM, ELNO and ELGA

        Returns:
            str: localization
        """

    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsReal object

        Returns:
            BaseMesh: Mesh object
        """

    def getNumberOfComponents(self):
        """Get number of components

        Returns:
            int: number of components
        """

    def getPhysicalQuantity(self):
        """Get physical quantity

        Returns:
            str: physical quantity
        """

    def getValues(self, *args, **kwargs):
        """Overloaded function.

        1. getValues(self: libaster.FieldOnCellsReal) -> JeveuxVector


                    Return a list of values as (x1, y1, z1, x2, y2, z2...)

                    Returns:
                        list[float]: List of values.


        2. getValues(self: libaster.FieldOnCellsReal, dofs: list[int] = []) -> list[float]


                    Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

                    Arguments:
                        dofs: dofs to extract

                    Returns:
                        list[float]: List of values.
        """

    def norm(self, arg0):
        """Return the euclidean norm of the field

        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"

        Returns:
            float: euclidean norm
        """

    def printMedFile(self, filename, local=True, version=""):
        """Print the field in MED format.

        Arguments:
            filename (Path|str): Path to the file to be printed.
            local (bool): Print local values only (relevant for ParallelMesh only,
                default: *True*)
            version (str): Version of MED file.

        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """

    def setDescription(self, arg0):
        pass

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.FieldOnCellsReal, value: float) -> None


                    Set values of the field

                    Arguments:
                        value (float): value to set


        2. setValues(self: libaster.FieldOnCellsReal, values: list[float]) -> None


                    Set values of the field

                    Arguments:
                        values (list[float]): list of values to set
        """

    def size(self):
        """Return the size of the field

        Returns:
            int: number of element in the field
        """

    def toFieldOnNodes(self):
        """Convert to FieldOnNodes

        Returns:
            FieldOnNodesReal: field converted
        """

    def toSimpleFieldOnCells(self):
        """Convert to SimpleFieldOnNodes

        Returns:
            SimpleFieldOnNodesReal: field converted
        """

    def toSimpleFieldOnNodes(self):
        """Convert to SimpleFieldOnNodes

        Returns:
            SimpleFieldOnNodesReal: field converted
        """

    def transform(self, func):
        """Apply a function to each value of the object.

        Arguments:
            func (*callable*): Callable Python object

        Returns:
            FieldOnCellsReal: New FieldOnCells object with the transformed values
        """


# class FieldOnCellsComplex in libaster


class FieldOnCellsComplex(DataField):
    pass

    # Method resolution order:
    #     FieldOnCellsComplex
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __iadd__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnCellsComplex) -> None

        2. __init__(self: libaster.FieldOnCellsComplex, arg0: libaster.FiniteElementDescriptor, arg1: str, arg2: str) -> None

        3. __init__(self: libaster.FieldOnCellsComplex, arg0: tuple) -> None

        4. __init__(self: libaster.FieldOnCellsComplex, arg0: libaster.FieldOnCellsComplex) -> None
        """

    def __isub__(self, arg0):
        pass

    def __len__(self):
        pass

    def __mul__(self, arg0):
        pass

    def __rmul__(self, arg0):
        pass

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def build(self, feds=[]):
        pass

    def copy(self):
        pass

    def getDescription(self):
        pass

    def getLocalization(self):
        """Get localization between ELEM, ELNO and ELGA

        Returns:
            str: localization
        """

    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsReal object

        Returns:
            BaseMesh: Mesh object
        """

    def getPhysicalQuantity(self):
        """Get physical quantity

        Returns:
            str: physical quantity
        """

    def getValues(self, *args, **kwargs):
        """Overloaded function.

        1. getValues(self: libaster.FieldOnCellsComplex) -> JeveuxVector


                    Return a list of values as (x1, y1, z1, x2, y2, z2...)

                    Returns:
                        list[complex]: List of values.


        2. getValues(self: libaster.FieldOnCellsComplex, dofs: list[int] = []) -> list[complex]


                    Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

                    Arguments:
                        dofs: dofs to extract

                    Returns:
                        list[complex]: List of values.
        """

    def printMedFile(self, filename, local=True, version=""):
        """Print the field in MED format.

        Arguments:
            filename (Path|str): Path to the file to be printed.
            local (bool): Print local values only (relevant for ParallelMesh only,
                default: *True*)
            version (str): Version of MED file.

        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """

    def setDescription(self, arg0):
        pass

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.FieldOnCellsComplex, value: complex) -> None


                    Set values of the field

                    Arguments:
                        value (complex): value to set


        2. setValues(self: libaster.FieldOnCellsComplex, values: list[complex]) -> None


                    Set values of the field

                    Arguments:
                        values (list[complex]): list of values to set
        """

    def size(self):
        """Return the size of the field

        Returns:
            int: number of element in the field
        """

    def toFieldOnNodes(self):
        """Convert to FieldOnNodes

        Returns:
            FieldOnCellsComplex: field converted
        """

    def transform(self, func):
        """Apply a function to each value of the object.

        Arguments:
            func (*callable*): Callable Python object

        Returns:
            FieldOnCellsComplex: New FieldOnCells object with the transformed values
        """


# class FieldOnCellsLong in libaster


class FieldOnCellsLong(DataField):
    pass

    # Method resolution order:
    #     FieldOnCellsLong
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __iadd__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnCellsLong) -> None

        2. __init__(self: libaster.FieldOnCellsLong, arg0: tuple) -> None

        3. __init__(self: libaster.FieldOnCellsLong, arg0: libaster.FieldOnCellsLong) -> None
        """

    def __isub__(self, arg0):
        pass

    def __len__(self):
        pass

    def __mul__(self, arg0):
        pass

    def __rmul__(self, arg0):
        pass

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def build(self, feds=[]):
        pass

    def copy(self):
        pass

    def getDescription(self):
        pass

    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsReal object

        Returns:
            BaseMesh: Mesh object
        """

    def getValues(self, *args, **kwargs):
        """Overloaded function.

        1. getValues(self: libaster.FieldOnCellsLong) -> JeveuxVector


                    Return a list of values as (x1, y1, z1, x2, y2, z2...)

                    Returns:
                        list[int]: List of values.


        2. getValues(self: libaster.FieldOnCellsLong, dofs: list[int] = []) -> list[int]


                    Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

                    Arguments:
                        dofs: dofs to extract

                    Returns:
                        list[int]: List of values.
        """

    def printMedFile(self, filename, local=True, version=""):
        """Print the field in MED format.

        Arguments:
            filename (Path|str): Path to the file to be printed.
            local (bool): Print local values only (relevant for ParallelMesh only,
                default: *True*)
            version (str): Version of MED file.

        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """

    def setDescription(self, arg0):
        pass

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.FieldOnCellsLong, value: int) -> None


                    Set values of the field

                    Arguments:
                        value (complex): value to set


        2. setValues(self: libaster.FieldOnCellsLong, values: list[int]) -> None


                    Set values of the field

                    Arguments:
                        values (list[complex]): list of values to set
        """

    def size(self):
        """Return the size of the field

        Returns:
            int: number of element in the field
        """


# class FieldOnCellsChar8 in libaster


class FieldOnCellsChar8(DataField):
    pass

    # Method resolution order:
    #     FieldOnCellsChar8
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnCellsChar8) -> None

        2. __init__(self: libaster.FieldOnCellsChar8, arg0: tuple) -> None

        3. __init__(self: libaster.FieldOnCellsChar8, arg0: libaster.FieldOnCellsChar8) -> None
        """

    def __setstate__(self, arg0):
        pass

    def build(self, feds=[]):
        pass

    def getDescription(self):
        """Return the description associated with the FieldOnCellsChar8 object

        Returns:
            FiniteElementDescriptor: FiniteElementDescriptor Object
        """

    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsChar8 object

        Returns:
            BaseMesh: Mesh object
        """

    def setDescription(self, arg0):
        pass


# class FieldOnNodesReal in libaster


class FieldOnNodesReal(DataField):
    pass

    # Method resolution order:
    #     FieldOnNodesReal
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __iadd__(self, arg0):
        pass

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnNodesReal) -> None

        2. __init__(self: libaster.FieldOnNodesReal, arg0: libaster.FieldOnNodesReal) -> None

        3. __init__(self: libaster.FieldOnNodesReal, arg0: Model) -> None

        4. __init__(self: libaster.FieldOnNodesReal, arg0: libaster.BaseDOFNumbering) -> None

        5. __init__(self: libaster.FieldOnNodesReal, arg0: libaster.EquationNumbering) -> None

        6. __init__(self: libaster.FieldOnNodesReal, arg0: tuple) -> None

        7. __init__(self: libaster.FieldOnNodesReal, mesh: libaster.BaseMesh, quantity: str, cmps: list[str]) -> None

        8. __init__(self: libaster.FieldOnNodesReal, mesh: libaster.BaseMesh, quantity: str, values: dict[str, float], groupsOfNodes: list[str] = [], groupsOfCells: list[str] = []) -> None
        """

    def __isub__(self, arg0):
        pass

    def __itruediv__(self, *args, **kwargs):
        """Overloaded function.

        1. __itruediv__(self: libaster.FieldOnNodesReal, arg0: float) -> libaster.FieldOnNodesReal

        2. __itruediv__(self: libaster.FieldOnNodesReal, arg0: libaster.FieldOnNodesReal) -> libaster.FieldOnNodesReal
        """

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def __truediv__(self, arg0):
        pass

    def applyLagrangeScaling(self, scaling):
        """Multiply in-place the Lagrange multipliers DOFs by the scaling value

        Arguments:
            scaling (float): scaling velue
        """

    def asPhysicalQuantity(self, physQuantity, map_cmps):
        """Return a new field with a new physical quantity and renamed components.

        Arguments:
            physQuantity [str]: name of the new physical quantity
            map_cmps[dict[str, str]]: dict to rename components
            (only renamed component will be keeped)

        Returns:
            FieldOnNodesReal: field with name physical quantity.
        """

    def build(self, mesh=None):
        pass

    def copy(self):
        pass

    def copyUsingDescription(self, desc, warn=True):
        """Return a new field using the description.
        Be careful, Lagrange DOFs are set to zero. Moreover, components that are
        not present in the field are also set to zero in the output field.

        Arguments:
            desc [EquationNumbering]: description of equations
            warn [bool]: If set to true, raises a warning if values are set to zero

        Returns:
            FieldOnNodesReal: field using new description.
        """

    def dot(self, other):
        """Return the dot product of two fields

        Arguments:
            other (FieldOnNodes): other field

        Returns:
            float: dot product
        """

    def fromPetsc(self, vec, scaling=1.0, local=False):
        """Import a PETSc vector into the field.

        Arguments:
            vec (Vec): The PETSc vector
            scaling (float) : The scaling of the Lagrange DOFs
            local (bool) : Only import the dof that are local to the subdomain
        """

    def getComponents(self):
        """Get list of components

        Returns:
            list[str]: list of components
        """

    def getDescription(self):
        pass

    def getImaginaryPart(self):
        """Extract the imaginary part of the real field (a 0-filled field is produced)

        Returns:
            FieldOnNodesReal: imaginary part
        """

    def getLocalization(self):
        """Get localization = NOEU

        Returns:
            str: "NOEU"
        """

    def getMesh(self, *args, **kwargs):
        """Overloaded function.

        1. getMesh(self: libaster.FieldOnNodesReal) -> libaster.BaseMesh

        2. getMesh(self: libaster.FieldOnNodesReal) -> libaster.BaseMesh
        """

    def getNumberOfComponents(self):
        """Get number of components

        Returns:
            int: number of components
        """

    def getPhysicalQuantity(self):
        pass

    def getRealPart(self):
        """Extract the real part of the real field (the field is duplicated)

        Returns:
            FieldOnNodesReal: real part
        """

    def getValues(self, *args, **kwargs):
        """Overloaded function.

        1. getValues(self: libaster.FieldOnNodesReal) -> JeveuxVector


                    Return a list of values as (x1, y1, z1, x2, y2, z2...)

                    Returns:
                        list[float]: List of values.


        2. getValues(self: libaster.FieldOnNodesReal, cmps: list[str] = [], groupsOfNodes: list[str] = []) -> list[float]


                    Return a list of values as (x1, y1, z1, x2, y2, z2...)

                    Arguments:
                        cmps[list[str]]: filter on list of components
                        groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
                        If empty, the full mesh is used

                    Returns:
                        list[double]: List of values.


        3. getValues(self: libaster.FieldOnNodesReal, dofs: list[int] = []) -> list[float]


                    Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

                    Arguments:
                        dofs: dofs to extract

                    Returns:
                        list[double]: List of values.
        """

    def norm(self, normType="NORM_INFINITY", list_cmp=[]):
        """Return the euclidean norm of the field

        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY" (default: "NORM_INFINITY")
            list_cmp (list[str]) : list of components used to compute norm (default: all)

        Returns:
            float: euclidean norm
        """

    def printMedFile(self, fileName, local=True, version=""):
        pass

    def scale(self, vect):
        """Scale in-place the field by a diagonal matrix stored as an array

        Arguments:
            vect (float): diagonal matrix stored as an array
        """

    def setDescription(self, arg0):
        pass

    def setMesh(self, arg0):
        pass

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.FieldOnNodesReal, value: float) -> None


                    Set values of the field

                    Arguments:
                        value (float): value to set


        2. setValues(self: libaster.FieldOnNodesReal, values: list[float]) -> None


                    Set values of the field

                    Arguments:
                        values (list[float]): list of values to set


        3. setValues(self: libaster.FieldOnNodesReal, value: dict[str, float], groupsOfNodes: list[str] = []) -> None


                    Set values of the field where components and values are given as a dict.
                    If the component is not present in the field then it is discarded
                    Example: { "X1" : 0.0, "X3" : 0.0 }

                    Arguments:
                        value (dict[str, float]): dict of values to set (key: str, value: float)
                        groupsOfNodes (list[str]): list of groups. If empty, the full mesh is considered
        """

    def size(self):
        """Return the size of the field

        Returns:
            int: number of element in the field
        """

    def toFieldOnCells(self, fed, loc):
        """Converts to FieldOnCells

        Arguments:
            fed [FiniteElementDescriptor]: finite element descriptor
            loc [str] : name of localization like 'ELGA'.

        Returns:
            FieldOnCellsReal: field converted.
        """

    def toSimpleFieldOnNodes(self):
        """Convert to SimpleFieldOnNodes

        Returns:
            SimpleFieldOnNodesReal: field converted
        """

    def transferFromConnectionToParallelMesh(self, arg0):
        """Transfer FieldOnNodes from a ConnectionMesh to a ParallelMesh

        Arguments:
            mesh [Mesh]: the target mesh

        Returns:
            FieldOnNodesReal: transfered field
        """

    def transfertToConnectionMesh(self, arg0):
        """Transfer SimpleFieldOnNodes to a ConnectionMesh

        Returns:
            FieldOnNodesReal: transfered field
        """

    def transform(self, func):
        """Apply a function to each value of the object.

        Arguments:
            func (*callable*): Callable Python object

        Returns:
            FieldOnNodesReal: New FieldOnNodes object with the transformed values
        """

    def updateGhostValues(self):
        """Communicates the values of the ghost DOFs on a FieldOnNodes."""

    def updateValuePointers(self):
        pass


# class FieldOnNodesComplex in libaster


class FieldOnNodesComplex(DataField):
    pass

    # Method resolution order:
    #     FieldOnNodesComplex
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnNodesComplex) -> None

        2. __init__(self: libaster.FieldOnNodesComplex, arg0: libaster.FieldOnNodesComplex) -> None

        3. __init__(self: libaster.FieldOnNodesComplex, arg0: Model) -> None

        4. __init__(self: libaster.FieldOnNodesComplex, arg0: libaster.BaseDOFNumbering) -> None

        5. __init__(self: libaster.FieldOnNodesComplex, arg0: libaster.EquationNumbering) -> None

        6. __init__(self: libaster.FieldOnNodesComplex, arg0: tuple) -> None
        """

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def build(self, mesh=None):
        pass

    def dot(self, other):
        """Return the dot product of two complex fields

        Arguments:
            other (FieldOnNodes): other field

        Returns:
            complex: dot product
        """

    def getComponents(self):
        """Get list of components

        Returns:
            list[str]: list of components
        """

    def getDescription(self):
        pass

    def getImaginaryPart(self):
        """Extract the imaginary part of the complex field

        Returns:
            FieldOnNodesReal: imaginary part
        """

    def getLocalization(self):
        """Get localization = NOEU

        Returns:
            str: "NOEU"
        """

    def getMesh(self, *args, **kwargs):
        """Overloaded function.

        1. getMesh(self: libaster.FieldOnNodesComplex) -> libaster.BaseMesh

        2. getMesh(self: libaster.FieldOnNodesComplex) -> libaster.BaseMesh
        """

    def getNumberOfComponents(self):
        """Get number of components

        Returns:
            int: number of components
        """

    def getPhysicalQuantity(self):
        pass

    def getRealPart(self):
        """Extract the real part of the complex field

        Returns:
            FieldOnNodesReal: real part
        """

    def getValues(self, *args, **kwargs):
        """Overloaded function.

        1. getValues(self: libaster.FieldOnNodesComplex) -> JeveuxVector


                    Return a list of values as (x1, y1, z1, x2, y2, z2...)

                    Returns:
                        list[complex]: List of values.


        2. getValues(self: libaster.FieldOnNodesComplex, cmps: list[str] = [], groupsOfNodes: list[str] = []) -> list[complex]


                    Return a list of values as (x1, y1, z1, x2, y2, z2...)

                    Arguments:
                        cmps[list[str]]: filter on list of components
                        groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
                        If empty, the full mesh is used

                    Returns:
                        list[complex]: List of values.


        3. getValues(self: libaster.FieldOnNodesComplex, dofs: list[int] = []) -> list[complex]


                    Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

                    Arguments:
                        dofs: dofs to extract

                    Returns:
                        list[complex]: List of values.
        """

    def norm(self, normType="NORM_INFINITY", list_cmp=[]):
        """Return the euclidean norm of the field

        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY" (default: "NORM_INFINITY")
            list_cmp (list[str]) : list of components used to compute norm (default: all)

        Returns:
            float: euclidean norm
        """

    def printMedFile(self, fileName, local=True, version=""):
        pass

    def scale(self, vect):
        """Scale in-place the field by a diagonal matrix stored as an array

        Arguments:
            vect (float): diagonal matrix stored as an array
        """

    def setDescription(self, arg0):
        pass

    def setMesh(self, arg0):
        pass

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.FieldOnNodesComplex, value: complex) -> None


                    Set values of the field

                    Arguments:
                        value (complex): value to set


        2. setValues(self: libaster.FieldOnNodesComplex, values: list[complex]) -> None


                    Set values of the field

                    Arguments:
                        values (list[complex]): list of values to set
        """

    def toSimpleFieldOnNodes(self):
        """Convert to SimpleFieldOnNodes

        Returns:
            SimpleFieldOnNodesComplex: field converted
        """

    def transform(self, func):
        """Apply a function to each value of the object.

        Arguments:
            func (*callable*): Callable Python object

        Returns:
            FieldOnNodesComplex: New FieldOnNodes object with the transformed values
        """

    def updateValuePointers(self):
        pass


# class FieldOnNodesLong in libaster


class FieldOnNodesLong(DataField):
    pass

    # Method resolution order:
    #     FieldOnNodesLong
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnNodesLong) -> None

        2. __init__(self: libaster.FieldOnNodesLong, arg0: libaster.FieldOnNodesLong) -> None

        3. __init__(self: libaster.FieldOnNodesLong, arg0: libaster.BaseDOFNumbering) -> None

        4. __init__(self: libaster.FieldOnNodesLong, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def build(self, mesh=None):
        pass

    def getDescription(self):
        pass

    def getMesh(self):
        pass

    def setDescription(self, arg0):
        pass

    def setMesh(self, arg0):
        pass


# class FieldOnNodesChar8 in libaster


class FieldOnNodesChar8(DataField):
    pass

    # Method resolution order:
    #     FieldOnNodesChar8
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FieldOnNodesChar8) -> None

        2. __init__(self: libaster.FieldOnNodesChar8, arg0: libaster.FieldOnNodesChar8) -> None

        3. __init__(self: libaster.FieldOnNodesChar8, arg0: libaster.BaseDOFNumbering) -> None

        4. __init__(self: libaster.FieldOnNodesChar8, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def build(self, mesh=None):
        pass

    def getDescription(self):
        pass

    def getMesh(self):
        pass

    def setDescription(self, arg0):
        pass

    def setMesh(self, arg0):
        pass


# class ConstantFieldValuesReal in libaster


class ConstantFieldValuesReal:
    pass

    # Method resolution order:
    #     ConstantFieldValuesReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getValues(self):
        """Return the field values

        Returns:
            list[float]: List of values
        """


# class ConstantFieldOnCellsReal in libaster


class ConstantFieldOnCellsReal(DataField):
    pass

    # Method resolution order:
    #     ConstantFieldOnCellsReal
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ConstantFieldOnCellsReal, arg0: libaster.BaseMesh) -> None

        2. __init__(self: libaster.ConstantFieldOnCellsReal, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getMesh(self):
        pass

    def getValues(self, arg0):
        """Return the field values

        Returns:
            list[float]: List of values
        """

    def setValueOnCells(self, arg0, arg1, arg2):
        pass

    def size(self):
        """Return the size of field

        Returns:
            int: size of field
        """

    def toSimpleFieldOnCells(self, arg0):
        """Convert to SimpleFieldOnCells

        Returns:
            SimpleFieldOnCellsReal: field converted
        """


# class ConstantFieldOnCellsChar16 in libaster


class ConstantFieldOnCellsChar16(DataField):
    pass

    # Method resolution order:
    #     ConstantFieldOnCellsChar16
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ConstantFieldOnCellsChar16, arg0: libaster.BaseMesh) -> None

        2. __init__(self: libaster.ConstantFieldOnCellsChar16, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getMesh(self):
        pass


# class ConstantFieldOnCellsLong in libaster


class ConstantFieldOnCellsLong(DataField):
    pass

    # Method resolution order:
    #     ConstantFieldOnCellsLong
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ConstantFieldOnCellsLong, arg0: libaster.BaseMesh) -> None

        2. __init__(self: libaster.ConstantFieldOnCellsLong, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getMesh(self):
        pass


# class SimpleFieldOnCellsReal in libaster


class SimpleFieldOnCellsReal(DataField):
    pass

    # Method resolution order:
    #     SimpleFieldOnCellsReal
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.SimpleFieldOnCellsReal, mesh: libaster.BaseMesh) -> None

        2. __init__(self: libaster.SimpleFieldOnCellsReal, mesh: libaster.BaseMesh, loc: str, quantity: str, cmps: list[str]) -> None

        3. __init__(self: libaster.SimpleFieldOnCellsReal, mesh: libaster.BaseMesh, loc: str, quantity: str, cmps: list[str], prol_zero: bool) -> None

        4. __init__(self: libaster.SimpleFieldOnCellsReal, mesh: libaster.BaseMesh, loc: str, quantity: str, cmps: list[str], nbPoints: int, nbSubPoints: int) -> None

        5. __init__(self: libaster.SimpleFieldOnCellsReal, mesh: libaster.BaseMesh, loc: str, quantity: str, cmps: list[str], nbPoints: int, nbSubPoints: int, prol_zero: bool) -> None

        6. __init__(self: libaster.SimpleFieldOnCellsReal, mesh: libaster.BaseMesh, loc: str, quantity: str, cmps: list[str], nbPoints: list[int], nbSubPoints: int, prol_zero: bool) -> None

        7. __init__(self: libaster.SimpleFieldOnCellsReal, arg0: tuple) -> None
        """

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def allocate(self, loc, quantity, cmps, nbPG, nbSP=1, zero=False):
        """Allocate the field.

        Arguments:
            loc [str]: localization like 'ELEM'
            quantity [str]: physical quantity like 'DEPL_R'
            cmps [list[str]]: list of components.
            nbPG [int]: number of Gauss Point by cell
            nbSP [int]: number of sub-point by point.
        """

    def asPhysicalQuantity(self, physQuantity, map_cmps):
        """Return a new field with a new physical quantity and renamed components.

        Arguments:
            physQuantity [str]: name of the new physical quantity
            map_cmps [dict[str, str]]: dict to rename components
            (only renamed component will be keeped)

        Returns:
            SimpleFieldOnCellsReal: field with name physical quantity.
        """

    def getCellsWithValues(self):
        """Returns the list of cells where the field is defined.

        Returns:
            tuple (int): Indexes of cells where the field is defined.
        """

    def getComponent(self, arg0):
        pass

    def getComponents(self):
        pass

    def getComponentsName2Index(self):
        pass

    def getLocalization(self):
        pass

    def getMaxNumberOfPoints(self):
        pass

    def getMaxNumberOfSubPoints(self):
        pass

    def getMesh(self):
        """Returns base mesh"""

    def getNumberOfCells(self):
        pass

    def getNumberOfComponents(self):
        pass

    def getNumberOfComponentsForSubpointsOfCell(self, arg0):
        pass

    def getNumberOfPointsOfCell(self, arg0):
        pass

    def getNumberOfSubPointsOfCell(self, arg0):
        pass

    def getPhysicalQuantity(self):
        pass

    def getValue(self, ima, icmp, ipt, ispt=0):
        """Returns the value of the `icmp` component of the field on the `ima` cell,
        at the `ipt` point, at the `ispt` sub-point.

        Args:
            ima  (int): Index of cells.
            icmp (int): Index of component.
            ipt  (int): Index of point.
            ispt (int): Index of sub-point (default = 0).

        Returns:
            float: Value of field at *ima*, of *icmp*, at *ipt*, at *ispt*;
                     NaN if the position is not allocated.
        """

    def getValuesWithDescription(self, *args, **kwargs):
        """Overloaded function.

        1. getValuesWithDescription(self: libaster.SimpleFieldOnCellsReal, cmps: list[str], cells: list[int]) -> tuple[list[float], tuple[list[int], list[str], list[int], list[int]]]


        Returns values and description corresponding to given cmp and given cells

        Arguments:
            cmps (list[str]): components to extract.
            cells (list[int]): list of cells.

        Returns:
            values[list[double],
            tuple[cells[list[int]], cmps[list[int]]  points[list[int]], subpoints[list[int]]]


        2. getValuesWithDescription(self: libaster.SimpleFieldOnCellsReal, cmps: list[str] = [], cells: list[str] = []) -> tuple[list[float], tuple[list[int], list[str], list[int], list[int]]]


        Returns values and description corresponding to given cmp and given cells

        Arguments:
            cmps (list[str]): components to extract.
            groupsOfCells (list[str]): list of groups of cells to use.

        Returns:
            values[list[double],
            tuple[cells[list[int]], cmps[list[int]]  points[list[int]], subpoints[list[int]]]
        """

    def hasValue(self, *args, **kwargs):
        """Overloaded function.

        1. hasValue(self: libaster.SimpleFieldOnCellsReal, ima: int, icmp: int, ipt: int, ispt: int = 0) -> bool


        Returns True  if the value of the `icmp` component of the field on the `ima` cell,
        at the `ipt` point, at the `ispt` sub-point is affected.

        Args:
            ima  (int): Index of cells.
            icmp (int): Index of component.
            ipt  (int): Index of point.
            ispt (int): Index of sub-point (default = 0).

        Returns:
            bool: True  if the value is affected


        2. hasValue(self: libaster.SimpleFieldOnCellsReal, ima: int, cmp: str, ipt: int, ispt: int = 0) -> bool


        Returns True  if the value of the `icmp` component of the field on the `ima` cell,
        at the `ipt` point, at the `ispt` sub-point is affected.

        Args:
            ima  (int): Index of cells.
            cmp (str): name of component.
            ipt  (int): Index of point.
            ispt (int): Index of sub-point (default = 0).

        Returns:
            bool: True  if the value is affected
        """

    def norm(self, arg0):
        """Return the norm of the field

        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"

        Returns:
            float: norm
        """

    def restrict(self, cmps=[], groupsOfCells=[]):
        """Return a new field restricted to the list of components and groups of cells given

        Arguments:
            cmps[list[str]]: filter on list of components
            If empty, all components are used
            groupsOfCells[list[str]]: filter on list of groups of cells (default=" ").
            If empty, the full mesh is used

        Returns:
            SimpleFieldOnCellsReal: field restricted.
        """

    def setValue(self, *args, **kwargs):
        """Overloaded function.

        1. setValue(self: libaster.SimpleFieldOnCellsReal, ima: int, icmp: int, ipt: int, ispt: int, val: float) -> None


        Set the value of the `icmp` component of the field on the `ima` cell,
        at the `ipt` point, at the `ispt` sub-point.

        Args:
            ima  (int): Index of cells.
            icmp (int): Index of component.
            ipt  (int): Index of point.
            ispt (int): Index of sub-point.
            val (float) : value to set


        2. setValue(self: libaster.SimpleFieldOnCellsReal, ima: int, icmp: int, ipt: int, val: float) -> None


        Set the value of the `icmp` component of the field on the `ima` cell,
        at the `ipt` point, at the `ispt=0` sub-point.

        Args:
            ima  (int): Index of cells.
            icmp (int): Index of component.
            ipt  (int): Index of point.
            val (float) : value to set
        """

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.SimpleFieldOnCellsReal, cells: list[int], cmps: list[str], npg: list[int], spt: list[int], values: list[float]) -> None


                    Set values for a given list of tuple (cell, cmp, ipg, isp, value).
                    Each value of the tuple is given as a separated list.

                    Arguments:
                        cells (list[int]): list of cells.
                        cmps (list[str)]: list of components
                        npg (list[int]): list of point
                        spt (list[int]): list of sub-point
                        values (list[float]): list of values to set.


        2. setValues(self: libaster.SimpleFieldOnCellsReal, values: list[float]) -> None


                     Set values for each cells and components as (cell_0_val_0, cell_0_val_1, ...)

                    Arguments:
                        values (list[float]): list of values to set.
        """

    def toFieldOnCells(self, fed, option="", nompar="", zeroExtension=True):
        """Converts to FieldOnCells

        Arguments:
            fed [FiniteElementDescriptor]: finite element descriptor
            option [str]: name of option like TOUT_INI_ELGA (default: " ")
            nompar [str]: name of parameter like DEPL_R (default: " ")
            zeroExtension [bool]: true if field could be extended to zero (default: True)

        Returns:
            FieldOnCellsReal: field converted.
        """

    def toFieldOnNodes(self):
        """Convert to FieldOnNodes

        Returns:
            FieldOnNodesReal: field converted
        """

    def toNumpy(self):
        """Returns two numpy arrays with shape ( number_of_cells_with_components, number_of_components )
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.

        Where the mask is `False` the corresponding value is set to zero.

        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """

    def toSimpleFieldOnNodes(self):
        """Convert to SimpleFieldOnNodes

        Returns:
            SimpleFieldOnNodesReal: field converted
        """

    def updateValuePointers(self):
        pass


# class SimpleFieldOnNodesReal in libaster


class SimpleFieldOnNodesReal(DataField):
    pass

    # Method resolution order:
    #     SimpleFieldOnNodesReal
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getitem__(self, *args, **kwargs):
        """Overloaded function.

        1. __getitem__(self: libaster.SimpleFieldOnNodesReal, arg0: tuple[int, int]) -> float

        2. __getitem__(self: libaster.SimpleFieldOnNodesReal, arg0: tuple[int, str]) -> float
        """

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.SimpleFieldOnNodesReal, mesh: libaster.BaseMesh) -> None

        2. __init__(self: libaster.SimpleFieldOnNodesReal, mesh: libaster.BaseMesh, quantity: str, cmps: list[str]) -> None

        3. __init__(self: libaster.SimpleFieldOnNodesReal, mesh: libaster.BaseMesh, quantity: str, cmps: list[str], prol_zero: bool) -> None

        4. __init__(self: libaster.SimpleFieldOnNodesReal, arg0: tuple) -> None
        """

    def __setitem__(self, *args, **kwargs):
        """Overloaded function.

        1. __setitem__(self: libaster.SimpleFieldOnNodesReal, arg0: tuple[int, int], arg1: float) -> float

        2. __setitem__(self: libaster.SimpleFieldOnNodesReal, arg0: tuple[int, str]) -> float
        """

    def __setstate__(self, arg0):
        pass

    def allocate(self, quantity, cmps, zero=False):
        """Allocate the field.

        Arguments:
            quantity [str]: physical quantity like 'DEPL_R'
            cmps [list[str]]: list of components.
        """

    def asPhysicalQuantity(self, physQuantity, map_cmps):
        """Return a new field with a new physical quantity and renamed components.

        Arguments:
            physQuantity [str]: name of the new physical quantity
            map_cmps [dict[str, str]]: dict to rename components
            (only renamed component will be keeped)

        Returns:
            SimpleFieldOnNodesReal: field with name physical quantity.
        """

    def getComponent(self, arg0):
        pass

    def getComponents(self):
        pass

    def getLocalization(self):
        pass

    def getMesh(self):
        """Returns base mesh"""

    def getNumberOfComponents(self):
        pass

    def getNumberOfNodes(self):
        pass

    def getPhysicalQuantity(self):
        pass

    def getValuesWithDescription(self, *args, **kwargs):
        """Overloaded function.

        1. getValuesWithDescription(self: libaster.SimpleFieldOnNodesReal, cmps: list[str] = [], groupsOfNodes: list[str] = []) -> tuple[list[float], tuple[list[int], list[str]]]


                    Return the values of components of the field.

                    Arguments:
                       cmps (list[str]) : Extracted components or all components if it is empty.
                       groups (list[str]): The extraction is limited to the given groups of nodes.

                    Returns:
                       tuple( values, description ): List of values and description.
                        The description provides a tuple with( nodes ids, components ).


        2. getValuesWithDescription(self: libaster.SimpleFieldOnNodesReal, cmps: list[str], nodes: list[int]) -> tuple[list[float], tuple[list[int], list[str]]]


                    Return the values of components of the field.

                    Arguments:
                       cmps (list[str]) : Extracted components or all components if it is empty.
                       nodes (list[int]): The extraction is limited to the given nodes.

                    Returns:
                       tuple( values, description ): List of values and description.
                        The description provides a tuple with( nodes ids, components ).
        """

    def hasComponent(self, arg0):
        pass

    def norm(self, arg0):
        """Return the norm of the field

        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"

        Returns:
            float: norm
        """

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.SimpleFieldOnNodesReal, nodes: list[int], cmps: list[str], values: list[float]) -> None


                    Set values for a given list of triplet (node, cmp, value).
                    Each value of the triplet is given as a separated list.

                    Arguments:
                        nodes (list[int]): list of nodes.
                        cmps (list[str]): list of comp components
                        values (list[float]): list of values to set.


        2. setValues(self: libaster.SimpleFieldOnNodesReal, values: list[float]) -> None


                     Set values for each nodes and components as (node_0_val_0, node_0_val_1, ...)

                    Arguments:
                        values (list[float]): list of values to set.



        3. setValues(self: libaster.SimpleFieldOnNodesReal, values: list[list[float]]) -> None


                    Set values for each nodes and components.

                    Arguments:
                        values (list[list[float]]): list of values to set.
                        For each node, give the values for all component is a list.


        4. setValues(self: libaster.SimpleFieldOnNodesReal, value: dict[str, float], nodes: list[int]) -> None


                    Set values of the field where components and values are given as a dict.
                    If the component is not present in the field then it is discarded
                    Example: { "X1" : 0.0, "X3" : 0.0 }

                    Arguments:
                        value (dict[str, float]): dict of values to set (key: str, value: float)
                        nodes (list[int]): list of nodes.


        5. setValues(self: libaster.SimpleFieldOnNodesReal, value: dict[str, float], groupsOfNodes: list[str] = []) -> None


                    Set values of the field where components and values are given as a dict.
                    If the component is not present in the field then it is discarded
                    Example: { "X1" : 0.0, "X3" : 0.0 }

                    Arguments:
                        value (dict[str, float]): dict of values to set (key: str, value: float)
                        groupsOfNodes (list[str]): list of groups. If empty, the full mesh is considered


        6. setValues(self: libaster.SimpleFieldOnNodesReal, value: float) -> None


                    Set the value everywhere.

                    Arguments:
                        value [float]: value to set everywhere.
        """

    def toFieldOnNodes(self, *args, **kwargs):
        """Overloaded function.

        1. toFieldOnNodes(self: libaster.SimpleFieldOnNodesReal) -> libaster.FieldOnNodesReal


        Convert to FieldOnNodes

        Returns:
            FieldOnNodesReal: field converted


        2. toFieldOnNodes(self: libaster.SimpleFieldOnNodesReal, dofNum: libaster.BaseDOFNumbering, zeroExtension: bool = <PythonBool.FALSE: 0>) -> libaster.FieldOnNodesReal


        Convert to FieldOnNodes

        Arguments:
            dofNum (BaseDOFNumbering): DOF numbering used to build FieldOnNodes
            zeroExtension (bool): true if field can be extended to zero (when missing values)

        Returns:
            FieldOnNodesReal: field converted
        """

    def toNumpy(self):
        """Returns two numpy arrays with shape ( number_of_components, space_dimension )
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.

        Where the mask is `False` the corresponding value is set to zero.

        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """

    def updateValuePointers(self):
        pass


# class SimpleFieldOnNodesComplex in libaster


class SimpleFieldOnNodesComplex(DataField):
    pass

    # Method resolution order:
    #     SimpleFieldOnNodesComplex
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getitem__(self, arg0):
        pass

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.SimpleFieldOnNodesComplex, arg0: libaster.BaseMesh, arg1: str, arg2: list[str], arg3: bool) -> None

        2. __init__(self: libaster.SimpleFieldOnNodesComplex, arg0: tuple) -> None
        """

    def __setitem__(self, arg0, arg1):
        pass

    def __setstate__(self, arg0):
        pass

    def getComponent(self, arg0):
        pass

    def getComponents(self):
        pass

    def getLocalization(self):
        pass

    def getMesh(self):
        """Returns base mesh"""

    def getNumberOfComponents(self):
        pass

    def getNumberOfNodes(self):
        pass

    def getPhysicalQuantity(self):
        pass

    def hasComponent(self, arg0):
        pass

    def toNumpy(self):
        """Returns two numpy arrays with shape ( number_of_components, space_dimension )
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.

        Where the mask is `False` the corresponding value is set to zero.

        Returns:
            ndarray (complex): Field values.
            ndarray (bool): Mask for the field values.
        """

    def updateValuePointers(self):
        pass


# class Table in libaster


class Table(DataStructure):
    pass

    # Method resolution order:
    #     Table
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Table) -> None

        2. __init__(self: libaster.Table, arg0: str) -> None
        """

    def getColumnType(self, param):
        """Return the type of values in a column.

        Arguments:
            param (str): Parameter name.

        Returns:
            str: "I" for integers, "R" for reals, "C" for complex, "Knn" for strings.
        """

    def getNumberOfLines(self):
        """Returns the number of lines of the table.

        Returns:
            int: Number of lines.
        """

    def getParameters(self):
        """Return the parameters names.

        Returns:
            list[str]: Names of the parameters.
        """

    def getValues(self, arg0):
        """For internal use only. See *get_column()*."""


# class TableOfFunctions in libaster


class TableOfFunctions(Table):
    pass

    # Method resolution order:
    #     TableOfFunctions
    #     Table
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.TableOfFunctions) -> None

        2. __init__(self: libaster.TableOfFunctions, arg0: str) -> None
        """

    def addFunction(self, arg0):
        """Add a function into the table."""

    def getFunction(self, pos):
        """Returns the function stored at a given position.

        Arguments:
            pos [int]: Index of the function to return (0-based).

        Returns:
            *Function*: Function stored.
        """

    def getNumberOfFunctions(self):
        """Returns the number of functions stored in the table.

        Returns:
            int: Number of functions.
        """


# class TableContainer in libaster


class TableContainer(Table):
    pass

    # Method resolution order:
    #     TableContainer
    #     Table
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.TableContainer) -> None

        2. __init__(self: libaster.TableContainer, arg0: str) -> None
        """

    def addObject(self, *args, **kwargs):
        """Overloaded function.

        1. addObject(self: libaster.TableContainer, arg0: str, arg1: ElementaryMatrix<double, (PhysicalQuantityEnum)4>) -> None

        2. addObject(self: libaster.TableContainer, arg0: str, arg1: ElementaryVector<double, (PhysicalQuantityEnum)4>) -> None

        3. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.FieldOnCellsReal) -> None

        4. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.FieldOnNodesReal) -> None

        5. addObject(self: libaster.TableContainer, arg0: str, arg1: FunctionComplex) -> None

        6. addObject(self: libaster.TableContainer, arg0: str, arg1: GeneralizedAssemblyMatrix<double>) -> None

        7. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.DataField) -> None

        8. addObject(self: libaster.TableContainer, arg0: str, arg1: ModeResult) -> None

        9. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.ConstantFieldOnCellsReal) -> None

        10. addObject(self: libaster.TableContainer, arg0: str, arg1: Function2D) -> None

        11. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.Table) -> None

        12. addObject(self: libaster.TableContainer, name: str, object: Function) -> None


                    Store a *DataStructure* in the table.

                    Arguments:
                        name (str): String to identify the object in the table.
                        object (misc): Object to be inserted, can be a Function, Mesh, Fields...
        """

    def getConstantFieldOnCellsReal(self, arg0):
        pass

    def getDataField(self, arg0):
        pass

    def getElementaryMatrixDisplacementReal(self, arg0):
        pass

    def getElementaryVectorDisplacementReal(self, arg0):
        pass

    def getFieldOnCellsReal(self, arg0):
        pass

    def getFieldOnNodesReal(self, arg0):
        pass

    def getFunction(self, arg0):
        pass

    def getFunction2D(self, arg0):
        pass

    def getFunctionComplex(self, arg0):
        pass

    def getGeneralizedAssemblyMatrix(self, arg0):
        pass

    def getModeResult(self, arg0):
        pass

    def getTable(self, arg0):
        pass


# class TimesList in libaster


class TimesList(DataStructure):
    pass

    # Method resolution order:
    #     TimesList
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.TimesList) -> None

        2. __init__(self: libaster.TimesList, arg0: str) -> None
        """

    def getValues(self):
        """Return the current times values.

        Returns:
            list[float]: Current values.
        """

    # ----------------------------------------------------------------------
    # Data descriptors defined here:

    @property
    def stepper(self):
        pass


# class GeneralizedDOFNumbering in libaster


class GeneralizedDOFNumbering(DataStructure):
    pass

    # Method resolution order:
    #     GeneralizedDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedDOFNumbering) -> None

        2. __init__(self: libaster.GeneralizedDOFNumbering, arg0: str) -> None
        """

    def getGeneralizedModel(self):
        pass

    def getModalBasis(self):
        pass

    def getMorseStorage(self):
        pass

    def setGeneralizedModel(self, arg0):
        pass

    def setModalBasis(self, *args, **kwargs):
        """Overloaded function.

        1. setModalBasis(self: libaster.GeneralizedDOFNumbering, arg0: ModeResult) -> bool

        2. setModalBasis(self: libaster.GeneralizedDOFNumbering, arg0: GeneralizedModeResult) -> bool
        """


# class FluidStructureInteraction in libaster


class FluidStructureInteraction(DataStructure):
    pass

    # Method resolution order:
    #     FluidStructureInteraction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FluidStructureInteraction) -> None

        2. __init__(self: libaster.FluidStructureInteraction, arg0: str) -> None
        """


# class TurbulentSpectrum in libaster


class TurbulentSpectrum(DataStructure):
    pass

    # Method resolution order:
    #     TurbulentSpectrum
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.TurbulentSpectrum) -> None

        2. __init__(self: libaster.TurbulentSpectrum, arg0: str) -> None
        """


# class GenericFunction in libaster


class GenericFunction(DataStructure):
    pass

    # Method resolution order:
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getProperties(self):
        """Returns the properties of the function.

        Returns:
            tuple[str]: Tuple containing: type of the function (same as `getType()`),
            type of interpolation, parameter name, result name,
            type of extrapolation, object name (same as `getName()`).
        """

    def setExtrapolation(self, type):
        """Define the type of extrapolation.

        Supported extrapolation types are: "L" for linear, "C" for constant and
        "E" for no extrapolation allowed.

        Arguments:
            type (str): Type of extrapolation on left and on right. Examples: "CC",
                "LE"...
        """


# class ListOfLoads in libaster


class ListOfLoads(DataStructure):
    pass

    # Method resolution order:
    #     ListOfLoads
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ListOfLoads) -> None

        2. __init__(self: libaster.ListOfLoads, arg0: str) -> None

        3. __init__(self: libaster.ListOfLoads, arg0: Model) -> None

        4. __init__(self: libaster.ListOfLoads, arg0: str, arg1: Model) -> None
        """

    def addDirichletBC(self, *args, **kwargs):
        """Overloaded function.

        1. addDirichletBC(self: libaster.ListOfLoads, arg0: DirichletBC) -> None

        2. addDirichletBC(self: libaster.ListOfLoads, arg0: DirichletBC, arg1: Function, arg2: str) -> None

        3. addDirichletBC(self: libaster.ListOfLoads, arg0: DirichletBC, arg1: Formula, arg2: str) -> None

        4. addDirichletBC(self: libaster.ListOfLoads, arg0: DirichletBC, arg1: Function2D, arg2: str) -> None

        5. addDirichletBC(self: libaster.ListOfLoads, arg0: DirichletBC, arg1: str) -> None
        """

    def addLoad(self, *args, **kwargs):
        """Overloaded function.

        1. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<double> >) -> None

        2. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None

        3. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<double> >, arg1: str) -> None

        4. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function, arg2: str) -> None

        5. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<double> >, arg1: Formula, arg2: str) -> None

        6. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function2D, arg2: str) -> None

        7. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: str) -> None

        8. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function, arg2: str) -> None

        9. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula, arg2: str) -> None

        10. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D, arg2: str) -> None

        11. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >) -> None

        12. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function) -> None

        13. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Formula) -> None

        14. addLoad(self: libaster.ListOfLoads, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function2D) -> None

        15. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<double> >) -> None

        16. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<double> >, arg1: Function) -> None

        17. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<double> >, arg1: Formula) -> None

        18. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<double> >, arg1: Function2D) -> None

        19. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None

        20. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function) -> None

        21. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula) -> None

        22. addLoad(self: libaster.ListOfLoads, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D) -> None

        23. addLoad(self: libaster.ListOfLoads, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >) -> None

        24. addLoad(self: libaster.ListOfLoads, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function) -> None

        25. addLoad(self: libaster.ListOfLoads, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Formula) -> None

        26. addLoad(self: libaster.ListOfLoads, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function2D) -> None

        27. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >) -> None

        28. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: str) -> None

        29. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function, arg2: str) -> None

        30. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: Formula, arg2: str) -> None

        31. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function2D, arg2: str) -> None

        32. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: str) -> None

        33. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function, arg2: str) -> None

        34. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula, arg2: str) -> None

        35. addLoad(self: libaster.ListOfLoads, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D, arg2: str) -> None

        36. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >) -> None

        37. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: Function) -> None

        38. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: Formula) -> None

        39. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: Function2D) -> None

        40. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None

        41. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function) -> None

        42. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula) -> None

        43. addLoad(self: libaster.ListOfLoads, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D) -> None
        """

    def getDirichletBCs(self):
        """Return list of DirichletBCs

        Returns:
            ListDiriBC: a list of DirichletBC
        """

    def getLoadNames(self):
        """Returns list of load's names.

        Returns:
            list[str]: list of load's names.
        """

    def getMechanicalLoadsFunction(self):
        """Return list of Function mechanical loads

        Returns:
            ListMecaLoadFunction: a list of Function mechanical loads
        """

    def getMechanicalLoadsReal(self):
        """Return list of real mechanical loads

        Returns:
            ListMecaLoadReal: a list of real mechanical loads
        """

    def getModel(self):
        """Return the model used

        Returns:
            Model: model used
        """

    def getParallelMechanicalLoadsFunction(self):
        """Return list of function parallel mechanical loads

        Returns:
            ListParaMecaLoadFunction: a list of function parallel mechanical loads
        """

    def getParallelMechanicalLoadsReal(self):
        """Return list of real parallel mechanical loads

        Returns:
            ListParaMecaLoadReal: a list of real parallel mechanical loads
        """

    def hasDifferential(self):
        """Return True if there are DIDI loads or DIDI Dirichlet BCs"""

    def hasDirichletBC(self):
        """Dirichlet BCs have been added or not ?

        Returns:
            bool: True if Dirichlet BCs have been added
        """

    def hasExternalLoad(self):
        """External load (= not Dirichlet BCs) have been added or not ?

        Returns:
            bool: True if External load have been added
        """

    def isBuilt(self):
        """The list of loads has been built or not.

        Returns:
            bool: True if has been built already.
        """

    def setDifferentialDisplacement(self, arg0):
        """Set differential displacement field for DIDI loads"""


# class BaseFunction in libaster


class BaseFunction(GenericFunction):
    pass

    # Method resolution order:
    #     BaseFunction
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getValues(self):
        """Return a list of the values of the function as (x1, x2, ..., y1, y2, ...)

        Returns:
            list[float]: List of values (size = 2 * *size()*).
        """

    def setAsConstant(self):
        """To be called for a constant function."""

    def setInterpolation(self, type):
        """Define the type of interpolation.

        Supported interpolation types are: "LIN" for linear, "LOG" for logarithmic and
        "NON" for no interpolation allowed.

        Arguments:
            type (str): Type of interpolation for abscissa and ordinates. Examples: "LIN LIN",
                "LIN LOG"...
        """

    def setParameterName(self, name):
        """Define the name of the abscissa.

        Arguments:
            name (str): Name of the abscissa.
        """

    def setResultName(self, name):
        """Define the name of the ordinates.

        Arguments:
            name (str): Name of the ordinates.
        """

    def setValues(self, absc, ordo):
        """Set the values of abscissa and ordinates.

        If the function already exists, its size can not be changed.

        Arguments:
            absc (list): List of abscissa.
            ordo (list): List of ordinates.
        """

    def size(self):
        """Return the number of points of the function.

        Returns:
            int: Number of points.
        """


# class Function in libaster


class Function(BaseFunction):
    pass

    # Method resolution order:
    #     Function
    #     BaseFunction
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Function) -> None

        2. __init__(self: libaster.Function, arg0: str) -> None
        """


# class FunctionComplex in libaster


class FunctionComplex(BaseFunction):
    pass

    # Method resolution order:
    #     FunctionComplex
    #     BaseFunction
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FunctionComplex) -> None

        2. __init__(self: libaster.FunctionComplex, arg0: str) -> None
        """

    def setValues(self, *args, **kwargs):
        """Overloaded function.

        1. setValues(self: libaster.FunctionComplex, arg0: list[float], arg1: list[float]) -> None

        2. setValues(self: libaster.FunctionComplex, arg0: list[float], arg1: list[complex]) -> None
        """


# class Formula in libaster


class Formula(GenericFunction):
    pass

    # Method resolution order:
    #     Formula
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Formula) -> None

        2. __init__(self: libaster.Formula, arg0: str) -> None
        """

    def evaluate(self, *val):
        """Evaluate the formula with the given variables values.

        Arguments:
            val (list[float]): List of the values of the variables.

        Returns:
            float/complex: Value of the formula for these values.
        """

    def getContext(self):
        """Return the context used to evaluate the formula.

        Returns:
            dict: Context used for evaluation.
        """

    def getExpression(self):
        """Return expression of the formula.

        Returns:
            str: Expression of the formula.
        """

    def getProperties(self):
        """Return the properties of the formula (for compatibility with function objects).

        Returns:
            list[str]: List of 6 strings as function objects contain.
        """

    def getVariables(self):
        """Return the variables names.

        Returns:
            list[str]: List of the names of the variables.
        """

    def setComplex(self):
        """Set the type of the formula as complex."""

    def setContext(self, context):
        """Define the context holding objects required to evaluate the expression.

        Arguments:
            context (dict): Context for the evaluation.
        """

    def setExpression(self, expression):
        """Define the expression of the formula.

        If the expression uses non builtin objects, the evaluation context must be
        defined using `:func:setContext`.

        Arguments:
            expression (str): Expression of the formula.
        """

    def setVariables(self, varnames):
        """Define the variables names.

        Arguments:
            varnames (list[str]): List of variables names.
        """


# built-in function jeveux_init in libaster


def jeveux_init(mpi_comm):
    """Initialize the memory manager (Jeveux).

    Arguments:
        mpi_comm (int): Identifier of MPI communicator (from ``py2f()``).
    """


# built-in function jeveux_finalize in libaster


def jeveux_finalize(arg0):
    """Finalize the memory manager (Jeveux)."""


# built-in function call_oper in libaster


def call_oper(syntax, jxveri):
    """Call a Fortran operator ('op' subroutine).

    Arguments:
        syntax (CommandSyntax): Object containing the user syntax.
        jxveri (int): If non null `JXVERI` is called after calling the operator.
    """


# built-in function call_oper_init in libaster


def call_oper_init():
    """Execute initializations before and after operator but without executing any
    operator.
    """


# built-in function cmd_ctxt_enter in libaster


def cmd_ctxt_enter():
    """Call Fortran 'cmd_ctxt_enter' subroutine."""


# built-in function cmd_ctxt_exit in libaster


def cmd_ctxt_exit():
    """Call Fortran 'cmd_ctxt_exit' subroutine."""


# built-in function write in libaster


def write(text):
    """Print a string using the fortran subroutine.

    Arguments:
        text (str): Text to be printed.
    """


# built-in function affich in libaster


def affich(code, text):
    """Print a string using the fortran subroutine on an internal file.

    Arguments:
        code (str): Code name of the internal file : 'MESSAGE' or 'CODE'.
        text (str): Text to be printed.
    """


# built-in function jeveux_status in libaster


def jeveux_status():
    """Return the status of the Jeveux memory manager.

    Returns:
        int: 0 after initialization and after shutdown, 1 during the execution.
    """


# built-in function jeveux_delete in libaster


def jeveux_delete(prefix):
    """Force manual deletion of Jeveux objects.

    Warning: Use only for debugging usages, it is dangerous for objects integrity
    and cpu consuming.

    Arguments:
        prefix (str): Root name of the Jeveux datastructure.
    """


# built-in function deleteTemporaryObjects in libaster


def deleteTemporaryObjects():
    """Delete temporary Jeveux objects"""


# built-in function deleteCachedObjects in libaster


def deleteCachedObjects():
    """Delete temporary and cached Jeveux objects (temporary, matrix, base, ...)"""


# built-in function onFatalError in libaster


def onFatalError(value=""):
    """Get/set the behavior in case of error.

    Arguments:
        value (str, optional): Set the new behavior in case of error (one of "ABORT",
            "EXCEPTION", "EXCEPTION+VALID" or "INIT"). If `value` is not provided,
            the current behavior is returned.

    Returns:
        str: Current value
    """


# built-in function matfpe in libaster


def matfpe(value):
    """Enable or disable floating point exceptions.

    Arguments:
        value (int): -1 to disable the FPE interception, 1 to enable FPE detection.
    """


# built-in function fe_invalid in libaster


def fe_invalid(value):
    """Enable or disable FE_INVALID exception.

    Arguments:
        value (int): -1 to disable the interception, 1 to enable detection.
    """


# built-in function reset_tpmax in libaster


def reset_tpmax():
    """Reset the cache value of tpmax (for internal use only)."""


# built-in function asmpi_set in libaster


def asmpi_set(arg0):
    """Set the current MPI communicator.

    Arguments:
        comm (int): id of the communicator.
    """


# built-in function asmpi_get in libaster


def asmpi_get():
    """Get the current MPI communicator.

    Returns:
        comm (int): id of the communicator.
    """


# built-in function asmpi_free in libaster


def asmpi_free(arg0):
    """Free the MPI communicator in argument.

    Arguments:
        comm (int): id of the communicator.
    """


# built-in function asmpi_split in libaster


def asmpi_split(arg0, arg1, arg2):
    """Split the MPI communicator in argument.

    Arguments:
        comm (int): id of the parent communicator to split.
        color (int): color to which the calling process will belong.
        name (str): name of the new communicator.

    Returns:
        comm (int) : id of the communicator.
    """


# built-in function asmpi_info in libaster


def asmpi_info(arg0):
    """Return the rank and size of the MPI communicator.

    Arguments:
        comm (int): id of the communicator.

    Returns:
        rank (int) : rank of the communicator.
        size (int) : size of the communicator.
    """


# class Function2D in libaster


class Function2D(GenericFunction):
    pass

    # Method resolution order:
    #     Function2D
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Function2D) -> None

        2. __init__(self: libaster.Function2D, arg0: str) -> None
        """

    def getParameters(self):
        """Return a list of the values of the parameter as (x1, x2, ...)

        Returns:
            list[float]: List of values (size = number of functions).
        """

    def getProperties(self):
        """Returns the properties of the function.

        Returns:
            tuple[str]: Tuple containing: type of the function (same as `getType()`),
            type of interpolation, parameter name, result name,
            type of extrapolation, object name (same as `getName()`),
            parameter name of functions + a list of dict for each functions that contain
            the type of interpolation and extrapolation.
        """

    def getValues(self):
        """Return a list of the values of the functions as [F1, F2, ...]
        where Fi is (x1, x2, ..., y1, y2, ...).

        Returns:
            list[list [float]]: List of values (size = number of functions).
        """


# class Contact in libaster


class Contact(DataStructure):
    pass

    # Method resolution order:
    #     Contact
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Contact, arg0: str, arg1: Model) -> None

        2. __init__(self: libaster.Contact, arg0: Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getModel(self):
        pass


# class ContactAlgo in libaster


class ContactAlgo:
    """Enumeration for contact algorithm."""

    # Method resolution order:
    #     ContactAlgo
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Lagrangian = 0

    Nitsche = 1

    Penalization = 2


# class ContactVariant in libaster


class ContactVariant:
    """Enumeration for contact variant."""

    # Method resolution order:
    #     ContactVariant
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Classic = 4

    Empty = 0

    Fast = 1

    Robust = 2

    Symetric = 3


# class ContactType in libaster


class ContactType:
    """Enumeration for contact type."""

    # Method resolution order:
    #     ContactType
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Bilateral = 1

    Unilateral = 0


# class FrictionAlgo in libaster


class FrictionAlgo:
    """Enumeration for friction algorithm."""

    # Method resolution order:
    #     FrictionAlgo
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Lagrangian = 0

    Nitsche = 1

    Penalization = 2


# class FrictionType in libaster


class FrictionType:
    """Enumeration for friction type."""

    # Method resolution order:
    #     FrictionType
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Coulomb = 2

    Stick = 3

    Tresca = 1

    Without = 0


# class PairingAlgo in libaster


class PairingAlgo:
    """Enumeration for pairing algorithm."""

    # Method resolution order:
    #     PairingAlgo
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Mortar = 0


# class InitialState in libaster


class InitialState:
    """Enumeration for initial state."""

    # Method resolution order:
    #     InitialState
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Interpenetrated = 0

    No = 1

    Yes = 2


# class JacobianType in libaster


class JacobianType:
    """Enumeration for jacobian type."""

    # Method resolution order:
    #     JacobianType
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Analytical = 0

    Perturbation = 1


# class IntegrationType in libaster


class IntegrationType:
    """Enumeration for integration type."""

    # Method resolution order:
    #     IntegrationType
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Elembased = 0

    Segbased = 1


# class ContactParameter in libaster


class ContactParameter:
    pass

    # Method resolution order:
    #     ContactParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def getAlgorithm(self):
        """Return the contact algorithm used. It is a value of an enum

        Returns:
            ContactAlgo: contact algorithm.
        """

    def getCoefficient(self):
        """Return the contact coefficient used. It is a value of a float

        Returns:
            float: contact coefficient.
        """

    def getIntegrationType(self):
        """Return how the integration is made. It is a value of an enum

        Returns:
            IntegrationType: Integration type.
        """

    def getJacobianType(self):
        """Return how the Jacobian is computed. It is a value of an enum

        Returns:
            JacobianType: Jacobian type.
        """

    def getType(self):
        """Return the contact type used. It is a value of an enum

        Returns:
            ContactType: contact type.
        """

    def getVariant(self):
        """Return the contact variant used. It is a value of an enum

        Returns:
            ContactVariant: contact variant.
        """

    def setAlgorithm(self, algo):
        """Set the contact algorithm used. It is a value of an enum

        Arguments:
            ContactAlgo: contact algorithm.
        """

    def setCoefficient(self, coeff):
        """Set the contact coefficient used. It is a value of a float

        Arguments:
            float: contact coefficient.
        """

    def setIntegrationType(self, type):
        """Set how the integration is made. It is a value of an enum

        Arguments:
            IntegrationType: Integration type.
        """

    def setJacobianType(self, type):
        """Set how the Jacobian is computed. It is a value of an enum

        Arguments:
            JacobianType: Jacobian type.
        """

    def setType(self, type):
        """Set the contact type used. It is a value of an enum

        Arguments:
            ContactType: contact type.
        """

    def setVariant(self, variant):
        """Set the contact variant used. It is a value of an enum

        Arguments:
            ContactVariant: contact variant.
        """


# class FrictionParameter in libaster


class FrictionParameter:
    pass

    # Method resolution order:
    #     FrictionParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def getAlgorithm(self):
        """Return the Friction algorithm used. It is a value of an enum

        Returns:
            FrictionAlgo: Friction algorithm.
        """

    def getCoefficient(self):
        """Return the Friction coefficient used. It is a value of a float

        Returns:
            float: Friction coefficient.
        """

    def getCoulomb(self):
        """Return the Coulomb coefficient used. It is a value of a float

        Returns:
            float: Coulomb coefficient.
        """

    def getTresca(self):
        """Return the Tresca coefficient used. It is a value of a float

        Returns:
            float: Tresca coefficient.
        """

    def getType(self):
        """Return the Friction type used. It is a value of an enum

        Returns:
            FrictionType: Friction type.
        """

    def setAlgorithm(self, algo):
        """Set the Friction algorithm used. It is a value of an enum

        Arguments:
            FrictionAlgo: Friction algorithm.
        """

    def setCoefficient(self, coeff):
        """Set the Friction coefficient used. It is a value of a float

        Arguments:
            float: Friction coefficient.
        """

    def setCoulomb(self, coulomb):
        """Set the Coulomb coefficient used. It is a value of a float

        Arguments:
            float: Coulomb coefficient.
        """

    def setTresca(self, tresca):
        """Set the Tresca coefficient used. It is a value of a float

        Arguments:
            float: Tresca coefficient.
        """

    def setType(self, type):
        """Set the Friction type used. It is a value of an enum

        Arguments:
            FrictionType: Friction type.
        """

    # ----------------------------------------------------------------------
    # Data descriptors defined here:

    @property
    def hasFriction(self):
        """bool: enable or disable the use of friction."""


# class PairingParameter in libaster


class PairingParameter:
    pass

    # Method resolution order:
    #     PairingParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def getAlgorithm(self):
        """Return the Pairing algorithm used. It is a value of an enum

        Returns:
            PairingAlgo: Pairing algorithm.
        """

    def getAreaIntersectionTolerance(self):
        """Return the tolerance used for the intersection area criterium. It is a value of a float

        Returns:
            float: area intersection tolerance.
        """

    def getDistanceFunction(self):
        """Return the fictive distance function. It is a value of a pointer

        Returns:
            GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        """

    def getDistanceRatio(self):
        """Return the pairing distance ratio used. It is a value of a float

        Returns:
            float: pairing distance.
        """

    def getElementaryCharacteristics(self):
        """Return the elementary characteristics. It is a value of a pointer

        Returns:
            ElementaryCharacteristicsPtr: cara_elel pointer.
        """

    def getInitialState(self):
        """Return the initial contact state. It is a value of an enum

        Returns:
            InitialState: Initial contact state.
        """

    def getPairingMethod(self):
        """Return the pairing method used. It is a value of an enum

        Returns:
            PairingMethod: pairing method.
        """

    def getPairingTolerance(self):
        """Return the pairing tolerance used. It is a value of a float

        Returns:
            float: pairing tolerance.
        """

    def setAlgorithm(self, algo):
        """Set the Pairing algorithm used. It is a value of an enum

        Arguments:
            PairingAlgo: Pairing algorithm.
        """

    def setAreaIntersectionTolerance(self, area_tole):
        """Return the tolerance used for the intersection area criterium. It is a value of a float

        Arguments:
            float: area intersection tolerance.
        """

    def setDistanceFunction(self, dist_supp):
        """Set the fictive distance function. It is a value of a pointer

        Arguments:
            GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        """

    def setDistanceRatio(self, dist_ratio):
        """Set the pairing distance ratio used. It is a value of a float

        Arguments:
            float: pairing distance ratio.
        """

    def setElementaryCharacteristics(self, cara):
        """Set the elementary characteristics. It is a value of a pointer

        Arguments:
            ElementaryCharacteristicsPtr: cara_elel pointer.
        """

    def setInitialState(self, cont_init):
        """Set the initial contact state. It is a value of an enum

        Arguments:
            InitialState: Initial contact state.
        """

    def setPairingMethod(self, pair_method):
        """Set the cpairing method used. It is a value of an enum

        Arguments:
            PairingMethod: pairing method.
        """

    def setPairingTolerance(self, pair_tole):
        """Set the pairing tolerance used. It is a value of a float

        Arguments:
            float: pairing tolerance.
        """

    # ----------------------------------------------------------------------
    # Data descriptors defined here:

    @property
    def hasBeamDistance(self):
        """bool: enable or disable the use of a fictive distance for beam."""

    @property
    def hasShellDistance(self):
        """bool: enable or disable the use of a fictive distance for shell."""


# class ContactNew in libaster


class ContactNew(DSWithCppPickling):
    pass

    # Method resolution order:
    #     ContactNew
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ContactNew, arg0: str, arg1: Model) -> None

        2. __init__(self: libaster.ContactNew, arg0: Model) -> None

        3. __init__(self: libaster.ContactNew, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def appendContactZone(self, zone):
        """Append a new contact zone to the contact definition

        Arguments:
            zone (ContactZone): contact zone to append
        """

    def build(self):
        """Build and check internal objects"""

    def getContactZone(self, zone_id):
        """Return the specified contact zone

        Arguments:
            zone_id (int): index of the contact zone (0-based)

        Returns:
            ContactZone: contact zone.
        """

    def getContactZones(self):
        """Return the list of contact zones

        Returns:
            List[ContactZone]: List of contact zones.
        """

    def getFiniteElementDescriptor(self):
        """Return the finite element descriptor to define virtual cells for Lagrange multipliers

        Returns:
            FiniteElementDescriptor: fed.
        """

    def getMesh(self):
        """Return the mesh used in the contact definition

        Returns:
            Mesh: mesh.
        """

    def getModel(self):
        """Return the model used in the contact definition

        Returns:
            Model: model.
        """

    def getNumberOfContactZones(self):
        """Return the number of contact zones used

        Returns:
            inter: number of contact zones.
        """

    def getVerbosity(self):
        """Get level of verbosity:*
              0- without
              1- normal
              2- detailled

        Returns:
            integer: level of verbosity
        """

    def isParallel(self):
        """bool: true if parallel contact."""

    def setVerbosity(self, level):
        """Set level of verbosity:
              0- without
              1- normal (default)
              2- detailled

        Arguments:
            level (int): level of verbosity
        """

    # ----------------------------------------------------------------------
    # Data descriptors defined here:

    @property
    def hasFriction(self):
        """bool: enable or disable the use of friction."""

    @property
    def hasSmoothing(self):
        """bool: enable or disable  the use of smoothing."""


# class FrictionNew in libaster


class FrictionNew(ContactNew):
    pass

    # Method resolution order:
    #     FrictionNew
    #     ContactNew
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FrictionNew, arg0: str, arg1: Model) -> None

        2. __init__(self: libaster.FrictionNew, arg0: Model) -> None

        3. __init__(self: libaster.FrictionNew, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass


# class ContactZone in libaster


class ContactZone(DSWithCppPickling):
    """Object to define a zone of contact."""

    # Method resolution order:
    #     ContactZone
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ContactZone, arg0: str) -> None

        2. __init__(self: libaster.ContactZone) -> None

        3. __init__(self: libaster.ContactZone, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def build(self, arg0):
        """Build and check internal objects

        Returns:
            bool: success or failure
        """

    def getContactParameter(self):
        """Get contact parameters defining method, coefficient...

        Returns:
            ContactParameter: contact parameters
        """

    def getFrictionParameter(self):
        """Get friction parameters defining method, coefficient...

        Returns:
            FrictionParameter: friction parameters
        """

    def getMesh(self):
        """Return the mesh used in the contact zone definition

        Returns:
            BaseMesh: mesh
        """

    def getMeshPairing(self):
        """Get pairing of surface meshes

        Returns:
            MeshPairing: mesh pairing
        """

    def getModel(self):
        """Return the model used in the contact zone definition

        Returns:
            Model: model
        """

    def getPairingParameter(self):
        """Get pairing parameters defining algorithm, distance...

        Returns:
            PairingParameter: pairing parameters
        """

    def getSlaveCells(self):
        """Get slave's cells index

        Returns:
            list[int]: slave's cells index
        """

    def getSlaveNodes(self):
        """Get slave's nodes index

        Returns:
            list[int]: slave's nodes index
        """

    def getVerbosity(self):
        """Get level of verbosity
        0- without
        1- normal
        2- detailled

        Returns:
            int: level of verbosity
        """

    def setContactParameter(self, contParam):
        """Set contact parameters defining method, coefficient...

        Arguments:
            contParam (ContactParameter) : contact parameters
        """

    def setExcludedSlaveGroupOfCells(self, cellGroupsName):
        """Set excluded groups of cells on slave side

        Arguments:
            cellGroupsName (str) : excluded groups' names
        """

    def setExcludedSlaveGroupOfNodes(self, nodeGroupsName):
        """Set excluded groups of nodes on slave side

        Arguments:
            nodeGroupsName (str) : excluded groups' names
        """

    def setFrictionParameter(self, fricParam):
        """Set friction parameters defining method, coefficient...

        Arguments:
            fricParam (FrictionParameter) : friction parameters
        """

    def setMasterGroupOfCells(self, master_name):
        """Set master's name of group of cells

        Arguments:
            master_name (str) : name of group for master cells
        """

    def setPairingParameter(self, pairParam):
        """Set pairing parameters defining algorithm, distance...

        Arguments:
            pairParam (PairingParameter) : pairing parameters
        """

    def setSlaveGroupOfCells(self, slave_name):
        """Set slave's name of group of cells

        Arguments:
            slave_name (str) : name of group for slave cells
        """

    def setVerbosity(self, level):
        """Set level of verbosity
        0- without
        1- normal (default)
        2- detailled

        Arguments:
            level (int) : level of verbosity
        """

    # ----------------------------------------------------------------------
    # Data descriptors defined here:

    @property
    def checkNormals(self):
        """bool: attribute that holds the checking of outwards normals."""

    @property
    def hasFriction(self):
        """Get status of friction

        Returns:
            bool: friction or not
        """

    @property
    def hasSmoothing(self):
        """Smoothing of normals

        Returns:
            bool: smoothing or not
        """


# class MeshPairing in libaster


class MeshPairing(DSWithCppPickling):
    """Object to create a pairing operator between two meshed surfaces."""

    # Method resolution order:
    #     MeshPairing
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MeshPairing, arg0: str) -> None

        2. __init__(self: libaster.MeshPairing) -> None

        3. __init__(self: libaster.MeshPairing, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def checkNormals(self, model):
        """Check orientation of normals

        Arguments:
            ModelPtr: a pointer to the model

        Returns:
            nothing
        """

    def compute(self, dist_pairing=-1.0, pair_tole=1e-08, area_tole=1e-08):
        """Compute pairing

        Arguments:
            dist_pairing (real): tolerance from DIST_RATIO (projection outside cell)
            pair_tole (real): tolerance for pairing
            area_tole (real): tolerance for removing intersection cells
        """

    def getIntersectionArea(self, *args, **kwargs):
        """Overloaded function.

        1. getIntersectionArea(self: libaster.MeshPairing, indexPair: int) -> float


        Compute intersection of area

        Arguments:
            indexPair (integer): index of pair

        Returns:
            double: area of intersection


        2. getIntersectionArea(self: libaster.MeshPairing, indexPair: int) -> float


        Get area of intersection for a given pair

        Arguments:
            indexPair (integer): index of pair

        Returns:
            real: area of intersection
        """

    def getIntersectionPoints(self, indexPair, CoordinatesSpace=1):
        """Get coordinates of intersection points for a given pair

        Arguments:
            indexPair (integer): index of pair
            CoordinatesSpace (CoordinatesSpace): space to describe coordinates

        Returns:
            list[list]: coordinates in given space
        """

    def getListOfPairs(self):
        """Get pairs

        Returns:
            list: pairs (slave-master)
        """

    def getMasterCellNeighbors(self, cell_number):
        """Get the master cells in the neighbor of a given master cell number

        Arguments:
            int: master cell number

        Returns:
            list: master neighbors cells
        """

    def getMasterCellsFromNode(self, node_number):
        """Get the master cells associated with a node number

        Arguments:
            int: node number

        Returns:
            list: master cells associated
        """

    def getMesh(self):
        """Return the mesh

        Returns:
            Mesh: mesh.
        """

    def getNumberOfIntersectionPoints(self, *args, **kwargs):
        """Overloaded function.

        1. getNumberOfIntersectionPoints(self: libaster.MeshPairing, indexPair: int) -> int


        Get number of intersection points

        Arguments:
            indexPair (integer): index of pair

        Returns:
            integer: number of intersection points


        2. getNumberOfIntersectionPoints(self: libaster.MeshPairing) -> list[int]


                Get number of intersection points of all pairs

                Returns:
                    list: number of intersection points
        """

    def getNumberOfPairs(self):
        """Get number of pairs

        Returns:
            integer: number of pairs
        """

    def getQuadraturePoints(self, indexPair):
        """Get coordinates of quadrature points for a given pair in global space

        Arguments:
            indexPair (integer): index of pair

        Returns:
            list: quadrature points
        """

    def getSlaveCellNeighbors(self, cell_number):
        """Get the slave cells in the neighbor of a given slave cell number

        Arguments:
            int: slave cell number

        Returns:
            list: slave neighbors cells
        """

    def getSlaveCellsFromNode(self, node_number):
        """Get the slave cells associated with a node number

        Arguments:
            int: node number

        Returns:
            list: slave cells associated
        """

    def getVerbosity(self):
        """Get level of verbosity

        Returns:
            integer: level of verbosity
        """

    def setExcludedSlaveGroupOfCells(self, groups):
        """Set excluded groups of cells on slave side

        Arguments:
            str: excluded groups' names
        """

    def setExcludedSlaveGroupOfNodes(self, groups):
        """Set excluded groups of nodes on slave side

        Arguments:
            str: excluded groups' names
        """

    def setMesh(self, mesh):
        """Set Mesh

        Arguments:
            mesh (BaseMesh): support mesh
        """

    def setMethod(self, method):
        """Set method of pairing

        Arguments:
            method (PairingMethod): method ("OLD", "Fast", "Robust)
        """

    def setPair(self, groupNameSlav, groupNameMast):
        """Set pair of meshed surfaces

        Arguments:
            groupNameSlav (str): slave's name
            groupNameMast (str): master's name
        """

    def setVerbosity(self, level):
        """Set level of verbosity
              0 - without
              1 - normal (default)
              2 - detailled (text)

        Arguments:
            level (integer): level of verbosity
        """


# class CoordinatesSpace in libaster


class CoordinatesSpace:
    """Type of coordinates: Slave or Global."""

    # Method resolution order:
    #     CoordinatesSpace
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Global = 1

    Slave = 0


# class PairingMethod in libaster


class PairingMethod:
    """Type of pairing: Fast, BrutForce and Legacy."""

    # Method resolution order:
    #     PairingMethod
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    BrutForce = 2

    Fast = 0

    Legacy = 1


# class ContactPairing in libaster


class ContactPairing(DataStructure):
    """Object to create contact pairing."""

    # Method resolution order:
    #     ContactPairing
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ContactPairing, arg0: str, arg1: libaster.ContactNew) -> None

        2. __init__(self: libaster.ContactPairing, arg0: libaster.ContactNew) -> None
        """

    def clearPairing(self, *args, **kwargs):
        """Overloaded function.

        1. clearPairing(self: libaster.ContactPairing) -> None


        Clean pairing for all zones

        Returns:
            bool: true if the pairing quantities are cleared


        2. clearPairing(self: libaster.ContactPairing, zone_index: int) -> None


        Clean pairing for a zone

        Arguments:
            zone_index(int) : index of zone

        Returns:
            bool: true if the pairing quantities are cleared
        """

    def compute(self, *args, **kwargs):
        """Overloaded function.

        1. compute(self: libaster.ContactPairing) -> bool


        Compute the pairing quantities on all zones

        Returns:
            bool: True if the pairing quantities are updated appropriately


        2. compute(self: libaster.ContactPairing, zone_index: int) -> bool


        Compute the pairing quantities on a zone

        Arguments:
            zone_index(int)

        Returns:
            bool: True if the pairing quantities are updated appropriately
        """

    def getCoordinates(self):
        """Coordinates of nodes used for pairing (almost always different from the intial mesh).

        Returns:
            MeshCoordinatesField: the coordinates field
        """

    def getFiniteElementDescriptor(self):
        """Return Finite Element Descriptor for virtual cells from pairing.

        Returns:
            FiniteElementDescriptor: finite element for virtual cells
        """

    def getIntersectionPoints(self, zone_index, CoordinatesSpace=1):
        """Get the intersection points between master and slave cells

        Arguments:
            zone_index(int) : index of zone
            CoordinatesSpace (CoordinatesSpace): space to describe coordinates

        Returns:
            list[pair]: list of pair of coordinates of intersection points
        """

    def getListOfPairs(self, *args, **kwargs):
        """Overloaded function.

        1. getListOfPairs(self: libaster.ContactPairing, zone_index: int) -> list[tuple[int, int]]


        Get list of contact pairs for a contact zone

        Arguments:
            zone_index(int)

        Returns:
            list[tuple[int, int]]: list of contact pairs


        2. getListOfPairs(self: libaster.ContactPairing) -> list[tuple[int, int]]


        Get list of contact pairs on all zones

        Returns:
            list[tuple[int, int]]: list of contact pairs
        """

    def getMesh(self):
        """Mesh

        Returns:
            BaseMesh: the mesh
        """

    def getNumberOfIntersectionPoints(self, zone_index):
        """Get list of the number of intersection points beetween a master and slave cells.

        Arguments:
            zone_index(int) : index of zone

        Returns:
            list: list of number of intersection points
        """

    def getNumberOfPairs(self, *args, **kwargs):
        """Overloaded function.

        1. getNumberOfPairs(self: libaster.ContactPairing) -> int


        Return number of pairs on all zones

        Returns:
            int: number of pairs


        2. getNumberOfPairs(self: libaster.ContactPairing, zone_index: int) -> int


        Return the number of pairs on a zone

        Arguments:
            zone_index(int)
        Returns:
            int: number of pairs
        """

    def getNumberOfZones(self):
        """Return the number of zones

        Returns:
            int: number of zones
        """

    def getVerbosity(self):
        """Get level of verbosity

        Returns:
            integer: level of verbosity
        """

    def setCoordinates(self, coordinates):
        """Set the mesh coordinates field

        Arguments:
            coordinates (MeshCoordinatesField) : coordinates to use for pairing
        """

    def setVerbosity(self, verbosity):
        """Set level of verbosity
              0 - without
              1 - normal (default)
              2 - detailled (text)

        Arguments:
            level (integer): level of verbosity
        """

    def updateCoordinates(self, disp):
        """Update the mesh coordinates given a displacement field

        Arguments:
            disp (FieldOnNodes) : field for displacement
        """


# class ContactComputation in libaster


class ContactComputation:
    pass

    # Method resolution order:
    #     ContactComputation
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, arg0):
        pass

    def __setstate__(self, arg0):
        pass

    def contactCoefficient(self):
        """Compute contact coefficients at the nodes of the slave surface based on values of COEF_CONT
        and COEF_FROT

        Returns:
            list[FieldOnNodesReal]: coefficients (COEF_CONT and COEF_FROT)
        """

    def contactData(self, pairing, material, initial_contact):
        """Compute contact data as input to compute contact forces and matrices.

        Arguments:
            pairing (ContactPairing): pairing object
            material (MaterialField): material field
            initial_contact (bool): True to use value in contact definition (CONTACT_INIT).

        Returns:
            FieldOnCellsReal: contact data
        """

    def getVerbosity(self):
        """Get level of verbosity
        0- without
        1- normal
        2- detailled

        Returns:
            int: level of verbosity
        """

    def setVerbosity(self, level):
        """Set level of verbosity
        0- without
        1- normal (default)
        2- detailled

        Arguments:
            level (int) : level of verbosity
        """


# class CouplingMethod in libaster


class CouplingMethod:
    """Enumeration for coupling method."""

    # Method resolution order:
    #     CouplingMethod
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Lagrangian = 3

    Nitsche = 1

    Penalization = 2

    Undefined = 0


# class CouplingZonePairing in libaster


class CouplingZonePairing(DataStructure):
    """Object to create contact pairing."""

    # Method resolution order:
    #     CouplingZonePairing
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.CouplingZonePairing, arg0: libaster.BaseMesh, arg1: int) -> None

        2. __init__(self: libaster.CouplingZonePairing, arg0: libaster.BaseMesh) -> None
        """

    def check(self, model):
        """Check common nodes and normals.

        Arguments:
            model [Model]: model.
        """

    def setCoefficient(self, coef_pena):
        """Set penalization's coefficient.

        Arguments:
            coef_pena [float]: penalization's coefficient.
        """

    def setMasterGroupsOfCells(self, groups_name):
        """Set master's side.

        Arguments:
            groups_name [list[str]]: list of groups.
        """

    def setMethod(self, method):
        """Set method.

        Returns:
            method [CouplingMethod]: method.
        """

    def setPairingParameters(self, parameters):
        """Set pairing parameters.

        Arguments:
            parameters [PairingParameter]: PairingParameterPtr.
        """

    def setSlaveGroupsOfCells(self, groups_name):
        """Set slave's side.

        Arguments:
            groups_name [list[str]]: list of groups.
        """

    def setVerbosity(self, verbosity):
        """Set verbosity.

        Arguments:
            verbosity [float]: verbosity level.
        """


# class CouplingPairing in libaster


class CouplingPairing(DataStructure):
    """Object to create contact pairing."""

    # Method resolution order:
    #     CouplingPairing
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0, arg1):
        pass

    def addZone(self, zone):
        """Add a new zone of coupling;

        Argument:
            zone [CouplingZonePairing]: zone.
        """

    def compute(self):
        """Compute the pairing quantitie

        Returns:
            bool: True if the pairing quantities are updated appropriately
        """

    def getFiniteElementDescriptor(self):
        """Return Finite Element Descriptor for virtual cells from pairing.

        Returns:
            FiniteElementDescriptor: finite element for virtual cells
        """

    def getHHOProjector(self, pairing):
        """Get HHO L2-projector

        Argument:
            pairing [FieldOnCellsReal]: pairing.

        Returns:
            FieldOnCellsReal: L2-projector.
        """

    def getMesh(self):
        """Mesh

        Returns:
            BaseMesh: the mesh
        """

    def getNumberOfPairs(self):
        """Return number of pairs

        Returns:
            int: number of pairs
        """

    def getPairingField(self):
        """Get intersection points

        Returns:
            FieldOnCellsReal: intersection points.
        """


# class BaseAssemblyMatrix in libaster


class BaseAssemblyMatrix(DataStructure):
    pass

    # Method resolution order:
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.BaseAssemblyMatrix, arg0: str) -> None

        2. __init__(self: libaster.BaseAssemblyMatrix, arg0: str, arg1: str) -> None

        3. __init__(self: libaster.BaseAssemblyMatrix, arg0: PhysicalProblem, arg1: str) -> None
        """

    def getCalculOption(self):
        """Return the option of CALCUL.

        Returns:
            str: Name of the option.
        """

    def getDOFNumbering(self):
        pass

    def getDirichletBCDOFs(self):
        """Return a vector with DOFs eliminated by Dirichlet boundaries conditions (if it exists)

        Returns:
            tuple(int): a vector with DOFs eliminated by Dirichlet boundaries conditions of
                size = neq + 1,
                tuple(ieq = 0, neq - 1) = 1 then DOF eliminated else 0,
                tuple(neq) = number of DOFs eliminated.
        """

    def getLagrangeScaling(self):
        """Return the scaling used for Lagrange multipliers. It returns 1 if no Lagrange.

        Returns:
            float: scaling used for Lagrange multipliers. It returns 1 if no Lagrange
            are present.
        """

    def getMesh(self):
        """Return the mesh.

        Returns:
            Mesh: a pointer to the mesh
        """

    def hasDirichletEliminationDOFs(self):
        """Tell if matrix has some DOFs eliminated by Dirichlet boundaries conditions.

        Returns:
            bool: *True* if matrix has some DOFs eliminated by Dirichlet boundaries conditions else *False*
        """

    def isBuilt(self):
        """Tell if the matrix has already been built.

        Returns:
            bool: *True* if the matrix has been built.
        """

    def isFactorized(self):
        """Tell if the matrix is factorized.

        Returns:
            bool: *True* if the matrix is factorized, *False* otherwise.
        """

    def isSymmetric(self):
        """Return True if matrix is symmetric"""

    def print(self, format="ASTER", unit=6):
        """Print the matrix in code_aster or matlab format (with information on the DOF).

        Arguments:
            format (str): 'ASTER' (default) or 'MATLAB'
        """

    def setDOFNumbering(self, arg0):
        pass

    def setSolverName(self, arg0):
        pass

    def symmetrize(self):
        """Make the assembly matrix symmetric in place"""

    def transpose(self):
        """Return the transpose of the matrix in place."""

    def updateDOFNumbering(self):
        pass

    def zeroDirichletBCDOFs(self, arg0):
        """Zero the Dirichlet BC DOFs in the given field

        Arguments:
            field (FieldOnNodesReal): the field to modify
        """


# class AssemblyMatrixDisplacementReal in libaster


class AssemblyMatrixDisplacementReal(BaseAssemblyMatrix):
    pass

    # Method resolution order:
    #     AssemblyMatrixDisplacementReal
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __call__(self, i_row, i_col):
        """Return one matrix coefficient.
        Only supports Morse storage. A coefficient absent from the profile is zero.
        UpdateValuePointers() has to be called first.

        Arguments:
            i_row [int] Local row DOF index, zero-based.
            j_col [int] Local column DOF index, zero-based.

        Returns:
            [float]: value of the matritx at (i_row, j_col)
        """

    def __iadd__(self, arg0):
        pass

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AssemblyMatrixDisplacementReal) -> None

        2. __init__(self: libaster.AssemblyMatrixDisplacementReal, arg0: str) -> None

        3. __init__(self: libaster.AssemblyMatrixDisplacementReal, arg0: PhysicalProblem) -> None

        4. __init__(self: libaster.AssemblyMatrixDisplacementReal, arg0: libaster.AssemblyMatrixDisplacementReal) -> None
        """

    def __isub__(self, arg0):
        pass

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def applyDirichletBC(self, arg0, arg1):
        """Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

        Arguments:
            DirichletBC [FieldOnNodes] the values on the DirichletBC.
            Rhs [FieldOnNodes] The residual to be modified.
        """

    def assemble(self, *args, **kwargs):
        """Overloaded function.

        1. assemble(self: libaster.AssemblyMatrixDisplacementReal, elemMatrix: ElementaryMatrix<double, (PhysicalQuantityEnum)4>, listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        2. assemble(self: libaster.AssemblyMatrixDisplacementReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)4>], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        3. assemble(self: libaster.AssemblyMatrixDisplacementReal, elemMatrix: list[tuple[ElementaryMatrix<double, (PhysicalQuantityEnum)4>, float]], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        4. assemble(self: libaster.AssemblyMatrixDisplacementReal, elemMatrix: ElementaryMatrix<double, (PhysicalQuantityEnum)4>, dirichlet: DirichletBC) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        5. assemble(self: libaster.AssemblyMatrixDisplacementReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)4>], dirichlet: DirichletBC) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        6. assemble(self: libaster.AssemblyMatrixDisplacementReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)4>], dirichlet: list[DirichletBC]) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (list[DirichletBC]) : dirichlet BC to impose.


        7. assemble(self: libaster.AssemblyMatrixDisplacementReal, elemMatrix: list[tuple[ElementaryMatrix<double, (PhysicalQuantityEnum)4>, float]], dirichlet: DirichletBC) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.
        """

    def copy(self):
        pass

    def defineSolver(self):
        pass

    def extractSubMatrix(self, fed, option, parameter):
        """Extract a sub-matrix.

        Arguments:
            fed [FiniteElementDescriptor] extract sub-matrix for element of fed.
            option [str] Option to known which dofs to extract.
            parameter [str] parameter of the option to known which dofs to extract.

        Returns:
            [ElementaryMatrixReal]: elementary matrix extracted.
        """

    def getLowerValues(self):
        pass

    def getUpperValues(self):
        pass

    def scale(self, arg0, arg1):
        """Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors

        Arguments:
            lvect (list[float]): List of the values.
            rvect (list[float]): List of the values.
        """

    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.

        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.

        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """

    def size(self, local=True):
        """Get the size of the matrix

        Arguments:
            local (bool) local or global size
        """

    def updateValuePointers(self):
        """Update pointers."""


# class AssemblyMatrixEliminatedReal in libaster


class AssemblyMatrixEliminatedReal(AssemblyMatrixDisplacementReal):
    pass

    # Method resolution order:
    #     AssemblyMatrixEliminatedReal
    #     AssemblyMatrixDisplacementReal
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AssemblyMatrixEliminatedReal) -> None

        2. __init__(self: libaster.AssemblyMatrixEliminatedReal, arg0: str) -> None
        """


# class AssemblyMatrixDisplacementComplex in libaster


class AssemblyMatrixDisplacementComplex(BaseAssemblyMatrix):
    pass

    # Method resolution order:
    #     AssemblyMatrixDisplacementComplex
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __call__(self, i_row, i_col):
        """Return one matrix coefficient.
        Only supports Morse storage. A coefficient absent from the profile is zero.
        UpdateValuePointers() has to be called first.

        Arguments:
            i_row [int] Local row DOF index, zero-based.
            j_col [int] Local column DOF index, zero-based.

        Returns:
            [cmplx]: value of the matritx at (i_row, j_col)
        """

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AssemblyMatrixDisplacementComplex) -> None

        2. __init__(self: libaster.AssemblyMatrixDisplacementComplex, arg0: str) -> None
        """

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def assemble(self, *args, **kwargs):
        """Overloaded function.

        1. assemble(self: libaster.AssemblyMatrixDisplacementComplex, elemMatrix: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>, listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        2. assemble(self: libaster.AssemblyMatrixDisplacementComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        3. assemble(self: libaster.AssemblyMatrixDisplacementComplex, elemMatrix: list[tuple[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>, float]], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        4. assemble(self: libaster.AssemblyMatrixDisplacementComplex, elemMatrix: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>, dirichlet: DirichletBC) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        5. assemble(self: libaster.AssemblyMatrixDisplacementComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>], dirichlet: DirichletBC) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        6. assemble(self: libaster.AssemblyMatrixDisplacementComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>], dirichlet: list[DirichletBC]) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (list[DirichletBC]) : dirichlet BC to impose.


        7. assemble(self: libaster.AssemblyMatrixDisplacementComplex, elemMatrix: list[tuple[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>, float]], dirichlet: DirichletBC) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.
        """

    def copy(self):
        pass

    def defineSolver(self):
        pass

    def extractSubMatrix(self, fed, option, parameter):
        """Extract a sub-matrix.

        Arguments:
            fed [FiniteElementDescriptor] extract sub-matrix for element of fed.
            option [str] Option to known which dofs to extract.
            parameter [str] parameter of the option to known which dofs to extract.

        Returns:
            [ElementaryMatrixComplex]: elementary matrix extracted.
        """

    def getLowerValues(self):
        pass

    def getUpperValues(self):
        pass

    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.

        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.

        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """

    def transposeConjugate(self):
        pass

    def updateValuePointers(self):
        """Update pointers."""


# class AssemblyMatrixTemperatureReal in libaster


class AssemblyMatrixTemperatureReal(BaseAssemblyMatrix):
    pass

    # Method resolution order:
    #     AssemblyMatrixTemperatureReal
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __call__(self, i_row, i_col):
        """Return one matrix coefficient.
        Only supports Morse storage. A coefficient absent from the profile is zero.
        UpdateValuePointers() has to be called first.

        Arguments:
            i_row [int] Local row DOF index, zero-based.
            j_col [int] Local column DOF index, zero-based.

        Returns:
            [float]: value of the matritx at (i_row, j_col)
        """

    def __iadd__(self, arg0):
        pass

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AssemblyMatrixTemperatureReal) -> None

        2. __init__(self: libaster.AssemblyMatrixTemperatureReal, arg0: str) -> None

        3. __init__(self: libaster.AssemblyMatrixTemperatureReal, arg0: PhysicalProblem) -> None
        """

    def __isub__(self, arg0):
        pass

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def applyDirichletBC(self, arg0, arg1):
        """Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

        Arguments:
            DirichletBC [FieldOnNodes] the values on the DirichletBC.
            Rhs [FieldOnNodes] The residual to be modified.
        """

    def assemble(self, *args, **kwargs):
        """Overloaded function.

        1. assemble(self: libaster.AssemblyMatrixTemperatureReal, elemMatrix: ElementaryMatrix<double, (PhysicalQuantityEnum)6>, listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        2. assemble(self: libaster.AssemblyMatrixTemperatureReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)6>], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        3. assemble(self: libaster.AssemblyMatrixTemperatureReal, elemMatrix: list[tuple[ElementaryMatrix<double, (PhysicalQuantityEnum)6>, float]], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        4. assemble(self: libaster.AssemblyMatrixTemperatureReal, elemMatrix: ElementaryMatrix<double, (PhysicalQuantityEnum)6>, dirichlet: DirichletBC) -> None


                       Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        5. assemble(self: libaster.AssemblyMatrixTemperatureReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)6>], dirichlet: DirichletBC) -> None


                         Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        6. assemble(self: libaster.AssemblyMatrixTemperatureReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)6>], dirichlet: list[DirichletBC]) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (list[DirichletBC]) : dirichlet BC to impose.


        7. assemble(self: libaster.AssemblyMatrixTemperatureReal, elemMatrix: list[tuple[ElementaryMatrix<double, (PhysicalQuantityEnum)6>, float]], dirichlet: DirichletBC) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.
        """

    def copy(self):
        pass

    def defineSolver(self):
        pass

    def extractSubMatrix(self, fed, option, parameter):
        """Extract a sub-matrix.

        Arguments:
            fed [FiniteElementDescriptor] extract sub-matrix for element of fed.
            option [str] Option to known which dofs to extract.
            parameter [str] parameter of the option to known which dofs to extract.

        Returns:
            [ElementaryMatrixReal]: elementary matrix extracted.
        """

    def getLowerValues(self):
        pass

    def getUpperValues(self):
        pass

    def scale(self, arg0, arg1):
        """Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors

        Arguments:
            lvect (list[float]): List of the values.
            rvect (list[float]): List of the values.
        """

    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.

        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.

        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """

    def size(self, local=True):
        """Get the size of the matrix

        Arguments:
            local (bool) local or global size
        """

    def updateValuePointers(self):
        """Update pointers."""


# class AssemblyMatrixTemperatureComplex in libaster


class AssemblyMatrixTemperatureComplex(BaseAssemblyMatrix):
    pass

    # Method resolution order:
    #     AssemblyMatrixTemperatureComplex
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AssemblyMatrixTemperatureComplex) -> None

        2. __init__(self: libaster.AssemblyMatrixTemperatureComplex, arg0: str) -> None
        """

    def assemble(self, *args, **kwargs):
        """Overloaded function.

        1. assemble(self: libaster.AssemblyMatrixTemperatureComplex, elemMatrix: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)6>, listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        2. assemble(self: libaster.AssemblyMatrixTemperatureComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)6>], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices added.

                        Arguments:
                            clean (bool) : Clean elementary matrices after building (default = true)


        3. assemble(self: libaster.AssemblyMatrixTemperatureComplex, elemMatrix: list[tuple[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)6>, float]], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        4. assemble(self: libaster.AssemblyMatrixTemperatureComplex, elemMatrix: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)6>, dirichlet: DirichletBC) -> None


                       Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        5. assemble(self: libaster.AssemblyMatrixTemperatureComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)6>], dirichlet: DirichletBC) -> None


                         Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        6. assemble(self: libaster.AssemblyMatrixTemperatureComplex, elemMatrix: list[tuple[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)6>, float]], dirichlet: DirichletBC) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.
        """

    def transposeConjugate(self):
        pass


# class AssemblyMatrixPressureReal in libaster


class AssemblyMatrixPressureReal(BaseAssemblyMatrix):
    pass

    # Method resolution order:
    #     AssemblyMatrixPressureReal
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __iadd__(self, arg0):
        pass

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AssemblyMatrixPressureReal) -> None

        2. __init__(self: libaster.AssemblyMatrixPressureReal, arg0: str) -> None
        """

    def __isub__(self, arg0):
        pass

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def applyDirichletBC(self, arg0, arg1):
        """Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

        Arguments:
            DirichletBC [FieldOnNodes] the values on the DirichletBC.
            Rhs [FieldOnNodes] The residual to be modified.
        """

    def assemble(self, *args, **kwargs):
        """Overloaded function.

        1. assemble(self: libaster.AssemblyMatrixPressureReal, elemMatrix: ElementaryMatrix<double, (PhysicalQuantityEnum)5>, listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        2. assemble(self: libaster.AssemblyMatrixPressureReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)5>], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        3. assemble(self: libaster.AssemblyMatrixPressureReal, elemMatrix: list[tuple[ElementaryMatrix<double, (PhysicalQuantityEnum)5>, float]], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        4. assemble(self: libaster.AssemblyMatrixPressureReal, elemMatrix: ElementaryMatrix<double, (PhysicalQuantityEnum)5>, dirichlet: DirichletBC) -> None


                       Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        5. assemble(self: libaster.AssemblyMatrixPressureReal, elemMatrix: list[ElementaryMatrix<double, (PhysicalQuantityEnum)5>], dirichlet: DirichletBC) -> None


                         Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        6. assemble(self: libaster.AssemblyMatrixPressureReal, elemMatrix: list[tuple[ElementaryMatrix<double, (PhysicalQuantityEnum)5>, float]], dirichlet: DirichletBC) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.
        """

    def copy(self):
        pass

    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.

        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.

        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """


# class AssemblyMatrixPressureComplex in libaster


class AssemblyMatrixPressureComplex(BaseAssemblyMatrix):
    pass

    # Method resolution order:
    #     AssemblyMatrixPressureComplex
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, arg0):
        pass

    def __iadd__(self, arg0):
        pass

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AssemblyMatrixPressureComplex) -> None

        2. __init__(self: libaster.AssemblyMatrixPressureComplex, arg0: str) -> None
        """

    def __isub__(self, arg0):
        pass

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def assemble(self, *args, **kwargs):
        """Overloaded function.

        1. assemble(self: libaster.AssemblyMatrixPressureComplex, elemMatrix: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>, listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        2. assemble(self: libaster.AssemblyMatrixPressureComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        3. assemble(self: libaster.AssemblyMatrixPressureComplex, elemMatrix: list[tuple[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>, float]], listOfLoads: libaster.ListOfLoads = None) -> None


                        Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            listOfLoads (ListOfLoads) : list of loads to assemble


        4. assemble(self: libaster.AssemblyMatrixPressureComplex, elemMatrix: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>, dirichlet: DirichletBC) -> None


                       Assembly matrix from elementar matrices and list of loads.

                        Arguments:
                            elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        5. assemble(self: libaster.AssemblyMatrixPressureComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>], dirichlet: DirichletBC) -> None


                         Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        6. assemble(self: libaster.AssemblyMatrixPressureComplex, elemMatrix: list[tuple[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>, float]], dirichlet: DirichletBC) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                            elementary matrix and the multiplicatif coefficent to assemble.
                            dirichlet (DirichletBC) : dirichlet BC to impose.


        7. assemble(self: libaster.AssemblyMatrixPressureComplex, elemMatrix: list[ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>], dirichlet: list[DirichletBC]) -> None


                        Arguments:
                            elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                            dirichlet (list[DirichletBC]) : dirichlet BC to impose.
        """

    def copy(self):
        pass

    def defineSolver(self):
        pass

    def getLowerValues(self):
        pass

    def getUpperValues(self):
        pass

    def transposeConjugate(self):
        pass


# class ElementaryTermReal in libaster


class ElementaryTermReal(DataField):
    pass

    # Method resolution order:
    #     ElementaryTermReal
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryTermReal) -> None

        2. __init__(self: libaster.ElementaryTermReal, arg0: str) -> None
        """

    def exists(self):
        """Return True if it exists

        Returns:
            bool: True if exist
        """

    def getFiniteElementDescriptor(self):
        """Return the finite element descriptor

        Returns:
            FiniteElementDescriptor: finite element descriptor
        """

    def getLocalMode(self):
        """Return the local mode.

        Returns:
            str: the local mode
        """

    def getMesh(self):
        """Return the mesh

        Returns:
            BaseMesh: a pointer to the mesh
        """

    def getOption(self):
        """Return the optior used to compute it

        Returns:
            str: name of the option
        """

    def getPhysicalQuantity(self):
        """Return the physical quantity

        Returns:
            str: name of the physical quantity
        """

    def getValues(self):
        """Return the values of the field.

        Returns:
            list[list[float]]: values
        """


# class ElementaryTermComplex in libaster


class ElementaryTermComplex(DataField):
    pass

    # Method resolution order:
    #     ElementaryTermComplex
    #     DataField
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryTermComplex) -> None

        2. __init__(self: libaster.ElementaryTermComplex, arg0: str) -> None
        """

    def exists(self):
        """Return True if it exists

        Returns:
            bool: True if exist
        """

    def getFiniteElementDescriptor(self):
        """Return the finite element descriptor

        Returns:
            FiniteElementDescriptor: finite element descriptor
        """

    def getLocalMode(self):
        """Return the local mode.

        Returns:
            str: the local mode
        """

    def getMesh(self):
        """Return the mesh

        Returns:
            BaseMesh: a pointer to the mesh
        """

    def getOption(self):
        """Return the optior used to compute it

        Returns:
            str: name of the option
        """

    def getPhysicalQuantity(self):
        """Return the physical quantity

        Returns:
            str: name of the physical quantity
        """


# class BaseElementaryMatrix in libaster


class BaseElementaryMatrix(DataStructure):
    pass

    # Method resolution order:
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getModel(self):
        pass

    def setModel(self, arg0):
        pass


# class ElementaryMatrixDisplacementReal in libaster


class ElementaryMatrixDisplacementReal(BaseElementaryMatrix):
    pass

    # Method resolution order:
    #     ElementaryMatrixDisplacementReal
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryMatrixDisplacementReal) -> None

        2. __init__(self: libaster.ElementaryMatrixDisplacementReal, arg0: str) -> None

        3. __init__(self: libaster.ElementaryMatrixDisplacementReal, arg0: Model, arg1: str) -> None
        """

    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.

        1. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementReal, arg0: libaster.ElementaryTermReal) -> None

        2. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementReal, arg0: list[libaster.ElementaryTermReal]) -> None
        """

    def build(self):
        pass

    def getElementaryTerms(self):
        pass

    def getFiniteElementDescriptors(self):
        pass

    def getMesh(self):
        pass

    def getNumberOfElementaryTerms(self):
        pass

    def hasElementaryTerms(self):
        pass


# class ElementaryMatrixDisplacementComplex in libaster


class ElementaryMatrixDisplacementComplex(BaseElementaryMatrix):
    pass

    # Method resolution order:
    #     ElementaryMatrixDisplacementComplex
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryMatrixDisplacementComplex) -> None

        2. __init__(self: libaster.ElementaryMatrixDisplacementComplex, arg0: str) -> None

        3. __init__(self: libaster.ElementaryMatrixDisplacementComplex, arg0: Model, arg1: str) -> None
        """

    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.

        1. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementComplex, arg0: libaster.ElementaryTermComplex) -> None

        2. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementComplex, arg0: list[libaster.ElementaryTermComplex]) -> None
        """

    def build(self):
        pass

    def getElementaryTerms(self):
        pass

    def getFiniteElementDescriptors(self):
        pass

    def getMesh(self):
        pass

    def getNumberOfElementaryTerms(self):
        pass

    def hasElementaryTerms(self):
        pass


# class ElementaryMatrixTemperatureReal in libaster


class ElementaryMatrixTemperatureReal(BaseElementaryMatrix):
    pass

    # Method resolution order:
    #     ElementaryMatrixTemperatureReal
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __imul__(self, arg0):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryMatrixTemperatureReal) -> None

        2. __init__(self: libaster.ElementaryMatrixTemperatureReal, arg0: str) -> None

        3. __init__(self: libaster.ElementaryMatrixTemperatureReal, arg0: Model, arg1: str) -> None
        """

    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.

        1. addElementaryTerm(self: libaster.ElementaryMatrixTemperatureReal, arg0: libaster.ElementaryTermReal) -> None

        2. addElementaryTerm(self: libaster.ElementaryMatrixTemperatureReal, arg0: list[libaster.ElementaryTermReal]) -> None
        """

    def build(self):
        pass

    def getElementaryTerms(self):
        pass

    def getFiniteElementDescriptors(self):
        pass

    def getMesh(self):
        pass

    def getNumberOfElementaryTerms(self):
        pass

    def hasElementaryTerms(self):
        pass


# class ElementaryMatrixPressureComplex in libaster


class ElementaryMatrixPressureComplex(BaseElementaryMatrix):
    pass

    # Method resolution order:
    #     ElementaryMatrixPressureComplex
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryMatrixPressureComplex) -> None

        2. __init__(self: libaster.ElementaryMatrixPressureComplex, arg0: str) -> None

        3. __init__(self: libaster.ElementaryMatrixPressureComplex, arg0: Model, arg1: str) -> None
        """

    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.

        1. addElementaryTerm(self: libaster.ElementaryMatrixPressureComplex, arg0: libaster.ElementaryTermComplex) -> None

        2. addElementaryTerm(self: libaster.ElementaryMatrixPressureComplex, arg0: list[libaster.ElementaryTermComplex]) -> None
        """

    def build(self):
        pass

    def getElementaryTerms(self):
        pass

    def getFiniteElementDescriptors(self):
        pass

    def getMesh(self):
        pass

    def getNumberOfElementaryTerms(self):
        pass

    def hasElementaryTerms(self):
        pass


# class BaseElementaryVector in libaster


class BaseElementaryVector(DSWithCppPickling):
    pass

    # Method resolution order:
    #     BaseElementaryVector
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.BaseElementaryVector, arg0: str, arg1: str, arg2: Model) -> None

        2. __init__(self: libaster.BaseElementaryVector, arg0: Model) -> None
        """

    def addSubstructuring(self, arg0):
        pass

    def assembleWithLoadFunctions(self, dofNume, loads, time=0.0):
        pass

    def assembleWithMask(self, arg0, arg1, arg2):
        pass

    def build(self, FED=[]):
        pass


# class ElementaryVectorReal in libaster


class ElementaryVectorReal(BaseElementaryVector):
    pass

    # Method resolution order:
    #     ElementaryVectorReal
    #     BaseElementaryVector
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryVectorReal, arg0: str, arg1: str, arg2: Model) -> None

        2. __init__(self: libaster.ElementaryVectorReal, arg0: Model) -> None
        """

    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.

        1. addElementaryTerm(self: libaster.ElementaryVectorReal, term: libaster.ElementaryTermReal) -> None


                    Add elementary term

                    Arguments:
                        term (ElementaryTermReal): elementary term


        2. addElementaryTerm(self: libaster.ElementaryVectorReal, terms: list[libaster.ElementaryTermReal]) -> None


                    Add vector of elementary term

                    Arguments:
                        terms (list[ElementaryTermReal]): vector of elementary term
        """

    def assemble(self, dofNume, minimum=False):
        pass

    def getElementaryTerms(self):
        pass

    def getFiniteElementDescriptor(self):
        pass

    def getNumberOfElementaryTerms(self):
        pass

    def getVeass(self):
        pass

    def setVeass(self, arg0, arg1):
        pass


# class ElementaryVectorComplex in libaster


class ElementaryVectorComplex(BaseElementaryVector):
    pass

    # Method resolution order:
    #     ElementaryVectorComplex
    #     BaseElementaryVector
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryVectorComplex, arg0: str, arg1: str, arg2: Model) -> None

        2. __init__(self: libaster.ElementaryVectorComplex, arg0: Model) -> None
        """

    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.

        1. addElementaryTerm(self: libaster.ElementaryVectorComplex, term: libaster.ElementaryTermComplex) -> None


                    Add elementary term

                    Arguments:
                        term (ElementaryTermComplex): elementary term


        2. addElementaryTerm(self: libaster.ElementaryVectorComplex, terms: list[libaster.ElementaryTermComplex]) -> None


                    Add vector of elementary term

                    Arguments:
                        terms (list[ElementaryTermComplex]): vector of elementary term
        """

    def assemble(self, dofNume, minimum=False):
        pass

    def getElementaryTerms(self):
        pass

    def getNumberOfElementaryTerms(self):
        pass

    def getVeass(self):
        pass

    def setVeass(self, arg0, arg1):
        pass


# class ElementaryVectorDisplacementReal in libaster


class ElementaryVectorDisplacementReal(ElementaryVectorReal):
    pass

    # Method resolution order:
    #     ElementaryVectorDisplacementReal
    #     ElementaryVectorReal
    #     BaseElementaryVector
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryVectorDisplacementReal) -> None

        2. __init__(self: libaster.ElementaryVectorDisplacementReal, arg0: tuple) -> None

        3. __init__(self: libaster.ElementaryVectorDisplacementReal, arg0: str, arg1: Model) -> None

        4. __init__(self: libaster.ElementaryVectorDisplacementReal, arg0: Model) -> None
        """

    def __setstate__(self, arg0):
        pass


# class ElementaryVectorTemperatureReal in libaster


class ElementaryVectorTemperatureReal(ElementaryVectorReal):
    pass

    # Method resolution order:
    #     ElementaryVectorTemperatureReal
    #     ElementaryVectorReal
    #     BaseElementaryVector
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryVectorTemperatureReal) -> None

        2. __init__(self: libaster.ElementaryVectorTemperatureReal, arg0: tuple) -> None

        3. __init__(self: libaster.ElementaryVectorTemperatureReal, arg0: str, arg1: Model) -> None

        4. __init__(self: libaster.ElementaryVectorTemperatureReal, arg0: Model) -> None
        """

    def __setstate__(self, arg0):
        pass


# class ElementaryVectorPressureComplex in libaster


class ElementaryVectorPressureComplex(ElementaryVectorComplex):
    pass

    # Method resolution order:
    #     ElementaryVectorPressureComplex
    #     ElementaryVectorComplex
    #     BaseElementaryVector
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElementaryVectorPressureComplex) -> None

        2. __init__(self: libaster.ElementaryVectorPressureComplex, arg0: tuple) -> None

        3. __init__(self: libaster.ElementaryVectorPressureComplex, arg0: str, arg1: Model) -> None

        4. __init__(self: libaster.ElementaryVectorPressureComplex, arg0: Model) -> None
        """

    def __setstate__(self, arg0):
        pass


# class GeneralizedAssemblyMatrix in libaster


class GeneralizedAssemblyMatrix(DataStructure):
    pass

    # Method resolution order:
    #     GeneralizedAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def exists(self):
        """Return True if the matrix exists

        Returns:
            bool: True if the matrix exists else False.
        """

    def getGeneralizedDOFNumbering(self):
        pass

    def getModalBasis(self):
        pass

    def isDense(self):
        """Return True if the matrix is dense

        Returns:
            bool: True if the matrix is dense else False.
        """

    def isDiagonal(self):
        """Return True if the matrix is diagonal

        Returns:
            bool: True if the matrix is diagonal else False.
        """

    def setGeneralizedDOFNumbering(self, arg0):
        pass

    def setModalBasis(self, *args, **kwargs):
        """Overloaded function.

        1. setModalBasis(self: libaster.GeneralizedAssemblyMatrix, arg0: GeneralizedModeResult) -> bool

        2. setModalBasis(self: libaster.GeneralizedAssemblyMatrix, arg0: ModeResult) -> bool
        """

    def size(self):
        """Return the size of the matrix

        Returns:
            int: size of the matrix.
        """


# class GeneralizedAssemblyMatrixReal in libaster


class GeneralizedAssemblyMatrixReal(GeneralizedAssemblyMatrix):
    pass

    # Method resolution order:
    #     GeneralizedAssemblyMatrixReal
    #     GeneralizedAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedAssemblyMatrixReal) -> None

        2. __init__(self: libaster.GeneralizedAssemblyMatrixReal, arg0: str) -> None
        """

    def allocate(self, isSymmetric=True):
        """Allocate the matrix"""

    def getLowerValues(self):
        """Return the lower part of the matrix.

        Returns:
            list[float]: lower part of the matrix.
        """

    def getUpperValues(self):
        """Return the upper part of the matrix.

        Returns:
            list[float]: upper part of the matrix.
        """

    def isSymmetric(self):
        """Return True if the matrix is symmetric

        Returns:
            bool: True if the matrix is symmetric else False.
        """

    def setLowerValues(self, values):
        """Set the lower part of the matrix.

        Arguments:
            values [list[float]]: set lower part of the matrix.
        """

    def setUpperValues(self, values):
        """Set the upper part of the matrix.

        Arguments:
            values [list[float]]: set upper part of the matrix.
        """


# class GeneralizedAssemblyMatrixComplex in libaster


class GeneralizedAssemblyMatrixComplex(GeneralizedAssemblyMatrix):
    pass

    # Method resolution order:
    #     GeneralizedAssemblyMatrixComplex
    #     GeneralizedAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedAssemblyMatrixComplex) -> None

        2. __init__(self: libaster.GeneralizedAssemblyMatrixComplex, arg0: str) -> None
        """

    def getLowerValues(self):
        """Return the lower part of the matrix.

        Returns:
            list[complex]: lower part of the matrix.
        """

    def getUpperValues(self):
        """Return the upper part of the matrix.

        Returns:
            list[complex]: upper part of the matrix.
        """

    def isSymmetric(self):
        """Return True if the matrix is symmetric

        Returns:
            bool: True if the matrix is symmetric else False.
        """

    def setLowerValues(self, values):
        """Set the lower part of the matrix.

        Arguments:
            values [list[complex]]: set lower part of the matrix.
        """

    def setUpperValues(self, values):
        """Set the upper part of the matrix.

        Arguments:
            values [list[complex]]: set upper part of the matrix.
        """


# class GeneralizedAssemblyVector in libaster


class GeneralizedAssemblyVector(DataStructure):
    pass

    # Method resolution order:
    #     GeneralizedAssemblyVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""


# class GeneralizedAssemblyVectorReal in libaster


class GeneralizedAssemblyVectorReal(GeneralizedAssemblyVector):
    pass

    # Method resolution order:
    #     GeneralizedAssemblyVectorReal
    #     GeneralizedAssemblyVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedAssemblyVectorReal) -> None

        2. __init__(self: libaster.GeneralizedAssemblyVectorReal, arg0: str) -> None
        """

    def getValues(self):
        """Return a list of values as (x1, y1, z1, x2, y2, z2...)

        Returns:
            list[float]: List of values.
        """

    def setValues(self, arg0):
        """Set values of vector.

        Arguments:
            values (list[float]): set vector.
        """


# class GeneralizedAssemblyVectorComplex in libaster


class GeneralizedAssemblyVectorComplex(GeneralizedAssemblyVector):
    pass

    # Method resolution order:
    #     GeneralizedAssemblyVectorComplex
    #     GeneralizedAssemblyVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedAssemblyVectorComplex) -> None

        2. __init__(self: libaster.GeneralizedAssemblyVectorComplex, arg0: str) -> None
        """

    def getValues(self):
        """Return a list of values as (x1, y1, z1, x2, y2, z2...)

        Returns:
            list[complex]: List of values.
        """

    def setValues(self, arg0):
        """Set values of vector.

        Arguments:
            values (list[complex]): set vector.
        """


# class GeneralizedLoad in libaster


class GeneralizedLoad(DataStructure):
    pass

    # Method resolution order:
    #     GeneralizedLoad
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedLoad) -> None

        2. __init__(self: libaster.GeneralizedLoad, arg0: str) -> None
        """

    def getAssemblyVector(self):
        pass

    def getDOFNumbering(self):
        pass

    def getMPCs(self):
        pass

    def setDOFNumbering(self, arg0):
        pass

    def setMPCs(self, arg0):
        pass


# class InterspectralMatrix in libaster


class InterspectralMatrix(DataStructure):
    pass

    # Method resolution order:
    #     InterspectralMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.InterspectralMatrix) -> None

        2. __init__(self: libaster.InterspectralMatrix, arg0: str) -> None
        """

    def getColumnComponents(self):
        pass

    def getColumnIndexes(self):
        pass

    def getColumnNodes(self):
        pass

    def getLineComponents(self):
        pass

    def getLineIndexes(self):
        pass

    def getLineNodes(self):
        pass

    def getNumberOfFrequencies(self):
        pass


# class LinearSolver in libaster


class LinearSolver(DataStructure):
    pass

    # Method resolution order:
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def build(self):
        """build internal objects of the solver

        Returns:
             bool: True if the building is a success, else False
        """

    def deleteFactorizedMatrix(self):
        """delete the factorized matrix and its preconditionner if created.
        This is the case for Mumps and Petsc.

        Returns:
             bool: True if success, else False
        """

    def enableXfem(self):
        """Enable preconditionning for XFEM modeling."""

    def factorize(self, matrix, raiseException=False):
        """Factorize the matrix.

        Arguments:
            matrix (BaseAssemblyMatrix) : matrix to factorize
            raiseException (bool): if *True* an exception is raised in case of error,
            otherwise it stops with an error (default: *False*).

        Returns:
            bool: *True* if factorization is a success, else *False*
        """

    def getMatrix(self):
        """return the factorized matrix

        Returns:
            BaseAssemblyMatrix: factorized matrix
        """

    def getPrecondMatrix(self):
        """return the preconditionning matrix

        Returns:
            BaseAssemblyMatrix: preconditionning matrix
        """

    def getSolverName(self):
        """Get the name of the solver used between 'MUMPS', 'PETSC', 'MULT_FRONT' and 'PETSC'

        Returns:
             str: name of the solver used
        """

    def setCataPath(self, path):
        """Set the path of the catalog that defines the solver keywords.
        It can be command name or a path as *code_aster.Cata.Commons.xxxx*.

        Arguments:
             path (str): command name or path of the catalog.
        """

    def setKeywords(self, arg0):
        pass

    def solve(self, *args, **kwargs):
        """Overloaded function.

        1. solve(self: libaster.LinearSolver, rhs: libaster.FieldOnNodesReal, dirichletBC: libaster.FieldOnNodesReal = None) -> libaster.FieldOnNodesReal

        2. solve(self: libaster.LinearSolver, rhs: libaster.FieldOnNodesComplex, dirichletBC: libaster.FieldOnNodesComplex = None) -> libaster.FieldOnNodesComplex
        """

    def supportParallelMesh(self):
        """tell if the solver is enable in HPC

        Returns:
             bool: True if the solver support ParallelMesh, else False
        """


# class MultFrontSolver in libaster


class MultFrontSolver(LinearSolver):
    pass

    # Method resolution order:
    #     MultFrontSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MultFrontSolver) -> None

        2. __init__(self: libaster.MultFrontSolver, arg0: str) -> None
        """


# class LdltSolver in libaster


class LdltSolver(LinearSolver):
    pass

    # Method resolution order:
    #     LdltSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.LdltSolver) -> None

        2. __init__(self: libaster.LdltSolver, arg0: str) -> None
        """


# class MumpsSolver in libaster


class MumpsSolver(LinearSolver):
    pass

    # Method resolution order:
    #     MumpsSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MumpsSolver) -> None

        2. __init__(self: libaster.MumpsSolver, arg0: str) -> None
        """


# class PetscSolver in libaster


class PetscSolver(LinearSolver):
    pass

    # Method resolution order:
    #     PetscSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.PetscSolver) -> None

        2. __init__(self: libaster.PetscSolver, arg0: str) -> None
        """

    def getPetscOptions(self):
        """return the petsc solver options

        Returns:
            string: the petsc solver options
        """


# class GcpcSolver in libaster


class GcpcSolver(LinearSolver):
    pass

    # Method resolution order:
    #     GcpcSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GcpcSolver) -> None

        2. __init__(self: libaster.GcpcSolver, arg0: str) -> None
        """


# class GenericModalBasis in libaster


class GenericModalBasis(DataStructure):
    pass

    # Method resolution order:
    #     GenericModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""


# class StandardModalBasis in libaster


class StandardModalBasis(GenericModalBasis):
    pass

    # Method resolution order:
    #     StandardModalBasis
    #     GenericModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.StandardModalBasis) -> None

        2. __init__(self: libaster.StandardModalBasis, arg0: str) -> None
        """


# class RitzBasis in libaster


class RitzBasis(GenericModalBasis):
    pass

    # Method resolution order:
    #     RitzBasis
    #     GenericModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.RitzBasis) -> None

        2. __init__(self: libaster.RitzBasis, arg0: str) -> None
        """


# class InterfaceType in libaster


class InterfaceType:
    """Enumeration of interface type."""

    # Method resolution order:
    #     InterfaceType
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    CraigBampton = 1

    HarmonicCraigBampton = 2

    MacNeal = 0


# class StructureInterface in libaster


class StructureInterface(DataStructure):
    pass

    # Method resolution order:
    #     StructureInterface
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.StructureInterface) -> None

        2. __init__(self: libaster.StructureInterface, arg0: str) -> None

        3. __init__(self: libaster.StructureInterface, arg0: libaster.DOFNumbering) -> None

        4. __init__(self: libaster.StructureInterface, arg0: str, arg1: libaster.DOFNumbering) -> None
        """


# class AcousticLoadComplex in libaster


class AcousticLoadComplex(DataStructure):
    pass

    # Method resolution order:
    #     AcousticLoadComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AcousticLoadComplex, arg0: Model) -> None

        2. __init__(self: libaster.AcousticLoadComplex, arg0: str, arg1: Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getMesh(self):
        pass

    def getModel(self):
        pass


# class DirichletBC in libaster


class DirichletBC(DataStructure):
    pass

    # Method resolution order:
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def build(self):
        pass

    def getModel(self):
        """Return the model

        Returns:
            ModelPtr: a pointer to the model
        """

    def getPhysics(self):
        """To know the physics supported by the model

        Returns:
            str: Mechanics or Thermal or Acoustic
        """

    def setSyntax(self, syntax):
        """Function to set the syntax used to build object

        Arguments:
            syntax: the syntax
        """


# class MechanicalDirichletBC in libaster


class MechanicalDirichletBC(DirichletBC):
    pass

    # Method resolution order:
    #     MechanicalDirichletBC
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MechanicalDirichletBC, arg0: Model) -> None

        2. __init__(self: libaster.MechanicalDirichletBC, arg0: str, arg1: Model) -> None
        """

    def addBCOnCells(self, *args, **kwargs):
        """Overloaded function.

        1. addBCOnCells(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool

        2. addBCOnCells(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: list[str]) -> bool
        """

    def addBCOnNodes(self, *args, **kwargs):
        """Overloaded function.

        1. addBCOnNodes(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool

        2. addBCOnNodes(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: list[str]) -> bool
        """


# class ThermalDirichletBC in libaster


class ThermalDirichletBC(DirichletBC):
    pass

    # Method resolution order:
    #     ThermalDirichletBC
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ThermalDirichletBC, arg0: Model) -> None

        2. __init__(self: libaster.ThermalDirichletBC, arg0: str, arg1: Model) -> None
        """

    def addBCOnCells(self, *args, **kwargs):
        """Overloaded function.

        1. addBCOnCells(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool

        2. addBCOnCells(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: list[str]) -> bool
        """

    def addBCOnNodes(self, *args, **kwargs):
        """Overloaded function.

        1. addBCOnNodes(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool

        2. addBCOnNodes(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: list[str]) -> bool

        3. addBCOnNodes(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: libaster.Function, arg2: list[str]) -> bool
        """


# class AcousticDirichletBC in libaster


class AcousticDirichletBC(DirichletBC):
    pass

    # Method resolution order:
    #     AcousticDirichletBC
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AcousticDirichletBC, arg0: Model) -> None

        2. __init__(self: libaster.AcousticDirichletBC, arg0: str, arg1: Model) -> None
        """


# class MechanicalLoadReal in libaster


class MechanicalLoadReal(DSWithCppPickling, DataStructure):
    pass

    # Method resolution order:
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MechanicalLoadReal, arg0: Model) -> None

        2. __init__(self: libaster.MechanicalLoadReal, arg0: str, arg1: Model) -> None

        3. __init__(self: libaster.MechanicalLoadReal, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getFiniteElementDescriptor(self):
        pass

    def getMechanicalLoadDescription(self):
        pass

    def getMesh(self):
        pass

    def getModel(self):
        pass

    def getPairingField(self):
        """Return pairing intersection.

        Returns:
            FieldOnCellsReal: pairing intersection.
        """

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """

    def hasLoadField(self, arg0):
        pass

    def setPairingField(self, pairs):
        """Set pairing intersection.

        Arguments:
            pairs (FieldOnCellsReal): pairing intersection.
        """

    def updateValuePointers(self):
        pass


# class MechanicalLoadFunction in libaster


class MechanicalLoadFunction(DSWithCppPickling, DataStructure):
    pass

    # Method resolution order:
    #     MechanicalLoadFunction
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MechanicalLoadFunction, arg0: Model) -> None

        2. __init__(self: libaster.MechanicalLoadFunction, arg0: str, arg1: Model) -> None

        3. __init__(self: libaster.MechanicalLoadFunction, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getFiniteElementDescriptor(self):
        pass

    def getMesh(self):
        pass

    def getModel(self):
        pass

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """

    def hasLoadField(self, arg0):
        pass

    def updateValuePointers(self):
        pass


# class MechanicalLoadComplex in libaster


class MechanicalLoadComplex(DSWithCppPickling, DataStructure):
    pass

    # Method resolution order:
    #     MechanicalLoadComplex
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MechanicalLoadComplex, arg0: Model) -> None

        2. __init__(self: libaster.MechanicalLoadComplex, arg0: str, arg1: Model) -> None

        3. __init__(self: libaster.MechanicalLoadComplex, arg0: tuple) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getFiniteElementDescriptor(self):
        pass

    def getMesh(self):
        pass

    def getModel(self):
        pass

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """

    def hasLoadField(self, arg0):
        pass

    def updateValuePointers(self):
        pass


# class MechanicalLoadDescriptionReal in libaster


class MechanicalLoadDescriptionReal(DataStructure):
    pass

    # Method resolution order:
    #     MechanicalLoadDescriptionReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getConstantLoadField(self, arg0):
        pass


# class Loads in libaster


class Loads:
    """Enumeration for type of load."""

    # Method resolution order:
    #     Loads
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    DistributedPressure = 9

    ForceOnBeam = 5

    ForceOnEdge = 1

    ForceOnFace = 2

    ForceOnShell = 6

    ImposedDoF = 8

    InternalForce = 4

    LineicForce = 3

    NodalForce = 0

    NormalSpeedOnFace = 10

    PressureOnPipe = 7

    THMFlux = 12

    WavePressureOnFace = 11


# class NodalForceReal in libaster


class NodalForceReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     NodalForceReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.NodalForceReal, arg0: Model) -> None

        2. __init__(self: libaster.NodalForceReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class NodalStructuralForceReal in libaster


class NodalStructuralForceReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     NodalStructuralForceReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.NodalStructuralForceReal, arg0: Model) -> None

        2. __init__(self: libaster.NodalStructuralForceReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class ForceOnFaceReal in libaster


class ForceOnFaceReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     ForceOnFaceReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ForceOnFaceReal, arg0: Model) -> None

        2. __init__(self: libaster.ForceOnFaceReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class ForceOnEdgeReal in libaster


class ForceOnEdgeReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     ForceOnEdgeReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ForceOnEdgeReal, arg0: Model) -> None

        2. __init__(self: libaster.ForceOnEdgeReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class StructuralForceOnEdgeReal in libaster


class StructuralForceOnEdgeReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     StructuralForceOnEdgeReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.StructuralForceOnEdgeReal, arg0: Model) -> None

        2. __init__(self: libaster.StructuralForceOnEdgeReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class LineicForceReal in libaster


class LineicForceReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     LineicForceReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.LineicForceReal, arg0: Model) -> None

        2. __init__(self: libaster.LineicForceReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class InternalForceReal in libaster


class InternalForceReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     InternalForceReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.InternalForceReal, arg0: Model) -> None

        2. __init__(self: libaster.InternalForceReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class StructuralForceOnBeamReal in libaster


class StructuralForceOnBeamReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     StructuralForceOnBeamReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.StructuralForceOnBeamReal, arg0: Model) -> None

        2. __init__(self: libaster.StructuralForceOnBeamReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class LocalForceOnBeamReal in libaster


class LocalForceOnBeamReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     LocalForceOnBeamReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.LocalForceOnBeamReal, arg0: Model) -> None

        2. __init__(self: libaster.LocalForceOnBeamReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class StructuralForceOnShellReal in libaster


class StructuralForceOnShellReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     StructuralForceOnShellReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.StructuralForceOnShellReal, arg0: Model) -> None

        2. __init__(self: libaster.StructuralForceOnShellReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class LocalForceOnShellReal in libaster


class LocalForceOnShellReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     LocalForceOnShellReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.LocalForceOnShellReal, arg0: Model) -> None

        2. __init__(self: libaster.LocalForceOnShellReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class PressureOnShellReal in libaster


class PressureOnShellReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     PressureOnShellReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.PressureOnShellReal, arg0: Model) -> None

        2. __init__(self: libaster.PressureOnShellReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class PressureOnPipeReal in libaster


class PressureOnPipeReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     PressureOnPipeReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.PressureOnPipeReal, arg0: Model) -> None

        2. __init__(self: libaster.PressureOnPipeReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class ImposedDisplacementReal in libaster


class ImposedDisplacementReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     ImposedDisplacementReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ImposedDisplacementReal, arg0: Model) -> None

        2. __init__(self: libaster.ImposedDisplacementReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class ImposedPressureReal in libaster


class ImposedPressureReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     ImposedPressureReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ImposedPressureReal, arg0: Model) -> None

        2. __init__(self: libaster.ImposedPressureReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class DistributedPressureReal in libaster


class DistributedPressureReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     DistributedPressureReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.DistributedPressureReal, arg0: Model) -> None

        2. __init__(self: libaster.DistributedPressureReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class NormalSpeedOnFaceReal in libaster


class NormalSpeedOnFaceReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     NormalSpeedOnFaceReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.NormalSpeedOnFaceReal, arg0: Model) -> None

        2. __init__(self: libaster.NormalSpeedOnFaceReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class WavePressureOnFaceReal in libaster


class WavePressureOnFaceReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     WavePressureOnFaceReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.WavePressureOnFaceReal, arg0: Model) -> None

        2. __init__(self: libaster.WavePressureOnFaceReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class DistributedHeatFluxReal in libaster


class DistributedHeatFluxReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     DistributedHeatFluxReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.DistributedHeatFluxReal, arg0: Model) -> None

        2. __init__(self: libaster.DistributedHeatFluxReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class DistributedHydraulicFluxReal in libaster


class DistributedHydraulicFluxReal(MechanicalLoadReal):
    pass

    # Method resolution order:
    #     DistributedHydraulicFluxReal
    #     MechanicalLoadReal
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.DistributedHydraulicFluxReal, arg0: Model) -> None

        2. __init__(self: libaster.DistributedHydraulicFluxReal, arg0: str, arg1: Model) -> None
        """

    def build(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class PhysicalQuantityComponent in libaster


class PhysicalQuantityComponent:
    """Enumeration for physical component."""

    # Method resolution order:
    #     PhysicalQuantityComponent
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Drx = 3

    Dry = 4

    Drz = 5

    Dx = 0

    Dy = 1

    Dz = 2

    F1 = 21

    F2 = 22

    F3 = 23

    Flun = 28

    FlunHydr1 = 29

    FlunHydr2 = 30

    Fx = 9

    Fy = 10

    Fz = 11

    Impe = 26

    Mf1 = 24

    Mf2 = 25

    Mfy = 19

    Mfz = 20

    MiddleTemp = 7

    Mt = 18

    Mx = 12

    My = 13

    Mz = 14

    N = 15

    Pres = 8

    Temp = 6

    Vnor = 27

    Vy = 16

    Vz = 17


# class ForceReal in libaster


class ForceReal:
    pass

    # Method resolution order:
    #     ForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class StructuralForceReal in libaster


class StructuralForceReal:
    pass

    # Method resolution order:
    #     StructuralForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class LocalBeamForceReal in libaster


class LocalBeamForceReal:
    pass

    # Method resolution order:
    #     LocalBeamForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class LocalShellForceReal in libaster


class LocalShellForceReal:
    pass

    # Method resolution order:
    #     LocalShellForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class DisplacementReal in libaster


class DisplacementReal:
    pass

    # Method resolution order:
    #     DisplacementReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class PressureReal in libaster


class PressureReal:
    pass

    # Method resolution order:
    #     PressureReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class ImpedanceReal in libaster


class ImpedanceReal:
    pass

    # Method resolution order:
    #     ImpedanceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class NormalSpeedReal in libaster


class NormalSpeedReal:
    pass

    # Method resolution order:
    #     NormalSpeedReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class HeatFluxReal in libaster


class HeatFluxReal:
    pass

    # Method resolution order:
    #     HeatFluxReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class HydraulicFluxReal in libaster


class HydraulicFluxReal:
    pass

    # Method resolution order:
    #     HydraulicFluxReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def debugPrint(self):
        pass

    def setValue(self, arg0, arg1):
        pass


# class ThermalLoadReal in libaster


class ThermalLoadReal(DataStructure):
    pass

    # Method resolution order:
    #     ThermalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ThermalLoadReal, arg0: Model) -> None

        2. __init__(self: libaster.ThermalLoadReal, arg0: str, arg1: Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getMesh(self):
        pass

    def getModel(self):
        pass

    def getThermalLoadDescription(self):
        pass

    def hasLoadField(self, arg0):
        """Return true if the wanted field exists

        Arguments:
            str: name of the load field

        Returns:
            bool: field exists
        """

    def hasLoadResult(self):
        """Return true if the LoadResult structure exists

        Returns:
            bool: field exists
        """


# class ThermalLoadFunction in libaster


class ThermalLoadFunction(DataStructure):
    pass

    # Method resolution order:
    #     ThermalLoadFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ThermalLoadFunction, arg0: Model) -> None

        2. __init__(self: libaster.ThermalLoadFunction, arg0: str, arg1: Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getMesh(self):
        pass

    def getModel(self):
        pass

    def hasLoadField(self, arg0):
        """Return true if the wanted field exists

        Arguments:
            str: name of the load field

        Returns:
            bool: field exists
        """

    def hasLoadResult(self):
        """Return true if the LoadResult structure exists

        Returns:
            bool: field exists
        """


# class ThermalLoadDescriptionReal in libaster


class ThermalLoadDescriptionReal(DataStructure):
    pass

    # Method resolution order:
    #     ThermalLoadDescriptionReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0, arg1):
        pass

    def getConstantLoadField(self, arg0):
        pass


# class BehaviourDefinition in libaster


class BehaviourDefinition(DataStructure):
    pass

    # Method resolution order:
    #     BehaviourDefinition
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.BehaviourDefinition) -> None

        2. __init__(self: libaster.BehaviourDefinition, arg0: str) -> None
        """


# class Material in libaster


class Material(DataStructure):
    pass

    # Method resolution order:
    #     Material
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Material) -> None

        2. __init__(self: libaster.Material, arg0: str) -> None

        3. __init__(self: libaster.Material, arg0: libaster.Material) -> None

        4. __init__(self: libaster.Material, arg0: libaster.Material, arg1: list[str]) -> None
        """

    def getFunction(self, materialName, propertyName):
        """Return the value of a property stored as a function.

        Raise an exception if the property does not exist.

        Arguments:
            materialName (str): Material name (without "_FO").
            propertyName (str): Property name.

        Returns:
            *Function*: Function object, *None* if the property does not exist or is not a function.
        """

    def getMaterialNames(self):
        """Return the list of the material names.

        Returns:
            list[str]: List of material names (without "_FO")
        """

    def getValueComplex(self, materialName, propertyName):
        """Return the value of a property stored as a complex.

        Raise an exception if the property does not exist.

        Arguments:
            materialName (str): Material name (without "_FO").
            propertyName (str): Property name.

        Returns:
            complex: Property value.
        """

    def getValueReal(self, materialName, propertyName):
        """Return the value of a property stored as a real.

        Raise an exception if the property does not exist.

        Arguments:
            materialName (str): Material name (without "_FO").
            propertyName (str): Property name.

        Returns:
            float: Property value.
        """

    def size(self):
        """Return the number of material names.

        Returns:
            int: Number of material names.
        """


# class PartOfMaterialField in libaster


class PartOfMaterialField:
    pass

    # Method resolution order:
    #     PartOfMaterialField
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.PartOfMaterialField) -> None

        2. __init__(self: libaster.PartOfMaterialField, arg0: list[libaster.Material], arg1: libaster.MeshEntity) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getMeshEntity(self):
        pass

    def getVectorOfMaterial(self):
        pass


# class MaterialField in libaster


class MaterialField(DataStructure):
    pass

    # Method resolution order:
    #     MaterialField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MaterialField, arg0: libaster.BaseMesh) -> None

        2. __init__(self: libaster.MaterialField, arg0: str, arg1: libaster.BaseMesh) -> None
        """

    def addBehaviourOnGroupOfCells(self, behaviour, nameOfGroups):
        """Add behaviour (from DEFI_COMPOR) on group of cells

        Arguments:
            behaviour (BehaviourDefinition): Behaviour (from DEFI_COMPOR)
            nameOfGroups (list(str)) : list of names of groups of cells
        """

    def addBehaviourOnMesh(self, behaviour):
        """Add behaviour (from DEFI_COMPOR) on mesh

        Arguments:
            behaviour (BehaviourDefinition): Behaviour (from DEFI_COMPOR)
        """

    def addExternalStateVariable(self, exteVari):
        """Add external state variable in material field

        Arguments:
            exteVari (ExternalStateVariablePtr): external state variable
        """

    def addMaterialOnGroupOfCells(self, material, nameOfGroups):
        """Add a material properties on list of groups of cells

        Arguments:
            material (Material): material properties
            nameOfGroups (list(str)) : list of names of groups of cells
        """

    def addMaterialOnMesh(self, material):
        """Add material properties on mesh

        Arguments:
            material (Material): material properties
        """

    def addMultipleMaterialOnGroupOfCells(self, material, nameOfGroups):
        """Add a vector of multiple material properties on group of cells

        Arguments:
            material (list(Material)): list of material properties
            nameOfGroups (list(str)) : list of names of groups of cells
        """

    def addMultipleMaterialOnMesh(self, material):
        """Add a vector of multiple material properties on mesh

        Arguments:
            material (list(Material)): list of material properties
        """

    def build(self):
        """Build material field"""

    def getExtStateVariablesOnMeshEntities(self):
        pass

    def getMaterialOnCell(self, arg0):
        """Get the material properties on a giver cell on the material field

        Returns:
            Material: material properties
        """

    def getMaterialsOnMeshEntities(self):
        pass

    def getMesh(self):
        """Get mesh of material field

        Returns:
            BaseMesh: mesh
        """

    def getVectorOfMaterial(self):
        """Get vector of all the material properties on the material field

        Returns:
            list(Material): vector of material properties
        """

    def getVectorOfPartOfMaterialField(self):
        """Get vector of all the material properties with mesh entities on the material field

        Returns:
            list(PartOfMaterial): vector of material properties with mesh entities
        """

    def hasExternalStateVariable(self, *args, **kwargs):
        """Overloaded function.

        1. hasExternalStateVariable(self: libaster.MaterialField, exteVariIden: externVarEnumInt) -> bool


                    Detects the presence of an external state variable

                    Arguments:
                        exteVariIden (externVarEnumInt or str): identifier for external state variable

                    Returns:
                        bool: True if has external state variables


        2. hasExternalStateVariable(self: libaster.MaterialField, exteVariIden: str) -> bool


                    Detects the presence of an external state variable

                    Returns:
                        bool: True if has external state variables


        3. hasExternalStateVariable(self: libaster.MaterialField) -> bool


                    Detects the presence of any external state variable

                    Returns:
                        bool: True if has external state variables
        """

    def hasExternalStateVariableForLoad(self):
        """Detects the presence of an external state variable for loads

        Returns:
            bool: True if has external state variables for loads
        """

    def hasExternalStateVariableWithReference(self):
        """Detects the presence of an external state variable with reference value

        Returns:
            bool: True if has external state variables with reference value
        """

    def setModel(self, model):
        """Set model of the material field

        Arguments:
            model (Model): model
        """

    def updateInternalState(self):
        """Update the internal state of the datastructure.

        Returns:
            bool: *True* in case of success, *False* otherwise.
        """


# class Grid in libaster


class Grid(Mesh):
    pass

    # Method resolution order:
    #     Grid
    #     Mesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Grid) -> None

        2. __init__(self: libaster.Grid, arg0: str) -> None
        """


# class MeshesMapping in libaster


class MeshesMapping(DataStructure):
    pass

    # Method resolution order:
    #     MeshesMapping
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MeshesMapping) -> None

        2. __init__(self: libaster.MeshesMapping, arg0: str) -> None
        """

    def getCoefficients(self):
        """Return the coefficients of the interpolation of the slave nodes on the master cells

        Returns:
            list[real] : interpolation coefficients for each slave node
        """

    def getFirstMesh(self):
        pass

    def getNodesIds(self):
        """Return the ids of the master nodes for the interpolation of the slave nodes

        Returns:
            list[int] : master nodes ids for each slave node
        """

    def getNumberOfMasterNodes(self):
        """Return the number of master nodes implied in the interpolation of the slave nodes

        Returns:
            list[int] : number of master nodes for each slave node
        """

    def getSecondMesh(self):
        pass

    def setFirstMesh(self, arg0):
        pass

    def setSecondMesh(self, arg0):
        pass


# class Skeleton in libaster


class Skeleton(BaseMesh):
    pass

    # Method resolution order:
    #     Skeleton
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Skeleton) -> None

        2. __init__(self: libaster.Skeleton, arg0: str) -> None
        """


# class DynamicMacroElement in libaster


class DynamicMacroElement(DataStructure):
    pass

    # Method resolution order:
    #     DynamicMacroElement
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.DynamicMacroElement) -> None

        2. __init__(self: libaster.DynamicMacroElement, arg0: str) -> None
        """

    def getDOFNumbering(self):
        pass

    def getDampingMatrix(self):
        pass

    def getGeneralizedDampingMatrix(self):
        pass

    def getGeneralizedMassMatrix(self):
        pass

    def getGeneralizedStiffnessMatrix(self):
        pass

    def getImpedanceDampingMatrix(self):
        pass

    def getImpedanceMassMatrix(self):
        pass

    def getImpedanceMatrix(self):
        pass

    def getImpedanceStiffnessMatrix(self):
        pass

    def getMassMatrix(self):
        pass

    def getMechanicalMode(self):
        pass

    def getNumberOfNodes(self):
        pass

    def getStiffnessMatrixComplex(self):
        pass

    def getStiffnessMatrixReal(self):
        pass

    def setDampingMatrix(self, arg0):
        pass

    def setImpedanceDampingMatrix(self, arg0):
        pass

    def setImpedanceMassMatrix(self, arg0):
        pass

    def setImpedanceMatrix(self, arg0):
        pass

    def setImpedanceStiffnessMatrix(self, arg0):
        pass

    def setMassMatrix(self, arg0):
        pass

    def setMechanicalMode(self, arg0):
        pass

    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.

        1. setStiffnessMatrix(self: libaster.DynamicMacroElement, arg0: libaster.AssemblyMatrixDisplacementComplex) -> bool

        2. setStiffnessMatrix(self: libaster.DynamicMacroElement, arg0: libaster.AssemblyMatrixDisplacementReal) -> bool
        """


# class StaticMacroElement in libaster


class StaticMacroElement(DataStructure):
    pass

    # Method resolution order:
    #     StaticMacroElement
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.StaticMacroElement) -> None

        2. __init__(self: libaster.StaticMacroElement, arg0: str) -> None
        """


# class SuperMesh in libaster


class SuperMesh(Mesh):
    pass

    # Method resolution order:
    #     SuperMesh
    #     Mesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.SuperMesh) -> None

        2. __init__(self: libaster.SuperMesh, arg0: str) -> None
        """

    def addDynamicMacroElement(self, arg0):
        """Add a dynamic macro element."""

    def addStaticMacroElement(self, arg0):
        """Add a static macro element."""

    def build(self):
        """Returns:
        bool: true if building is ok
        """

    def getDynamicMacroElements(self):
        """Return all dynamic macro elements."""

    def getNodeLabels(self):
        """Get node labels."""

    def getStaticMacroElements(self):
        """Return all static macro elements."""


# class CrackShape in libaster


class CrackShape:
    pass

    # Method resolution order:
    #     CrackShape
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self):
        pass

    def __setstate__(self, arg0):
        pass

    def getCenter(self):
        pass

    def getCrackSide(self):
        pass

    def getEndPoint(self):
        pass

    def getFilletRadius(self):
        pass

    def getHalfLength(self):
        pass

    def getNormal(self):
        pass

    def getSemiMajorAxis(self):
        pass

    def getSemiMinorAxis(self):
        pass

    def getShape(self):
        pass

    def getShapeName(self):
        pass

    def getStartingPoint(self):
        pass

    def getTangent(self):
        pass

    def getVectX(self):
        pass

    def getVectY(self):
        pass

    def setCylinderCrackShape(self, arg0, arg1, arg2, arg3, arg4):
        pass

    def setEllipseCrackShape(self, arg0, arg1, arg2, arg3, arg4, arg5):
        pass

    def setHalfLineCrackShape(self, arg0, arg1):
        pass

    def setHalfPlaneCrackShape(self, arg0, arg1, arg2):
        pass

    def setLineCrackShape(self, arg0, arg1):
        pass

    def setNotchCrackShape(self, arg0, arg1, arg2, arg3, arg4):
        pass

    def setSegmentCrackShape(self, arg0, arg1):
        pass

    def setSquareCrackShape(self, arg0, arg1, arg2, arg3, arg4, arg5, arg6):
        pass


# class Crack in libaster


class Crack(DataStructure):
    pass

    # Method resolution order:
    #     Crack
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Crack) -> None

        2. __init__(self: libaster.Crack, arg0: str) -> None
        """

    def getConfigInit(self):
        """Return the crack initial configuration"""

    def getCrackFrontAbsCurv(self):
        """Return the crack front absc curv

        Returns:
            list[float]: the crack front absc curv
        """

    def getCrackFrontBasis(self):
        """Return the crack front basis

        Returns:
            list[float]: the crack front basis
        """

    def getCrackFrontNodeBasis(self):
        """Return the basis at each crack front node

        Returns:
            FieldOnNodesReal: field of the crack front basis
        """

    def getCrackFrontNodes(self):
        """Return the crack front nodes

        Returns:
            list[str]: the crack nodes
        """

    def getCrackFrontPosition(self):
        """Return the crack front Position

        Returns:
            list[float]: the crack front Position
        """

    def getCrackFrontRadius(self):
        """Return the crack front Radius

        Returns:
            float: the crack front Radius
        """

    def getCrackTipCellsType(self):
        """Return the crack front cell type

        Returns:
            str: the crack front cell type
        """

    def getLowerLipGroupName(self):
        """Return the group name used to define lower side of cracklip

        Returns:
            str: group name
        """

    def getLowerNormNodes(self):
        """Return the names for nodes on the lower side of cracklip

        Returns:
            list[str]: node names
        """

    def getLowerNormNodes2(self):
        """Return the names for nodes on the lower side of cracklip (for POST_JMOD)

        Returns:
            list[str]: node names
        """

    def getNormal(self):
        """Return vector normal of the crack surface

        Returns:
            list[float]: normal to the crack surface
        """

    def getUpperLipGroupName(self):
        """Return the group name used to define upper side of cracklip

        Returns:
            str: group name
        """

    def getUpperNormNodes(self):
        """Return the names for nodes on the upper side of cracklip

        Returns:
            list[str]: node names
        """

    def getUpperNormNodes2(self):
        """Return the names for nodes on the upper side of cracklip (for POST_JMOD)

        Returns:
            list[str]: node names
        """

    def isSymmetric(self):
        """Return true if crack is symeric"""


# class GeneralizedModel in libaster


class GeneralizedModel(DataStructure):
    pass

    # Method resolution order:
    #     GeneralizedModel
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedModel) -> None

        2. __init__(self: libaster.GeneralizedModel, arg0: str) -> None
        """

    def addDynamicMacroElement(self, arg0, arg1):
        pass

    def getDynamicMacroElementFromName(self, arg0):
        pass

    def getDynamicMacroElementNames(self):
        pass

    def getDynamicStructureLinks(self):
        pass


# class Physics in libaster


class Physics:
    """Enumeration physics."""

    # Method resolution order:
    #     Physics
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Acoustic = 2

    Mechanics = 0

    Thermal = 1


# class Modelings in libaster


class Modelings:
    """Enumeration of modelings."""

    # Method resolution order:
    #     Modelings
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    AXIS_FLUIDE = 67

    AXIS_FLUI_ABSO = 68

    AXIS_FLUI_STRU = 69

    AXIS_FOURIER = 70

    AXIS_GRAD_INCO = 71

    AXIS_GRAD_VARI = 73

    AXIS_HH2D = 74

    AXIS_HH2MD = 75

    AXIS_HH2MS = 76

    AXIS_HH2S = 77

    AXIS_HHD = 78

    AXIS_HHM = 79

    AXIS_HHMD = 80

    AXIS_HHMS = 81

    AXIS_HHS = 83

    AXIS_HM = 84

    AXIS_HMD = 85

    AXIS_HMS = 86

    AXIS_INCO_UP = 87

    AXIS_INCO_UPG = 88

    AXIS_INCO_UPO = 89

    AXIS_INTERFACE = 90

    AXIS_INTERFACE_S = 91

    AXIS_JHMS = 92

    AXIS_JOINT = 93

    AXIS_SECH = 94

    AXIS_SECH_DIAG = 95

    AXIS_SI = 96

    AXIS_THH2D = 97

    AXIS_THH2MD = 98

    AXIS_THH2MS = 99

    AXIS_THH2S = 100

    AXIS_THHD = 101

    AXIS_THHMD = 102

    AXIS_THHMS = 103

    AXIS_THHS = 104

    AXIS_THM = 105

    AXIS_THMD = 106

    AXIS_THMS = 107

    AXIS_THVD = 108

    AXIS_THVS = 109

    Axisymmetrical = 65

    BARRE = 110

    BARRE_2D_NL = 1

    BARRE_3D_NL = 111

    CABLE = 112

    CABLE_GAINE = 113

    CABLE_POULIE = 114

    COQUE_3D = 116

    COQUE_AXIS = 117

    COQUE_SOLIDE = 119

    C_PLAN_SI = 121

    DIL_3D = 11

    DIS_T = 122

    DIS_TR = 123

    DIS_TR_2D = 3

    DIS_T_2D = 2

    DKT = 124

    DKTG = 125

    DST = 126

    D_PLAN_2DG = 128

    D_PLAN_ABSO = 129

    D_PLAN_DIL = 130

    D_PLAN_GRAD_HHO = 131

    D_PLAN_GRAD_INCO = 132

    D_PLAN_GRAD_SIGM = 133

    D_PLAN_GRAD_VARI = 134

    D_PLAN_HH2D = 135

    D_PLAN_HH2MD = 136

    D_PLAN_HH2MS = 137

    D_PLAN_HH2MS_DIL = 138

    D_PLAN_HH2M_SI = 139

    D_PLAN_HH2S = 140

    D_PLAN_HH2SUDA = 141

    D_PLAN_HHD = 142

    D_PLAN_HHM = 143

    D_PLAN_HHMD = 144

    D_PLAN_HHMS = 145

    D_PLAN_HHO = 146

    D_PLAN_HHS = 147

    D_PLAN_HM = 148

    D_PLAN_HMD = 149

    D_PLAN_HMS = 150

    D_PLAN_HMS_DIL = 151

    D_PLAN_HM_SI = 152

    D_PLAN_HM_SI_DIL = 153

    D_PLAN_HS = 154

    D_PLAN_INCO_UP = 155

    D_PLAN_INCO_UPG = 156

    D_PLAN_INCO_UPO = 157

    D_PLAN_MIX_STA = 158

    D_PLAN_SI = 159

    D_PLAN_THH2D = 160

    D_PLAN_THH2MD = 161

    D_PLAN_THH2MS = 162

    D_PLAN_THH2S = 163

    D_PLAN_THHD = 164

    D_PLAN_THHMD = 165

    D_PLAN_THHMS = 166

    D_PLAN_THHS = 167

    D_PLAN_THM = 168

    D_PLAN_THMD = 169

    D_PLAN_THMS = 170

    D_PLAN_THMS_DIL = 171

    D_PLAN_THVD = 172

    D_PLAN_THVS = 173

    FAISCEAU_3D = 12

    FLUIDE_2D = 4

    FLUIDE_3D = 13

    FLUI_ABSO_2D = 5

    FLUI_ABSO_3D = 14

    FLUI_PESA_2D = 6

    FLUI_STRU = 175

    FLUI_STRU_2D = 7

    GRAD_HHO_3D = 15

    GRAD_INCO_3D = 16

    GRAD_VARI_3D = 17

    GRILLE_EXCENTRE = 176

    GRILLE_MEMBRANE = 177

    HH2D_3D = 18

    HH2MD_3D = 19

    HH2MS_3D = 20

    HH2MS_DIL_3D = 21

    HH2M_SI_3D = 22

    HH2SUDA_3D = 24

    HH2S_3D = 23

    HHD_3D = 25

    HHMD_3D = 27

    HHMS_3D = 28

    HHM_3D = 26

    HHO_3D = 29

    HHS_3D = 30

    HMD_3D = 32

    HMS_3D = 33

    HMS_DIL_3D = 34

    HM_3D = 31

    HM_SI_3D = 35

    HM_SI_DIL_3D = 36

    HS_3D = 37

    INCO_UPG_3D = 39

    INCO_UPO_3D = 40

    INCO_UP_3D = 38

    INTERFACE_3D = 41

    INTERFACE_S_3D = 42

    JOINT_3D = 44

    JOINT_HYME_3D = 45

    MEMBRANE = 178

    MIX_STA_3D = 48

    PLAN_INTERFACE = 183

    PLAN_INTERFACE_S = 184

    PLAN_JHMS = 185

    PLAN_JOINT = 186

    PLAN_JOINT_HYME = 187

    PLAQ_MITC = 174

    POU_D_E = 188

    POU_D_EM = 189

    POU_D_SQUE = 190

    POU_D_T = 191

    POU_D_TG = 192

    POU_D_TGM = 193

    POU_D_T_GD = 194

    POU_FLUI_STRU = 195

    Planar = 179

    PlanarBar = 0

    PlaneStrain = 127

    PlaneStress = 120

    Q4G = 196

    Q4GG = 197

    SECH_3D = 46

    SECH_3D_DIAG = 47

    SI_3D = 49

    THH2D_3D = 50

    THH2MD_3D = 51

    THH2MS_3D = 52

    THH2S_3D = 53

    THHD_3D = 54

    THHMD_3D = 56

    THHMS_3D = 57

    THHM_3D = 55

    THHS_3D = 58

    THMD_3D = 60

    THMS_3D = 61

    THMS_DIL_3D = 62

    THM_3D = 59

    THVD_3D = 63

    THVS_3D = 64

    TUYAU_3M = 198

    TUYAU_6M = 199

    Tridimensional = 8

    TridimensionalAbsorbingBoundary = 9


# class Formulation in libaster


class Formulation:
    """Enumeration of formulation."""

    # Method resolution order:
    #     Formulation
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Dil = 6

    DilInco = 7

    Linear = 1

    NoFormulation = 0

    Quadratic = 2

    Sta = 8

    StaInco = 9

    UP = 4

    UPPhi = 3

    UPsi = 5


# class ModelSplitingMethod in libaster


class ModelSplitingMethod:
    """Enumeration for model split method ."""

    # Method resolution order:
    #     ModelSplitingMethod
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Centralized = 0

    GroupOfCells = 2

    SubDomain = 1


# class GraphPartitioner in libaster


class GraphPartitioner:
    """Enumeration for graph partitionner."""

    # Method resolution order:
    #     GraphPartitioner
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    Metis = 1

    Scotch = 0


# class Model in libaster


class Model(DataStructure):
    pass

    # Method resolution order:
    #     Model
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Model, arg0: ConnectionMesh) -> None

        2. __init__(self: libaster.Model, arg0: str, arg1: ConnectionMesh) -> None

        3. __init__(self: libaster.Model, arg0: libaster.BaseMesh) -> None

        4. __init__(self: libaster.Model, arg0: libaster.BaseMesh, arg1: bool) -> None

        5. __init__(self: libaster.Model, arg0: str, arg1: libaster.FiniteElementDescriptor) -> None

        6. __init__(self: libaster.Model, arg0: str, arg1: libaster.FiniteElementDescriptor, arg2: bool) -> None
        """

    def addModelingOnGroupOfCells(self, physics, modeling, grpma, formulation=0):
        """Add modeling on all mesh

        Arguments:
            physics (Physics): Physics
            modeling (Modelings): Modeling
            grpma (str): Name of element group
            formulation (Formulation): Formulation (optional)
        """

    def addModelingOnMesh(self, physics, modeling, formulation=0):
        """Add modeling on all mesh

        Arguments:
            physics (Physics): Physics
            modeling (Modelings): Modeling
            formulation (Formulation): Formulation (optional)
        """

    def banBalancing(self):
        """Prohibit model balancing"""

    def build(self):
        pass

    def existsHHO(self):
        pass

    def existsMultiFiberBeam(self):
        pass

    def existsPartition(self):
        pass

    def existsRdM(self):
        """To know if the model has RdM elements

        Returns:
            bool: *True* if the model contains beam, shell or discret elements, else *False*
        """

    def existsThm(self):
        pass

    def getConnectionMesh(self):
        """Return the ConnectionMesh

        Returns:
            ConnectionMesh: a pointer to the ConnectionMesh
        """

    def getFiniteElementDescriptor(self):
        pass

    def getGeometricDimension(self):
        """To know the geometric dimension supported by the model

        Returns:
            int: geometric dimension
        """

    def getGraphPartitioner(self):
        pass

    def getMesh(self):
        """Return the mesh

        Returns:
            Mesh: a pointer to the mesh
        """

    def getModelisationName(self):
        """Get modelisation name used in model

        Returns:
            str: modelisation name if single modelisation, else '#PLUSIEURS'
        """

    def getPartitionMethod(self):
        """Get partition method

        Returns:
            str: partition method
        """

    def getPhysics(self):
        """To know the physics supported by the model

        Returns:
            str: Mechanics or Thermal or Acoustic
        """

    def getSaneModel(self):
        pass

    def getSplittingMethod(self):
        pass

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """

    def getXfemContact(self):
        pass

    def isAcoustic(self):
        """To know if the model is acoustic or not

        Returns:
            bool: True - if the model is acoustic
        """

    def isAxis(self):
        """To know if the model is Axisymmetric

        Returns:
            bool: *True* if the model is axisymmetric, else *False*
        """

    def isMechanical(self):
        """To know if the model is mechanical or not

        Returns:
            bool: True - if the model is mechanical
        """

    def isThermal(self):
        """To know if the model is thermal or not

        Returns:
            bool: True - if the model is thermal
        """

    def isXfem(self):
        pass

    def setFrom(self, model):
        """Set a model defined on a ConnectionMesh from an other model

        Arguments:
            model (Model): Table identifier.
        """

    def setSaneModel(self, arg0):
        pass

    def setSplittingMethod(self, *args, **kwargs):
        """Overloaded function.

        1. setSplittingMethod(self: libaster.Model, arg0: libaster.ModelSplitingMethod, arg1: libaster.GraphPartitioner) -> None

        2. setSplittingMethod(self: libaster.Model, arg0: libaster.ModelSplitingMethod) -> None
        """

    def xfemPreconditioningEnable(self):
        pass


# class PrestressingCable in libaster


class PrestressingCable(DataStructure):
    pass

    # Method resolution order:
    #     PrestressingCable
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.PrestressingCable, arg0: libaster.Model, arg1: libaster.MaterialField, arg2: libaster.ElementaryCharacteristics) -> None

        2. __init__(self: libaster.PrestressingCable, arg0: str, arg1: libaster.Model, arg2: libaster.MaterialField, arg3: libaster.ElementaryCharacteristics) -> None
        """

    def getElementaryCharacteristics(self):
        pass

    def getMaterialField(self):
        pass

    def getModel(self):
        """Return the Model.

        Returns:
            *Model*: Model object.
        """


# class XfemCrack in libaster


class XfemCrack(DataStructure):
    pass

    # Method resolution order:
    #     XfemCrack
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.XfemCrack, arg0: libaster.Mesh) -> None

        2. __init__(self: libaster.XfemCrack, arg0: str, arg1: libaster.Mesh) -> None
        """

    def build(self):
        pass

    def enrichModelWithXfem(self, arg0):
        pass

    def getAuxiliaryGrid(self):
        pass

    def getCohesiveCrackTipForPropagation(self):
        pass

    def getCrackFrontRadius(self):
        pass

    def getCrackLipsEntity(self):
        pass

    def getCrackShape(self):
        pass

    def getCrackTipBasis(self):
        pass

    def getCrackTipCoords(self):
        pass

    def getCrackTipEntity(self):
        pass

    def getCrackTipMultiplicity(self):
        pass

    def getCrackTipNodeFacesField(self):
        pass

    def getDiscontinuityType(self):
        pass

    def getDiscontinuousField(self):
        pass

    def getEnrichedCells(self):
        pass

    def getEnrichedLayersNumber(self):
        pass

    def getEnrichmentRadiusZone(self):
        pass

    def getEnrichmentType(self):
        pass

    def getExistingCrackWithGrid(self):
        pass

    def getJunctingCracks(self):
        pass

    def getMesh(self):
        pass

    def getNormalLevelSetField(self):
        pass

    def getNormalLevelSetFunction(self):
        pass

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """

    def getTangentialLevelSetField(self):
        pass

    def getTangentialLevelSetFunction(self):
        pass

    def getTipType(self):
        pass

    def insertJunctingCracks(self, arg0):
        pass

    def setAuxiliaryGrid(self, arg0):
        pass

    def setCohesiveCrackTipForPropagation(self, arg0):
        pass

    def setCrackLipsEntity(self, arg0):
        pass

    def setCrackShape(self, arg0):
        pass

    def setCrackTipEntity(self, arg0):
        pass

    def setDiscontinuityType(self, arg0):
        pass

    def setDiscontinuousField(self, arg0):
        pass

    def setEnrichedCells(self, arg0):
        pass

    def setEnrichedLayersNumber(self, arg0):
        pass

    def setEnrichmentRadiusZone(self, arg0):
        pass

    def setEnrichmentType(self, arg0):
        pass

    def setExistingCrackWithGrid(self, arg0):
        pass

    def setMesh(self, arg0):
        pass

    def setNormalLevelSetField(self, arg0):
        pass

    def setNormalLevelSetFunction(self, arg0):
        pass

    def setPointForJunction(self, arg0):
        pass

    def setTangentialLevelSetField(self, arg0):
        pass

    def setTangentialLevelSetFunction(self, arg0):
        pass

    def updateInternalState(self):
        pass


# class Result in libaster


class Result(DataStructure):
    pass

    # Method resolution order:
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.Result, arg0: str) -> None

        2. __init__(self: libaster.Result, arg0: str, arg1: str) -> None
        """

    def addEquationNumbering(self, arg0):
        pass

    def addFiniteElementDescriptor(self, arg0):
        pass

    def allocate(self, nb_index):
        """Allocate result

        Arguments:
            nb_index (int):  number of index to allocate
        """

    def clear(self, *args, **kwargs):
        """Overloaded function.

        1. clear(self: libaster.Result) -> None


        Clear fields, models, parameters, ... in result


        2. clear(self: libaster.Result, index: int) -> None


        Clear fields, models, parameters, ... in result from the given index

        Arguments:
            index (int): index from begin cleaning
        """

    def createIndexFromParameter(self, para_name, para_value):
        """Create an index in the result

        Arguments:
            para_name (str): parameter name to store
            para_value (str): parameter value to store
        """

    def exists(self):
        """The result exists or nor

        Returns:
            bool: True if the result exists else False
        """

    def getAccessParameters(self, only_access=True):
        """Return the access parameters of the result as Python dict.

        Returns:
            dict{str : list[int,float,str]}: Dict of values for each access variable.
        """

    def getAllElementaryCharacteristics(self):
        """Return the list of all elementary characteristics used in the result

        Returns:
            list[ElementaryCharacteristics]: list of ElementaryCharacteristics.
        """

    def getConstantFieldsOnCellsChar16Names(self):
        """Return the names of the contant char16 fields on cells as Python list.

        Returns:
            list(str): List of names of the contant fields on cells.
        """

    def getConstantFieldsOnCellsRealNames(self):
        """Return the names of the contant real fields on cells as Python list.

        Returns:
            list(str): List of names of the contant fields on cells.
        """

    def getElementaryCharacteristics(self, *args, **kwargs):
        """Overloaded function.

        1. getElementaryCharacteristics(self: libaster.Result, index: int) -> libaster.ElementaryCharacteristics


        Get elementary characterictics at the specfied index

        Arguments:
            index (int): index

        Returns:
            ElementaryCharacteristics: a pointer to elementary characterictics.


        2. getElementaryCharacteristics(self: libaster.Result) -> libaster.ElementaryCharacteristics
        """

    def getEquationNumberings(self):
        """Get list of field's description to build internal FieldOnNodes

        Returns:
            list[EquationNumbering]: list of field's description
        """

    def getFieldsNames(self, *args, **kwargs):
        """Overloaded function.

        1. getFieldsNames(self: libaster.Result) -> list[str]


        Return the list of names of stored fields

        Returns:
            list[str]: List of names of stored fields.


        2. getFieldsNames(self: libaster.Result) -> list[str]


        Return the list of names of stored fields

        Returns:
            list[str]: List of names of stored fields.
        """

    def getFieldsOnCellsComplexNames(self):
        """Return the names of the complex fields on cells as Python list.

        Returns:
            list(str): List of names of the complex fields on cells.
        """

    def getFieldsOnCellsLongNames(self):
        """Return the names of the integer fields on cells as Python list.

        Returns:
            list(str): List of names of the integer fields on cells.
        """

    def getFieldsOnCellsRealNames(self):
        """Return the names of the real fields on cells as Python list.

        Returns:
            list(str): List of names of the real fields on cells.
        """

    def getFieldsOnNodesComplexNames(self):
        """Return the names of the complex fields on nodes as Python list.

        Returns:
            list(str): List of names of the complex fields on nodes.
        """

    def getFieldsOnNodesRealNames(self):
        """Return the names of the real fields on nodes as Python list.

        Returns:
            list(str): List of names of the real fields on nodes.
        """

    def getFiniteElementDescriptors(self):
        """Get list of finite element descriptor to build internal FieldOnCells

        Returns:
            list[FiniteElementDescriptor]: list of finite element descriptor
        """

    def getFirstIndex(self):
        """Get the first index stored in the result

        Returns:
            int: first index stored.
        """

    def getGeneralizedVectorComplexNames(self):
        """Return the names of the complex generalized vectors as Python list.

        Returns:
            list(str): List of names of the complex generalized vectors.
        """

    def getGeneralizedVectorRealNames(self):
        """Return the names of the real generalized vectors as Python list.

        Returns:
            list(str): List of names of the real generalized vectors.
        """

    def getIndexes(self):
        """Return the list of indexes used to store fields

        Returns:
            list[int]: List of indexs used to store fields.
        """

    def getIndexesForFieldName(self, arg0):
        """Returns the list of indexes used to store a specific field,
         indicated by its name.

        Returns:
            list[int]: List of indexs used to store fields.
        """

    def getLastIndex(self):
        """Get the last index stored in the result

        Returns:
            int: last index stored.
        """

    def getLastTime(self):
        """Get the last time value stored in the result

        Returns:
            float: last time value.
        """

    def getListOfLoads(self, index):
        """Get list of loads on the specified index

        Arguments:
            index (int): index to get

        Returns:
            ListOfLoads: a pointer to list of loads.
        """

    def getMaterialField(self, *args, **kwargs):
        """Overloaded function.

        1. getMaterialField(self: libaster.Result, index: int) -> libaster.MaterialField


        Return the material field for the given index.

        Arguments:
            index (int): index

        Returns:
            MaterialField: Material field.


        2. getMaterialField(self: libaster.Result) -> libaster.MaterialField
        """

    def getMaterialFields(self):
        """Return the list of all material fields used in the result

        Returns:
            list[MaterialField]: list of material field.
        """

    def getMesh(self):
        """Return a pointer to mesh

        Returns:
            mesh (Mesh): a pointer to the mesh.
        """

    def getModel(self, *args, **kwargs):
        """Overloaded function.

        1. getModel(self: libaster.Result, index: int) -> libaster.Model


        Return the model for the given index.

        Arguments:
            index (int): index

        Returns:
            Model: Model object.


        2. getModel(self: libaster.Result) -> libaster.Model
        """

    def getModels(self):
        """Return the list of all models used in the result

        Returns:
            list[Model]: list of models.
        """

    def getNumberOfIndexes(self):
        """Get the number of index stored in the result

        Returns:
            int: number of index stored.
        """

    def getParameters(self, only_access=False):
        """Return the parameters of the result as Python dict.

        Returns:
            dict{str : list[int,float,str]}: Dict of values for each parameter variable.
        """

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """

    def getTime(self, index):
        """Get time at the specified index

        Arguments:
            index (int):  index where to save time value

        Returns
            float: time value
        """

    def hasElementaryCharacteristics(self, *args, **kwargs):
        """Overloaded function.

        1. hasElementaryCharacteristics(self: libaster.Result, index: int) -> bool


        Test if a elementary characterictics is used at the specfied index

        Arguments:
            index (int): index

        Returns:
            bool: *True* if at least one elementary characterictics used else *False*.


        2. hasElementaryCharacteristics(self: libaster.Result) -> bool
        """

    def hasListOfLoads(self, *args, **kwargs):
        """Overloaded function.

        1. hasListOfLoads(self: libaster.Result, index: int) -> bool


        Test if a list of loads is used at the specfied index

        Arguments:
            index (int): index

        Returns:
            bool: *True* if at least one list of loads is used else *False*.


        2. hasListOfLoads(self: libaster.Result) -> bool
        """

    def hasMaterialField(self, index):
        """Test if a material field is used at the specfied index

        Arguments:
            index (int): index

        Returns:
            bool: *True* if at a material field used else *False*.
        """

    def hasModel(self, index):
        """Test if a model is used at the specfied index

        Arguments:
            index (int): index

        Returns:
            bool: *True* if at a model used else *False*.
        """

    def printInfo(self):
        pass

    def printListOfFields(self):
        """Print the names of all fields (real, complex, ...) stored in the result."""

    def printMedFile(
        self, filename, medname="", local=True, internalVar=True, fields=[], version=""
    ):
        """Print the result in a MED file.

        Args:
            filename (Path|str): Path to the output file.
            medname (str): Name of the result in the MED file. (default: "")
            local (bool): Print only the local domain if *True*. (default: True)
            fields (list[str]): Name of fields to save. (default: all)
            version (str): Version of MED file.
        """

    def resize(self, nbIndexes):
        """Resize the object.

        Arguments:
            nbIndexes (int): new expected size. Should be greater than the current size,
                otherwise the size is unchanged.
        """

    def setElementaryCharacteristics(self, *args, **kwargs):
        """Overloaded function.

        1. setElementaryCharacteristics(self: libaster.Result, cara_elem: libaster.ElementaryCharacteristics, exists_ok: bool = False) -> None


        Set elementary characterictics on all indexs

        Arguments:
            cara_elem (ElementaryCharacteristics): elementary characterictics to set.
            exists_ok (bool): If *True*, pass silently if a Model is already defined. *False* by default.


        2. setElementaryCharacteristics(self: libaster.Result, cara_elem: libaster.ElementaryCharacteristics, index: int) -> None


        Set elementary characterictics on the specified index

        Arguments:
            cara_elem (ElementaryCharacteristics): elementary characterictics to set.
            index (int): index to set
        """

    def setListOfLoads(self, load, index):
        """Set list of loads on the specified index

        Arguments:
            load (ListOfLoads): list of loads to set.
            index (int): index to set
        """

    def setMaterialField(self, *args, **kwargs):
        """Overloaded function.

        1. setMaterialField(self: libaster.Result, mater: libaster.MaterialField, exists_ok: bool = False) -> None


        Set material field on all indexs

        Arguments:
            mater (MaterialField): material field to set.
            exists_ok (bool): If *True*, pass silently if a Model is already defined. *False* by default.


        2. setMaterialField(self: libaster.Result, mater: libaster.MaterialField, index: int) -> None


        Set material field on the specified index

        Arguments:
            mater (MaterialField): material field to set.
            index (int): index to set
        """

    def setMesh(self, mesh):
        """Set the mesh used by the result.

        Arguments:
            mesh (BaseMesh): mesh to set
        """

    def setModel(self, *args, **kwargs):
        """Overloaded function.

        1. setModel(self: libaster.Result, model: libaster.Model, exists_ok: bool = False) -> None


        Set model on all indexs

        Arguments:
            model (Model): Model to be assigned.
            exists_ok (bool): If *True*, pass silently if a Model is already defined. *False* by default.


        2. setModel(self: libaster.Result, model: libaster.Model, index: int) -> None


        Set model on the specified index

        Arguments:
            model (Model): model to set
            index (int): index to set
        """

    def setParameterValue(self, *args, **kwargs):
        """Overloaded function.

        1. setParameterValue(self: libaster.Result, para_name: str, para_value: float, index: int) -> None


        Add parameter at the specified index

        Arguments:
            para_name (float): parameter name to store
            para_value (float): parameter value to store
            index (int):  index where to save value of parameter


        2. setParameterValue(self: libaster.Result, para_name: str, para_value: str, index: int) -> None


        Add parameter at the specified index

        Arguments:
            para_name (float): parameter name to store
            para_value (str): parameter value to store
            index (int):  index where to save value of parameter
        """

    def setTable(self, identifier, table):
        """Set a Table to the datastructure.

        Arguments:
            identifier (str): Table identifier.
            Table: Table to be stored with the given identifier.
        """

    def setTime(self, time, index):
        """Add time at the specified index

        Arguments:
            time (float): time value to save
            index (int):  index where to save time value
        """


# class TransientResult in libaster


class TransientResult(Result):
    pass

    # Method resolution order:
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.TransientResult) -> None

        2. __init__(self: libaster.TransientResult, arg0: str, arg1: str) -> None
        """


# class LoadResult in libaster


class LoadResult(TransientResult):
    pass

    # Method resolution order:
    #     LoadResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.LoadResult) -> None

        2. __init__(self: libaster.LoadResult, arg0: str) -> None
        """


# class ThermalResult in libaster


class ThermalResult(TransientResult):
    pass

    # Method resolution order:
    #     ThermalResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ThermalResult) -> None

        2. __init__(self: libaster.ThermalResult, arg0: str) -> None
        """


# class DryingResult in libaster


class DryingResult(TransientResult):
    pass

    # Method resolution order:
    #     DryingResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.DryingResult) -> None

        2. __init__(self: libaster.DryingResult, arg0: str) -> None
        """


# class CombinedFourierResult in libaster


class CombinedFourierResult(Result):
    pass

    # Method resolution order:
    #     CombinedFourierResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.CombinedFourierResult) -> None

        2. __init__(self: libaster.CombinedFourierResult, arg0: str) -> None
        """


# class ElasticFourierResult in libaster


class ElasticFourierResult(Result):
    pass

    # Method resolution order:
    #     ElasticFourierResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElasticFourierResult) -> None

        2. __init__(self: libaster.ElasticFourierResult, arg0: str) -> None
        """


# class ThermalFourierResult in libaster


class ThermalFourierResult(Result):
    pass

    # Method resolution order:
    #     ThermalFourierResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ThermalFourierResult) -> None

        2. __init__(self: libaster.ThermalFourierResult, arg0: str) -> None
        """


# class MultipleElasticResult in libaster


class MultipleElasticResult(Result):
    pass

    # Method resolution order:
    #     MultipleElasticResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MultipleElasticResult) -> None

        2. __init__(self: libaster.MultipleElasticResult, arg0: str) -> None
        """


# class NonLinearResult in libaster


class NonLinearResult(TransientResult):
    pass

    # Method resolution order:
    #     NonLinearResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.NonLinearResult) -> None

        2. __init__(self: libaster.NonLinearResult, arg0: str) -> None
        """

    def getTangentMatrix(self):
        pass

    def printMedFile(
        self, filename, medname="", local=False, internalVar=True, fields=[], version=""
    ):
        """Print the result in a MED file.

        Args:
            filename (Path|str): Path to the output file.
            medname (str): Name of the result in the MED file. (default: "")
            local (bool): Print only the local domain if *True*. (default: True)
            internalVarl (bool); Save VARI_ELGA_NOMME or not. (default: True)
            fields (list[str]): Name of fields to save. (default: all)
            version (str): Version of MED file.
        """

    def setContact(self, *args, **kwargs):
        """Overloaded function.

        1. setContact(self: libaster.NonLinearResult, arg0: libaster.Contact) -> None

        2. setContact(self: libaster.NonLinearResult, arg0: libaster.Contact, arg1: int) -> None
        """


# class PhysicalProblem in libaster


class PhysicalProblem:
    pass

    # Method resolution order:
    #     PhysicalProblem
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.PhysicalProblem, arg0: libaster.BaseDOFNumbering) -> None

        2. __init__(self: libaster.PhysicalProblem, arg0: libaster.Model, arg1: libaster.MaterialField) -> None

        3. __init__(self: libaster.PhysicalProblem, arg0: libaster.Model, arg1: libaster.MaterialField, arg2: libaster.ElementaryCharacteristics) -> None
        """

    def __setstate__(self, arg0):
        pass

    def addDirichletBC(self, *args, **kwargs):
        """Overloaded function.

        1. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC) -> None

        2. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC, arg1: libaster.Function, arg2: str) -> None

        3. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC, arg1: libaster.Formula, arg2: str) -> None

        4. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC, arg1: libaster.Function2D, arg2: str) -> None

        5. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC, arg1: str) -> None
        """

    def addLoad(self, *args, **kwargs):
        """Overloaded function.

        1. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal) -> None

        2. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction) -> None

        3. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal, arg1: str) -> None

        4. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal, arg1: libaster.Function, arg2: str) -> None

        5. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal, arg1: libaster.Formula, arg2: str) -> None

        6. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal, arg1: libaster.Function2D, arg2: str) -> None

        7. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction, arg1: str) -> None

        8. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Function, arg2: str) -> None

        9. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Formula, arg2: str) -> None

        10. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Function2D, arg2: str) -> None

        11. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex) -> None

        12. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex, arg1: libaster.Function) -> None

        13. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex, arg1: libaster.Formula) -> None

        14. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex, arg1: libaster.Function2D) -> None

        15. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >) -> None

        16. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: str) -> None

        17. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function, arg2: str) -> None

        18. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Formula, arg2: str) -> None

        19. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function2D, arg2: str) -> None

        20. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: str) -> None

        21. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function, arg2: str) -> None

        22. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Formula, arg2: str) -> None

        23. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function2D, arg2: str) -> None

        24. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >) -> None

        25. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function) -> None

        26. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Formula) -> None

        27. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function2D) -> None

        28. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None

        29. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function) -> None

        30. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Formula) -> None

        31. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function2D) -> None

        32. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal) -> None

        33. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal, arg1: libaster.Function) -> None

        34. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal, arg1: libaster.Formula) -> None

        35. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal, arg1: libaster.Function2D) -> None

        36. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction) -> None

        37. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction, arg1: libaster.Function) -> None

        38. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction, arg1: libaster.Formula) -> None

        39. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction, arg1: libaster.Function2D) -> None

        40. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex) -> None

        41. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex, arg1: libaster.Function) -> None

        42. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex, arg1: libaster.Formula) -> None

        43. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex, arg1: libaster.Function2D) -> None
        """

    def computeBehaviourProperty(self, COMPORTEMENT, SIGM_INIT="NON", INFO=1):
        """Create constant fields on cells for behaviour (COMPOR, CARCRI and MULCOM)

        Arguments:
            COMPORTEMENT (list[dict]): keywords as provided to STAT_NON_LINE/COMPORTEMENT
            SIGM_INIT (str): "OUI" if there is an initial stress field
            INFO (int): level of verbosity, 1 to have description of behaviour or 0 to be quiet
        """

    def computeDOFNumbering(self):
        """Build DOF numbering from the model and loads

        Returns:
            Bool: True if success
        """

    def computeListOfLoads(self, command_name=""):
        """Build the list of loads from the added loads

        Arguments:
            command_name (str): It is possible to add a command name to add more checking (default: "")

        Returns:
            Bool: True if success
        """

    def computeReferenceExternalStateVariables(self):
        """Compute field for external state variables reference value

        Returns:
            FieldOnCells: field for external state variables reference values
        """

    def getBehaviourProperty(self):
        """Return the behaviour properties

        Returns:
            BehaviourProperty: a pointer to the behaviour properties
        """

    def getCodedMaterial(self):
        """Return the coded material

        Returns:
            CodedMaterial: a pointer to the coded material
        """

    def getDOFNumbering(self):
        """Return the DOF numbering

        Returns:
            BaseDOFNumbering: a pointer to the DOF numbering
        """

    def getDirichletBCDOFs(self):
        """Return a vector with DOFs eliminated by Dirichlet boundaries conditions (if it exists)

        Returns:
            tuple(int): a vector with DOFs eliminated by Dirichlet boundaries conditions of
                size = neq + 1,
                tuple(ieq = 0, neq - 1) = 1 then DOF eliminated else 0,
                tuple(neq) = number of DOFs eliminated.
        """

    def getElementaryCharacteristics(self):
        """Return the elementary charateristics

        Returns:
            ElementaryCharacteristics: a pointer to the elementary charateristics
        """

    def getExternalStateVariables(self, time):
        """Get the field for external state variables

        Arguments:
            time [float] : time value to evaluate values

        Returns:
            FieldOnCellsReal : external values
        """

    def getListOfLoads(self):
        """Return list of loads.

        Returns:
            ListOfLoads: a pointer to list of loads
        """

    def getMaterialField(self):
        """Return the material field

        Returns:
            MaterialField: a pointer to the material field
        """

    def getMesh(self):
        """Return the mesh

        Returns:
            Mesh: a pointer to the mesh
        """

    def getModel(self):
        """Return the model

        Returns:
            Model: a pointer to the model
        """

    def getReferenceExternalStateVariables(self):
        """Get the field of reference values for external state variables

        Returns:
            FieldOnCellsReal : field of reference values
        """

    def isAcoustic(self):
        """To know if the probleme is acoustic or not

        Returns:
            bool: True - if the model is acoustic
        """

    def isMechanical(self):
        """To know if the problem is mechanical or not

        Returns:
            bool: True - if the model is mechanical
        """

    def isThermal(self):
        """To know if the problem is thermal or not

        Returns:
            bool: True - if the model is thermal
        """

    def setContactFED(self, virtualCell):
        """Set virtual cells from contact pairing

        Arguments:
            virtualCell (FiniteElementDescriptor)): a pointer to the FED
        """

    def setContactSlaveFED(self, contact):
        """Set virtual cells from contact definition

        Arguments:
            virtualCell (FiniteElementDescriptor)): a pointer to the FED
        """

    def setDOFNumbering(self, dofNum):
        """Set the DOF numbering

        Arguments:
            dofNum (BaseDOFNumbering): a pointer to the DOF numbering
        """

    def setListOfLoads(self, loads):
        """Set list of loads

        Arguments:
            loads (ListOfLoads): a pointer to the list of loads
        """

    def zeroDirichletBCDOFs(self, arg0):
        """Set in-place to zero the DOFs with DirichletBC (aka not assigned by Lagrange multipliers)

        Returns:
            field(FieldOnNodes): the modified field
        """


# class Glossary in libaster


class Glossary:
    pass

    # Method resolution order:
    #     Glossary
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getComponent(self, arg0):
        pass

    def getModeling(self, arg0):
        pass

    def getPhysics(self, arg0):
        pass


# built-in function getGlossary in libaster


def getGlossary():
    pass


# class CyclicSymmetryMode in libaster


class CyclicSymmetryMode(DataStructure):
    pass

    # Method resolution order:
    #     CyclicSymmetryMode
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.CyclicSymmetryMode) -> None

        2. __init__(self: libaster.CyclicSymmetryMode, arg0: str) -> None
        """


# class FullResult in libaster


class FullResult(Result):
    pass

    # Method resolution order:
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FullResult, arg0: str, arg1: str) -> None

        2. __init__(self: libaster.FullResult, arg0: str) -> None
        """

    def getDOFNumbering(self):
        pass

    def setDOFNumbering(self, arg0):
        pass


# class ModeResult in libaster


class ModeResult(FullResult):
    pass

    # Method resolution order:
    #     ModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ModeResult) -> None

        2. __init__(self: libaster.ModeResult, arg0: str) -> None
        """

    def getDOFNumbering(self):
        pass

    def getMassMatrix(self):
        pass

    def getNumberOfDynamicModes(self):
        pass

    def getNumberOfStaticModes(self):
        pass

    def getStiffnessMatrix(self):
        pass

    def setMassMatrix(self, *args, **kwargs):
        """Overloaded function.

        1. setMassMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementReal) -> None

        2. setMassMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixTemperatureReal) -> None

        3. setMassMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementComplex) -> None

        4. setMassMatrix(self: libaster.ModeResult, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> None
        """

    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.

        1. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementReal) -> None

        2. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixTemperatureReal) -> None

        3. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementComplex) -> None

        4. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixPressureReal) -> None

        5. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixPressureReal) -> None

        6. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.GeneralizedAssemblyMatrixReal) -> None
        """

    def setStructureInterface(self, arg0):
        pass


# class ModeResultComplex in libaster


class ModeResultComplex(ModeResult):
    pass

    # Method resolution order:
    #     ModeResultComplex
    #     ModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ModeResultComplex) -> None

        2. __init__(self: libaster.ModeResultComplex, arg0: str) -> None
        """

    def setDampingMatrix(self, arg0):
        pass

    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.

        1. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixDisplacementReal) -> bool

        2. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixDisplacementComplex) -> bool

        3. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixTemperatureReal) -> bool

        4. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixPressureReal) -> bool

        5. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.GeneralizedAssemblyMatrixReal) -> bool

        6. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> bool
        """

    def setStructureInterface(self, arg0):
        pass


# class AcousticModeResult in libaster


class AcousticModeResult(FullResult):
    pass

    # Method resolution order:
    #     AcousticModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.AcousticModeResult) -> None

        2. __init__(self: libaster.AcousticModeResult, arg0: str) -> None
        """

    def setStiffnessMatrix(self, arg0):
        pass


# class BucklingModeResult in libaster


class BucklingModeResult(FullResult):
    pass

    # Method resolution order:
    #     BucklingModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.BucklingModeResult) -> None

        2. __init__(self: libaster.BucklingModeResult, arg0: str) -> None
        """

    def getStiffnessMatrix(self):
        pass

    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.

        1. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixDisplacementReal) -> bool

        2. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixDisplacementComplex) -> bool

        3. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixTemperatureReal) -> bool

        4. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixPressureReal) -> bool

        5. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.GeneralizedAssemblyMatrixReal) -> bool

        6. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> bool
        """


# class GeneralizedResultReal in libaster


class GeneralizedResultReal(DataStructure):
    pass

    # Method resolution order:
    #     GeneralizedResultReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""


# class GeneralizedResultComplex in libaster


class GeneralizedResultComplex(DataStructure):
    pass

    # Method resolution order:
    #     GeneralizedResultComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""


# class TransientGeneralizedResult in libaster


class TransientGeneralizedResult(GeneralizedResultReal):
    pass

    # Method resolution order:
    #     TransientGeneralizedResult
    #     GeneralizedResultReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.TransientGeneralizedResult) -> None

        2. __init__(self: libaster.TransientGeneralizedResult, arg0: str) -> None
        """

    def build(self):
        """Builds C++ arguments associated to attributes stored by blocks of time indices"""

    def getAccelerationValues(self, *args, **kwargs):
        """Overloaded function.

        1. getAccelerationValues(self: libaster.TransientGeneralizedResult) -> list[float]


        Return generalized accelerations values for all time indices.

        Returns:
            list[double]: generalized accelerations values.


        2. getAccelerationValues(self: libaster.TransientGeneralizedResult, idx: int) -> list[float]


        Return generalized accelerations values at a given time index.

        Arguments:
            idx (int): time index

        Returns:
            list[double]: generalized accelerations values.
        """

    def getDOFNumbering(self):
        """Get DOF numbering

        Returns:
            DOFNumbering: DOF numbering
        """

    def getDisplacementValues(self, *args, **kwargs):
        """Overloaded function.

        1. getDisplacementValues(self: libaster.TransientGeneralizedResult) -> list[float]


        Return generalized displacements values for all time indices.

        Returns:
            list[double]: generalized displacements values.


        2. getDisplacementValues(self: libaster.TransientGeneralizedResult, idx: int) -> list[float]


        Return generalized displacements values at a given time index.

        Arguments:
            idx (int): time index

        Returns:
            list[double]: generalized displacements values.
        """

    def getGeneralizedDOFNumbering(self):
        """Get generalized DOF numbering

        Returns:
            GeneralizedDOFNumbering: generalized DOF numbering
        """

    def getIndexes(self):
        """Returns time indices of the transient calculation

        Returns:
            list[int]: time indices
        """

    def getNumberOfModes(self):
        """Returns the number of vectors in the generalized basis

        Returns:
            int: number of vectors in the generalized basis
        """

    def getTimes(self):
        """Returns values of instants of the transient calculation

        Returns:
            list[float]: instants values
        """

    def getVelocityValues(self, *args, **kwargs):
        """Overloaded function.

        1. getVelocityValues(self: libaster.TransientGeneralizedResult) -> list[float]


        Return generalized velocities values for all time indices.

        Returns:
            list[double]: generalized velocities values.


        2. getVelocityValues(self: libaster.TransientGeneralizedResult, idx: int) -> list[float]


        Return generalized velocities values at a given time index.

        Arguments:
            idx (int): time index

        Returns:
            list[double]: generalized velocities values.
        """

    def setAccelerationValues(self, idx, val):
        """Set generalized acceleration values at a given time index.

        Arguments:
            idx (int): time index

            val (list[double]): generalized acceleration values.
        """

    def setDOFNumbering(self, dofn):
        """Set DOF numbering

        Arguments:
            dofn (DOFNumbering): DOF numbering
        """

    def setDisplacementValues(self, idx, val):
        """Set generalized displacement values at a given time index.

        Arguments:
            idx (int): time index

            val (list[double]): generalized displacement values.
        """

    def setGeneralizedDOFNumbering(self, dofg):
        """Set generalized DOF numbering

        Arguments:
            dofg (GeneralizedDOFNumbering): generalized DOF numbering
        """

    def setVelocityValues(self, idx, val):
        """Set generalized velocity values at a given time index.

        Arguments:
            idx (int): time index

            val (list[double]): generalized velocity values.
        """


# class HarmoGeneralizedResult in libaster


class HarmoGeneralizedResult(GeneralizedResultComplex):
    pass

    # Method resolution order:
    #     HarmoGeneralizedResult
    #     GeneralizedResultComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.HarmoGeneralizedResult) -> None

        2. __init__(self: libaster.HarmoGeneralizedResult, arg0: str) -> None
        """

    def getDOFNumbering(self):
        pass

    def getDisplacement(self):
        pass

    def getGeneralizedDOFNumbering(self):
        pass

    def setAcceleration(self, arg0):
        pass

    def setDOFNumbering(self, arg0):
        pass

    def setDisplacement(self, arg0):
        pass

    def setGeneralizedDOFNumbering(self, arg0):
        pass

    def setVelocity(self, arg0):
        pass


# class ElasticResult in libaster


class ElasticResult(Result):
    pass

    # Method resolution order:
    #     ElasticResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ElasticResult) -> None

        2. __init__(self: libaster.ElasticResult, arg0: str) -> None
        """


# class MeshCoordinatesField in libaster


class MeshCoordinatesField(DataStructure):
    pass

    # Method resolution order:
    #     MeshCoordinatesField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __add__(self, *args, **kwargs):
        """Overloaded function.

        1. __add__(self: libaster.MeshCoordinatesField, arg0: libaster.MeshCoordinatesField) -> libaster.MeshCoordinatesField

        2. __add__(self: libaster.MeshCoordinatesField, arg0: libaster.FieldOnNodesReal) -> libaster.MeshCoordinatesField

        3. __add__(self: libaster.FieldOnNodesReal, arg0: libaster.MeshCoordinatesField) -> libaster.MeshCoordinatesField
        """

    def __getitem__(self, node_id):
        """Return the coordinates (x,y,z) at of Node node_id in the vector.

        The value is the same as *getValues()[3*node_id:3*node_id+2]* without creating the entire vector.

        Returns:
            tuple[float]: coordinates (x,y,z).
        """

    def __iadd__(self, arg0):
        pass

    def __imul__(self, arg0):
        pass

    def __init__(self, arg0):
        pass

    def __isub__(self, arg0):
        pass

    def __mul__(self, arg0):
        pass

    def __neg__(self):
        pass

    def __rmul__(self, arg0):
        pass

    def __sub__(self, arg0):
        pass

    def copy(self):
        """Return a copy of MeshCoordinatesField object

        Returns:
            MeshCoordinatesField : MeshCoordinatesField object
        """

    def getNode(self, node_id):
        """Return a node

        Arguments:
            node_id [int] : node id

        Returns:
            Node: Node object.
        """

    def getValues(self):
        """Return a list of values of the coordinates as (x1, y1, z1, x2, y2, z2...)

        Returns:
            list[float]: List of coordinates (size = 3 * number of nodes).
        """

    def setNode(self, node):
        """Set a node

        Arguments:
            node [Node] : node to set.
        """

    def size(self):
        """Return the size of the field

        Returns:
            int : number of values of MeshCoordinatesField object
        """

    def toFieldOnNodes(self, mesh):
        """Convert to FieldOnNodes

        Arguments:
            mesh[Mesh]: the mesh where the coordinates come from

        Returns:
            FieldOnNodesReal: the corresponding field
        """

    def toNumpy(self):
        """Return a numpy array view (no-copy) of values of the coordinates with shape (number of nodes, 3).

        Returns:
            np.ndarray: Array view of coordinates with shape=(number of nodes, 3).
        """

    def updateValuePointers(self):
        """Update values of internal pointer."""


# class FullTransientResult in libaster


class FullTransientResult(FullResult):
    pass

    # Method resolution order:
    #     FullTransientResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FullTransientResult, arg0: str) -> None

        2. __init__(self: libaster.FullTransientResult) -> None
        """


# class FullHarmonicResult in libaster


class FullHarmonicResult(FullResult):
    pass

    # Method resolution order:
    #     FullHarmonicResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FullHarmonicResult, arg0: str) -> None

        2. __init__(self: libaster.FullHarmonicResult) -> None
        """


# class FullHarmonicAcousticResult in libaster


class FullHarmonicAcousticResult(FullResult):
    pass

    # Method resolution order:
    #     FullHarmonicAcousticResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FullHarmonicAcousticResult, arg0: str) -> None

        2. __init__(self: libaster.FullHarmonicAcousticResult) -> None
        """


# class FluidStructureModalBasis in libaster


class FluidStructureModalBasis(DataStructure):
    pass

    # Method resolution order:
    #     FluidStructureModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.FluidStructureModalBasis) -> None

        2. __init__(self: libaster.FluidStructureModalBasis, arg0: str) -> None
        """

    def getTable(self, identifier):
        """Extract a Table from the datastructure.

        Arguments:
            identifier (str): Table identifier.

        Returns:
            Table: Table stored with the given identifier.
        """


# class GeneralizedModeResult in libaster


class GeneralizedModeResult(FullResult):
    pass

    # Method resolution order:
    #     GeneralizedModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.GeneralizedModeResult, arg0: str) -> None

        2. __init__(self: libaster.GeneralizedModeResult) -> None
        """

    def getDampingMatrix(self):
        pass

    def getGeneralizedDOFNumbering(self):
        pass

    def getGeneralizedVectorComplex(self, arg0, arg1):
        pass

    def getGeneralizedVectorReal(self, arg0, arg1):
        pass

    def getStiffnessMatrix(self):
        pass

    def setDampingMatrix(self, arg0):
        pass

    def setGeneralizedDOFNumbering(self, arg0):
        pass

    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.

        1. setStiffnessMatrix(self: libaster.GeneralizedModeResult, arg0: libaster.GeneralizedAssemblyMatrixReal) -> bool

        2. setStiffnessMatrix(self: libaster.GeneralizedModeResult, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> bool
        """


# class MGISBehaviour in libaster


class MGISBehaviour(DataStructure):
    pass

    # Method resolution order:
    #     MGISBehaviour
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MGISBehaviour) -> None

        2. __init__(self: libaster.MGISBehaviour, arg0: str) -> None
        """

    def setBehaviourName(self, name):
        """Define the name of the behaviour to be used from the MFront library.

        Arguments:
            name: Name of the behaviour.
        """

    def setLibPath(self, path):
        """Set the path to the MFront library.

        Arguments:
            path: Library path.
        """


# class ParallelMesh in libaster


class ParallelMesh(BaseMesh):
    pass

    # Method resolution order:
    #     ParallelMesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelMesh) -> None

        2. __init__(self: libaster.ParallelMesh, arg0: str) -> None
        """

    def convertToBiQuadratic(self, info=1):
        """Convert the mesh to a bi-quadratic one.
        For cells that have no bi-quadratic version, the quadratic version is used.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            ParallelMesh: the bi-quadratic mesh.
        """

    def convertToCubic(self, info=1):
        """Convert the mesh to a cubic one.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            Mesh: the cubic mesh.
        """

    def convertToLinear(self, info=1):
        """Convert the mesh to a linear one.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            ParallelMesh: the linearized mesh.
        """

    def convertToQuadratic(self, info=1):
        """Convert the mesh to a quadratic one.

        Arguments:
            info (int) : verbosity mode (1 or 2). Default 1.

        Returns:
            ParallelMesh: the quadratic mesh.
        """

    def fix(
        self,
        remove_orphan=True,
        positive_measure=True,
        outward_normal=True,
        double_nodes=True,
        double_cells=True,
        tole=1e-07,
        info=1,
    ):
        """Fix potential problems.

        Arguments:
            remove_orphan (bool) : remove orphelan nodes.
            positive_measure (bool) : reorder nodes to have a positive measure of cells.
            outward_normal (bool) : reorder nodes to have an outward normal for boundary faces.
            double_nodes (bool) : merge double nodes with almost same coordinates.
            double_cells (bool) : merge double cells with same nodes.
            tole (float) : tolerance for double nodes
            info (int) : verbosity mode (0 or 1 or 2).

        Returns:
            Mesh: fixed mesh
        """

    def getAllMedCellsTypes(self):
        """Return all Med types available in mesh (for all processors).

        Returns:
            list[int]: List of Med types.
        """

    def getCells(self, *args, **kwargs):
        """Overloaded function.

        1. getCells(self: libaster.ParallelMesh, group_name: str) -> list[int]


        Return the list of the indexes of the cells that belong to a group of cells.

        Arguments:
            group_name (str): Name of the local group.

        Returns:
            list[int]: Indexes of the cells of the local group.


        2. getCells(self: libaster.ParallelMesh, groups_name: list[str] = []) -> list[int]


        Return the list of the indexes of the cells that belong to the groups of cells.

        Arguments:
            groups_name (str): Name of the local groups.

        Returns:
            list[int]: Indexes of the cells of the local groups.
        """

    def getCellsOwner(self):
        """Return the rank of the processor which owns the cells

        Returns:
            list[int]: MPI-Rank of the owners of the cells
        """

    def getCellsRanks(self):
        """Return the rank of the sub-domains which have the cells.
        The first subdomain given for a cell is its owner.

        Returns:
            list[list[int]]: MPI-Rank of of the subdomains
        """

    def getGlobalToLocalNodeIds(self):
        """Returns global to local IDs mapping for nodes

        Returns:
            dict[int]: global to local IDs mapping.
        """

    def getGroupsOfCells(self, local=False):
        """Return the list of the existing (local or global) groups of cells.

        Arguments:
            local (bool): search in local or global groups

        Returns:
            list[str]: List of (local or global) groups names (stripped).
        """

    def getGroupsOfNodes(self, local=False):
        """Return the list of the existing (local or global) groups of nodes.

        Arguments:
            local (bool): search in local or global groups

        Returns:
            list[str]: List of (local or global) groups names (stripped).
        """

    def getInnerCells(self):
        """Return the list of the indexes of the inner cells in the mesh

        Returns:
            list[int]: Indexes of the cells.
        """

    def getInnerNodes(self):
        """Return the list of the indexes of the inner nodes in the mesh

        Returns:
            list[int]: Indexes of the nodes.
        """

    def getLastGhostsLayer(self):
        """Return ids in local numbering of ghost nodes on the last layer

        Returns:
            list[int]: List of Nodes ids.
        """

    def getNodesOwner(self):
        """Return the rank of the processor which owns the nodes

        Returns:
            list[int]: MPI-Rank of the owner of the nodes
        """

    def getNodesRanks(self):
        """Return the rank of the sub-domains which have the nodes.
        The first subdomain given for a node is its owner.

        Returns:
            list[list[int]]: MPI-Rank of of the subdomains
        """

    def getOppositeDomains(self):
        """Returns the list of opposite domains of local process"""

    def getOuterCells(self):
        """Return the list of the indexes of the outer cells in the mesh

        Returns:
            list[int]: Indexes of the cells.
        """

    def getOuterNodes(self):
        """Return the list of the indexes of the outer nodes in the mesh

        Returns:
            list[int]: Indexes of the nodes.
        """

    def getReceiveJoint(self, rank):
        """Returns ids of nodes in joint (inner nodes) for an opposite process

        Arguments:
            rank: Rank of opposite domain
        """

    def getSendJoint(self, rank):
        """Returns ids of nodes in joint (inner nodes) for an opposite process

        Arguments:
            rank: Rank of opposite domain
        """

    def hasGroupOfCells(self, group_name, local=False):
        """The global group exists in the mesh

        Arguments:
            group_name (str): Name of the global group.
            local (bool): search in local or global groups

        Returns:
            bool: *True* if exists, *False* otherwise.
        """

    def hasGroupOfNodes(self, group_name, local=False):
        """The (local or global) group exists in the mesh

        Arguments:
            group_name (str): Name of the (local or global) group.
            local (bool): search local or global groups

        Returns:
            bool: *True* if exists, *False* otherwise.
        """

    def isQuadratic(self, local=False):
        """Tells if the mesh contains quadratic cells.

        Arguments:
            local (bool): if *True* only local cells are checked.

        Returns:
            bool: *True* if the mesh contains quadratic cells, *False* otherwise.
        """

    def printMedFile(self, fileName, local=True, version=[0, 0, 0], mode="a"):
        """Print the mesh in the MED format

        Arguments:
            filename (Path|str): Name of the file
            local (bool=True) : print local values only (relevant for a ParallelMesh only)
            version (list): list of size 3 ([major, minor, release])
            mode (str): 'a' or 'w' ('w' to overwrite file)

        Returns:
            Bool: True if of
        """

    def setGroupOfCells(self, group_name, cell_ids):
        """Set new group of cells in the mesh

        Arguments:
            group_name (str): Name of the new group.
            cell_ids (list[int]) : cell ids which are in the group
        """

    def setGroupOfNodes(self, group_name, node_ids, localNumbering=False):
        """Set new group of nodes in the mesh

        Arguments:
            group_name (str): Name of the new group.
            node_ids (list[int]) : node ids which are in the group
            localNumbering=false (bool): ids are given in the local numbering ?
        """

    def setLastGhostsLayer(self, node_ids):
        """Set ids in local numbering of ghost nodes on the last layer

        Arguments:
            list[int]: List of ghost nodes ids.
        """


# class ParallelEquationNumbering in libaster


class ParallelEquationNumbering(EquationNumbering):
    pass

    # Method resolution order:
    #     ParallelEquationNumbering
    #     EquationNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelEquationNumbering) -> None

        2. __init__(self: libaster.ParallelEquationNumbering, arg0: str) -> None
        """

    def getDOFsWithDescription(self, *args, **kwargs):
        """Overloaded function.

        1. getDOFsWithDescription(self: libaster.ParallelEquationNumbering, cmps: list[str] = [], groupNames: list[str] = [], local: bool = True, same_rank: int = <PythonBool.NONE: -1>) -> tuple[tuple[list[int], list[str]], list[int]]


        Get the dofs associated to the given component restricted to the given group.

        Arguments:
            cmps (list[str]): components to extract.
            groupNames (list[str]): group names to filter.
            local (bool): if True use local dof index else use global index in HPC
            same_rank : - None: keep all nodes (default: None)
                        - True: keep the nodes which are owned by the current MPI-rank
                        - False: keep the nodes which are not owned by the current MPI-rank

        Returns:
            pair[list[int], list[str]]: list of nodes and list of components
            list[int]: list of dofs


        2. getDOFsWithDescription(self: libaster.ParallelEquationNumbering, cmps: list[str] = [], nodes: list[int] = [], local: bool = True, same_rank: int = <PythonBool.NONE: -1>) -> tuple[tuple[list[int], list[str]], list[int]]


        Get the dofs associated to the given component restricted to the given nodes.

        Arguments:
            cmps (list[str]): components to extract.
            nodes (list[int]): list of nodes to filter.
            local (bool): if True use local dof index else use global index in HPC
            same_rank : - None: keep all nodes (default: None)
                        - True: keep the nodes which are owned by the current MPI-rank
                        - False: keep the nodes which are not owned by the current MPI-rank

        Returns:
            pair[list[int], list[str]]: list of nodes and list of components.
            list[int]: list of dofs.
        """

    def getGhostDOFs(self, local=True, lastLayerOnly=False):
        """Returns the indexes of the ghost DOFs.

        Arguments:
            local (bool): local or global numbering
            lastLayerOnly (bool): last ghosts layer or all

        Returns:
            int: indexes of the ghost DOFs.
        """

    def getLocalToGlobalMapping(self):
        """Returns the mapping from the local to the global number of the DOFs.

        Returns:
            int: global number of the DOF.
        """

    def getNoGhostDOFs(self, local=True):
        """Returns the indexes of the DOFs owned locally (aka not ghost).

        Returns:
            int: indexes of the DOFs owned locally.
        """

    def getNodeAndComponentFromDOF(self, *args, **kwargs):
        """Overloaded function.

        1. getNodeAndComponentFromDOF(self: libaster.ParallelEquationNumbering, local: bool = True) -> list[tuple[int, str]]


        Return the list of node id and name of component for each dofs

        Arguments:
          local (bool) = True: if True use local node index else use global index in HPC
        Returns:
          list[tuple[int, str]] : node id and name of component for each dofs


        2. getNodeAndComponentFromDOF(self: libaster.ParallelEquationNumbering, dof: int, local: bool = True) -> tuple[int, str]


        Return the node id and name of component for given DOF

        Arguments:
          dof (int): DOF index
          local (bool) = True: if True use local node index else use global index in HPC
        Returns:
          tuple[int, str] : node id and name of component
        """

    def getNumberOfDOFs(self, local=False):
        """Returns the number of DOFs.

        Arguments:
            local (bool): local or parallel request

        Returns:
            int: number of DOFs.
        """

    def globalToLocalDOF(self, glob):
        """Returns the local number of a global DOF.

        Arguments:
            glob (int): global DOF number

        Returns:
            int: local number of the DOF.
        """


# class ParallelDOFNumbering in libaster


class ParallelDOFNumbering(BaseDOFNumbering):
    pass

    # Method resolution order:
    #     ParallelDOFNumbering
    #     BaseDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelDOFNumbering) -> None

        2. __init__(self: libaster.ParallelDOFNumbering, arg0: str) -> None

        3. __init__(self: libaster.ParallelDOFNumbering, arg0: str, arg1: libaster.ParallelEquationNumbering, arg2: libaster.Model) -> None
        """

    def getComponentFromDOF(self, dof, local=False):
        """Returns the component name associated to a dof index.

        - If the dof is associated to a physical DOF, the name of the component is returned.

        - If the dof is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, the name of the component which is constrained by the multiplier is
          returned, precedeed by 'LAGR:', e.g. 'LAGR:DX'.

        - If the dof is associated to a Lagrange multiplier DOF for a multipoint-constraint
          (MPC) implying several DOF, 'LAGR:MPC' is returned (since no component can be
          identified).

        Arguments:
            node (int): Index of the node.
            local (bool): dof in local or global numbering

        Returns:
            str: component names.
        """

    def getComponentFromNode(self, node, local=False):
        """Returns the components name associated to a node index.

        Arguments:
            node (int): Index of the node.
            local (bool, optional): local or global numbering of nodes (default: false).

        Returns:
            str: component names.
        """

    def getComponents(self):
        """Returns all the component names assigned in the numbering.

        Returns:
            str: component names.
        """

    def getDictOfLagrangeDOFs(self, local=False):
        """Returns the Rows Associated to the first and second Lagrange Multipliers Dof

        Arguments:
            local (bool, optional): local or global numbering of DOFs (default: false).

        Returns:
            [dict]: {1 : indexes of the first Lagrange multipliers dof,
                     2 : indexes of the second Lagrange multipliers dof }
        """

    def getGhostDOFs(self, local=True):
        """Returns the indexes of the ghost DOFs.

        Arguments:
            local (bool): local or global numbering

        Returns:
            int: indexes of the ghost DOFs.
        """

    def getLagrangeDOFs(self, local=False):
        """Returns the indexes of the Lagrange multipliers dof.

        Arguments:
            local (bool, optional): local or global numbering of DOFs (default: false).

        Returns:
            int: indexes of the Lagrange multipliers dof.
        """

    def getLocalToGlobalMapping(self):
        """Returns the mapping from the local to the global number of the DOFs.

        Returns:
            int: global number of the DOF.
        """

    def getNoGhostDOFs(self, local=True):
        """Returns the indexes of the DOFs owned locally (aka not ghost).

        Arguments:
            local (bool): local or global numbering

        Returns:
            int: indexes of the DOFs owned locally.
        """

    def getNodeAndComponentFromDOF(self, *args, **kwargs):
        """Overloaded function.

        1. getNodeAndComponentFromDOF(self: libaster.ParallelDOFNumbering, local: bool = True) -> list[tuple[int, str]]


        Return the list of node id and name of component for each dofs

        Arguments:
            local (bool) = True: if True use local node index else use global index in HPC
        Returns:
            list[tuple[int, str]] : node id and name of component for each dofs


        2. getNodeAndComponentFromDOF(self: libaster.ParallelDOFNumbering, dof: int, local: bool = True) -> tuple[int, str]


        Return the node id and name of component for given DOF

        Arguments:
            dof (int): DOF index
            local (bool) = True: if True use local node index else use global index in HPC
        Returns:
            tuple[int, str] : node id and name of component
        """

    def getNodeFromDOF(self, dof, local=False):
        """Returns the node index associated to a dof index.

        Arguments:
            dof (int): Index of the dof.
            local (bool, optional): local or global numbering of DOFs (default: false).

        Returns:
            int: index of the dof.
        """

    def getNumberOfDOFs(self, local=False):
        """Returns the number of DOFs.

        Arguments:
            local (bool): local or parallel request

        Returns:
            int: number of DOFs.
        """

    def getPhysicalDOFs(self, local=False):
        """Returns the indexes of the physical dof.

        Arguments:
            local (bool, optional): local or global numbering of DOFs (default: false).

        Returns:
            int: indexes of the physical dof.
        """

    def globalToLocalDOF(self, glob):
        """Returns the local number of a global DOF.

        Arguments:
            glob (int): global DOF number

        Returns:
            int: local number of the DOF.
        """

    def isPhysicalDOF(self, dof, local=False):
        """If the dof is associated to a physical DOF, return True

        If the dof is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, return False

        Arguments:
            dof (int): Index of the dof.
            local (bool, optional): local or global numbering of DOFs (default: false).

        Returns:
            int: index of the dof.
        """

    def localToGlobalDOF(self, loc):
        """Returns the global number of a local DOF.

        Arguments:
            loc (int): local DOF number

        Returns:
            int: global number of the DOF.
        """

    def useLagrangeDOF(self):
        """Lagrange multipliers are used for BC or MPC.

        Returns:
            bool: *True* if used, *False* otherwise.
        """

    def useSingleLagrangeDOF(self):
        """Single Lagrange multipliers are used for BC or MPC.

        Returns:
            bool: *True* if used, *False* otherwise.
        """


# class ParallelMechanicalLoadReal in libaster


class ParallelMechanicalLoadReal(DataStructure):
    pass

    # Method resolution order:
    #     ParallelMechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelMechanicalLoadReal, arg0: libaster.MechanicalLoadReal, arg1: libaster.Model) -> None

        2. __init__(self: libaster.ParallelMechanicalLoadReal, arg0: str, arg1: libaster.MechanicalLoadReal, arg2: libaster.Model) -> None

        3. __init__(self: libaster.ParallelMechanicalLoadReal, arg0: str, arg1: ParallelFiniteElementDescriptor, arg2: libaster.Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getModel(self):
        pass

    def setRebuildParameters(self, syntax, grpNo, grpMa):
        """Set parameters to be able to rebuild object in case of balancing

        Arguments:
            syntax (SyntaxSaver): syntax used to build object
            grpNo (list of strings): list of node groups
            grpMa (list of strings): list of cell groups
        """


# class ParallelMechanicalLoadFunction in libaster


class ParallelMechanicalLoadFunction(DataStructure):
    pass

    # Method resolution order:
    #     ParallelMechanicalLoadFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelMechanicalLoadFunction, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Model) -> None

        2. __init__(self: libaster.ParallelMechanicalLoadFunction, arg0: str, arg1: libaster.MechanicalLoadFunction, arg2: libaster.Model) -> None

        3. __init__(self: libaster.ParallelMechanicalLoadFunction, arg0: str, arg1: ParallelFiniteElementDescriptor, arg2: libaster.Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getModel(self):
        pass


# class ParallelThermalLoadReal in libaster


class ParallelThermalLoadReal(DataStructure):
    pass

    # Method resolution order:
    #     ParallelThermalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelThermalLoadReal, arg0: libaster.ThermalLoadReal, arg1: libaster.Model) -> None

        2. __init__(self: libaster.ParallelThermalLoadReal, arg0: str, arg1: libaster.ThermalLoadReal, arg2: libaster.Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getModel(self):
        pass

    def setRebuildParameters(self, syntax, grpNo, grpMa):
        """Set parameters to be able to rebuild object in case of balancing

        Arguments:
            syntax (SyntaxSaver): syntax used to build object
            grpNo (list of strings): list of node groups
            grpMa (list of strings): list of cell groups
        """


# class ParallelThermalLoadFunction in libaster


class ParallelThermalLoadFunction(DataStructure):
    pass

    # Method resolution order:
    #     ParallelThermalLoadFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelThermalLoadFunction, arg0: libaster.ThermalLoadFunction, arg1: libaster.Model) -> None

        2. __init__(self: libaster.ParallelThermalLoadFunction, arg0: str, arg1: libaster.ThermalLoadFunction, arg2: libaster.Model) -> None
        """

    def getFiniteElementDescriptor(self):
        pass

    def getModel(self):
        pass


# class ParallelFiniteElementDescriptor in libaster


class ParallelFiniteElementDescriptor(FiniteElementDescriptor):
    pass

    # Method resolution order:
    #     ParallelFiniteElementDescriptor
    #     FiniteElementDescriptor
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0, arg1, arg2):
        pass

    def getJointObjectName(self):
        pass

    def getJoints(self):
        """Return the vector of joints between the curent domain and the others subdomains.

        Returns:
            list: joints between subdomains.
        """


# class ParallelContactFEDescriptor in libaster


class ParallelContactFEDescriptor(FiniteElementDescriptor):
    pass

    # Method resolution order:
    #     ParallelContactFEDescriptor
    #     FiniteElementDescriptor
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0, arg1, arg2, arg3, arg4, arg5):
        pass

    def getJointObjectName(self):
        pass

    def getJoints(self):
        """Return the vector of joints between the curent domain and the others subdomains.

        Returns:
            list: joints between subdomains.
        """


# class ParallelContactNew in libaster


class ParallelContactNew(ContactNew):
    pass

    # Method resolution order:
    #     ParallelContactNew
    #     ContactNew
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelContactNew, arg0: str, arg1: libaster.Model, arg2: libaster.ParallelMesh) -> None

        2. __init__(self: libaster.ParallelContactNew, arg0: libaster.Model, arg1: libaster.ParallelMesh) -> None
        """

    def build(self):
        pass

    def getConnectionModel(self):
        pass

    def getParallelFiniteElementDescriptor(self):
        """Return ParallelFiniteElementDescriptor"""

    def isParallel(self):
        """bool: true if parallel contact."""


# class ParallelFrictionNew in libaster


class ParallelFrictionNew(ParallelContactNew):
    pass

    # Method resolution order:
    #     ParallelFrictionNew
    #     ParallelContactNew
    #     ContactNew
    #     DSWithCppPickling
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelFrictionNew, arg0: str, arg1: libaster.Model, arg2: libaster.ParallelMesh) -> None

        2. __init__(self: libaster.ParallelFrictionNew, arg0: libaster.Model, arg1: libaster.ParallelMesh) -> None
        """


# class ParallelContactPairing in libaster


class ParallelContactPairing(ContactPairing):
    pass

    # Method resolution order:
    #     ParallelContactPairing
    #     ContactPairing
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ParallelContactPairing, arg0: str, arg1: libaster.ParallelContactNew) -> None

        2. __init__(self: libaster.ParallelContactPairing, arg0: libaster.ParallelContactNew) -> None
        """

    def buildFiniteElementDescriptor(self):
        pass

    def getParallelFiniteElementDescriptor(self):
        """Return ParallelFiniteElementDescriptor"""


# built-in function applyFactorOnSubBlocks in libaster


def applyFactorOnSubBlocks(pyFctMat, pyRHS, pyISet, pyJSet=None):
    """Given a MUMPS factor matrix and a sparse PETSc Mat as the right-hand side,
            this routine solves for all columns. If an index set I is provided, it first
            determines the set J of row indices where the RHS has nonzero entries, and
            returns a sparse MATAIJ solution containing only the IxJ block. Supplying I
            can speed up the solve, but is most effective when the RHS column count is
            comparable to or larger than |J|.
    Arguments:
        pyFctMat (PETSc.Mat): the MUMPS factor matrix.
        pyRHS (PETSc.Mat): the sparse MATSEQAIJ holding the multiple right hand sides
        pyISet (PETSc.IS): the index set for entries to compute
        pyJSet (PETSc.IS): optional off diagonal layout index set
    Outputs:
        outMat: the petsc4py aij matrix matrix holding the solutions
    """


# class ConnectionMesh in libaster


class ConnectionMesh(BaseMesh):
    pass

    # Method resolution order:
    #     ConnectionMesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ConnectionMesh, arg0: libaster.ParallelMesh, arg1: list[str], arg2: list[str]) -> None

        2. __init__(self: libaster.ConnectionMesh, arg0: str, arg1: libaster.ParallelMesh, arg2: list[str], arg3: list[str]) -> None
        """

    def getCells(self, group_name=""):
        """Return the list of the indexes of the cells that belong to a group of cells.

        Arguments:
            group_name (str): Name of the local group.

        Returns:
            list[int]: Indexes of the cells of the local group.
        """

    def getGroupsOfCells(self, local=False):
        """Return the list of the existing groups of cells.

        Returns:
            list[str]: List of groups names (stripped).
        """

    def getGroupsOfNodes(self, local=False):
        """Return the list of the existing groups of nodes.

        Returns:
            list[str]: List of groups names (stripped).
        """

    def getNodesGlobalNumbering(self):
        """Return a tuple of the nodes of the mesh with a global numbering

        Returns:
            tuple[int]: list of nodes with global numbering
        """

    def getNodesLocalNumbering(self):
        """Return a tuple of the nodes of the mesh with a local numbering.
        The local numbering is the one coming from the owner of the node,
        hence some nodes can have the same local numbering

        Returns:
            tuple[int]: list of nodes with local numbering
        """

    def getParallelMesh(self):
        """Return a pointer to the ParallelMesh used to built it.

        Returns:
            ParallelMeshPtr: pointer to the ParallelMesh
        """

    def hasGroupOfCells(self, name, local=False):
        """Allows to know if the given group of cells is present in the mesh

        Arguments:
            name (str): name of the group of cell

        Returns:
            bool: True if the group is present
        """

    def hasGroupOfNodes(self, name, local=False):
        """Allows to know if the given group of nodes is present in the mesh

        Arguments:
            name (str): name of the group of nodes

        Returns:
            bool: True if the group is present
        """

    def isConnection(self):
        """Function to know if a mesh is a ConnectionMesh"""


# class ResultNaming in libaster


class ResultNaming:
    pass

    # Method resolution order:
    #     ResultNaming
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    # ----------------------------------------------------------------------
    # Static methods defined here:


# class ListOfFloats in libaster


class ListOfFloats(DataStructure):
    pass

    # Method resolution order:
    #     ListOfFloats
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ListOfFloats) -> None

        2. __init__(self: libaster.ListOfFloats, arg0: str) -> None
        """

    def getValues(self):
        pass

    def setVectorValues(self, arg0):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def size(self):
        pass


# class ListOfIntegers in libaster


class ListOfIntegers(DataStructure):
    pass

    # Method resolution order:
    #     ListOfIntegers
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ListOfIntegers) -> None

        2. __init__(self: libaster.ListOfIntegers, arg0: str) -> None
        """

    def getValues(self):
        pass

    def setVectorValues(self, arg0):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def size(self):
        pass


# class EmpiricalModeResult in libaster


class EmpiricalModeResult(Result):
    pass

    # Method resolution order:
    #     EmpiricalModeResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.EmpiricalModeResult) -> None

        2. __init__(self: libaster.EmpiricalModeResult, arg0: str) -> None
        """


# class EvolutionParameter in libaster


class EvolutionParameter:
    pass

    # Method resolution order:
    #     EvolutionParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, result, fieldName):
        """Constructor of object

        Arguments:
            result (TransientResult): transient result to define external state variable
            fieldName (str): field in transient result to define external state variable
        """

    def __setstate__(self, arg0):
        pass

    def getFieldName(self):
        pass

    def getLeftExtension(self):
        pass

    def getRightExtension(self):
        pass

    def getTimeFormula(self):
        pass

    def getTimeFunction(self):
        pass

    def getTransientResult(self):
        pass

    def setLeftExtension(self, typeExtension):
        """Set type of the extension to the left of the function to shift the results

        Arguments:
            typeExtension (str): type of extension ('CONSTANT', 'EXCLU', 'LINEAIRE')
        """

    def setRightExtension(self, typeExtension):
        """Set type of the extension to the right of the function to shift the results

        Arguments:
            typeExtension (str): type of extension ('CONSTANT', 'EXCLU', 'LINEAIRE')
        """

    def setTimeFunction(self, *args, **kwargs):
        """Overloaded function.

        1. setTimeFunction(self: libaster.EvolutionParameter, formula: libaster.Formula) -> None


                    Set function to shift results

                    Arguments:
                        formula (Formula): formula


        2. setTimeFunction(self: libaster.EvolutionParameter, function: libaster.Function) -> None


                    Set function to shift results

                    Arguments:
                        function (Function): function
        """


# class ExternalStateVariable in libaster


class ExternalStateVariable:
    pass

    # Method resolution order:
    #     ExternalStateVariable
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ExternalStateVariable, arg0: str, arg1: libaster.BaseMesh) -> None

        2. __init__(self: libaster.ExternalStateVariable, arg0: str, arg1: libaster.BaseMesh, arg2: str) -> None

        3. __init__(self: libaster.ExternalStateVariable, arg0: externVarEnumInt, arg1: libaster.BaseMesh) -> None

        4. __init__(self: libaster.ExternalStateVariable, arg0: externVarEnumInt, arg1: libaster.BaseMesh, arg2: str) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getEvolutionParameter(self):
        pass

    def getField(self):
        """Get the field of values"""

    def getReferenceValue(self):
        pass

    def getTransientResult(self):
        """Get the transient result"""

    def getType(self):
        pass

    def isSetRefe(self):
        pass

    def setEvolutionParameter(self, evolutionParameter):
        """Define evolution parameters for values of external state variable

        Arguments:
            evolutionParameter (EvolutionParameter): object EvolutionParameter to define
        """

    def setField(self, field):
        """Define constant value in time for external state variable

        Arguments:
            field (field): field to define value
        """

    def setReferenceValue(self, value):
        """Set reference value for external state variable

        Arguments:
            value (float): reference value
        """


# class ExternalVariableTraits in libaster


class ExternalVariableTraits:
    pass

    # Method resolution order:
    #     ExternalVariableTraits
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def getExternVarTypeStr(self):
        pass


# class externVarEnumInt in libaster


class externVarEnumInt:
    """Enumeration for external variable."""

    # Method resolution order:
    #     externVarEnumInt
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    ConcreteDrying = 11

    ConcreteHydration = 4

    Corrosion = 2

    Geometry = 1

    Irradiation = 5

    IrreversibleStrain = 3

    Neutral1 = 8

    Neutral2 = 9

    Neutral3 = 10

    NumberOfExternVarTypes = 14

    SteelPhases = 6

    Temperature = 0

    TotalFluidPressure = 12

    Unknown = -1

    VolumetricStrain = 13

    ZircaloyPhases = 7


# class ExternalStateVariablesResult in libaster


class ExternalStateVariablesResult(TransientResult):
    pass

    # Method resolution order:
    #     ExternalStateVariablesResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.ExternalStateVariablesResult) -> None

        2. __init__(self: libaster.ExternalStateVariablesResult, arg0: str) -> None
        """


# built-in function createEnthalpy in libaster


def createEnthalpy(rho_cp_func, beta_func):
    """Integrate the rho_cp function by adding a point at T=0 K to be sure \\
    to always manipulate a positive enthalpy.

    Arguments:
        rhoc_cp_func[Function]: Function of RHO_CP
        beta_func[Function]: Function of BETA to modify (add value at T=0K)
    """


# built-in function petscFinalize in libaster


def petscFinalize():
    """Stops the PETSc interface."""


# built-in function petscInitialize in libaster


def petscInitialize(options=""):
    """Starts the PETSc interface with options.

    Arguments:
        options[str]: PETSc options
    """


# built-in function assemblyMatrixToPetsc in libaster


def assemblyMatrixToPetsc(*args, **kwargs):
    """Overloaded function.

    1. assemblyMatrixToPetsc(matr: libaster.AssemblyMatrixDisplacementReal, local: bool) -> object


    Convert a *AssemblyMatrix* object to a PETSc *Mat* object.

    Arguments:
        matr (*AssemblyMatrix*): code_aster matrix.
        local (*bool*): extract only the sequential matrix of the subdomain or the global parallel
                        matrix

    Returns:
        *Mat*: PETSc matrix.


    2. assemblyMatrixToPetsc(matr: libaster.AssemblyMatrixTemperatureReal, local: bool) -> object


    Convert a *AssemblyMatrix* object to a PETSc *Mat* object.

    Arguments:
        matr (*AssemblyMatrix*): code_aster matrix.
        local (*bool*): extract only the sequential matrix of the subdomain or the global parallel
                        matrix

    Returns:
        *Mat*: PETSc matrix.
    """


# class BehaviourProperty in libaster


class BehaviourProperty(DataStructure):
    pass

    # Method resolution order:
    #     BehaviourProperty
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.BehaviourProperty) -> None

        2. __init__(self: libaster.BehaviourProperty, arg0: str) -> None

        3. __init__(self: libaster.BehaviourProperty, arg0: libaster.Model, arg1: libaster.MaterialField) -> None

        4. __init__(self: libaster.BehaviourProperty, arg0: str, arg1: libaster.Model, arg2: libaster.MaterialField) -> None
        """

    def getBehaviourField(self):
        """Return a pointer to the field for behaviour.

        Returns:
            ConstantFieldOnCellsChar16Ptr: behaviour.
        """

    def getConvergenceCriteria(self):
        """Return a pointer to the field for convergence criteria.

        Returns:
            ConstantFieldOnCellsRealPtr: convergence criteria.
        """

    def getMaterialField(self):
        """Return a pointer to the material field.

        Returns:
            MaterialFieldPtr: material field setted.
        """

    def getModel(self):
        """Return a pointer to the model.

        Returns:
            ModelPtr: model setted.
        """

    def getMultipleBehaviourField(self):
        """Return a pointer to the field for multiple behaviour like cristals.

        Returns:
            ConstantFieldOnCellsChar16Ptr: multiple behaviour.
        """

    def hasAnnealing(self):
        """Returns a flag if annealing post-processing is enabled

        Returns:
            bool: *True* if annealing is enabled, *False* otherwise.
        """

    def hasBehaviour(self, behaviour):
        """Return True if the given behaviour name is present.

        Arguments:
            behaviour (str): behaviour name

        Returns:
            bool: *True* if present, *False* otherwise.
        """

    def variNameToCmp(self, cells, variNames):
        """Return the components ("V1", "V1", "V3", etc.) corresponding to variNames for a given set of cells

        Arguments:
            cells (list[int]): the cells
            variNames (list[string]): the variNames

        Returns:
            list[string]: return the components
        """


# class CodedMaterial in libaster


class CodedMaterial:
    pass

    # Method resolution order:
    #     CodedMaterial
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.CodedMaterial, arg0: libaster.MaterialField, arg1: libaster.Model) -> None

        2. __init__(self: libaster.CodedMaterial, arg0: str, arg1: libaster.MaterialField, arg2: libaster.Model) -> None
        """

    def allocate(self, force=False):
        pass

    def constant(self):
        pass

    def getCodedMaterialField(self):
        pass


# built-in function setFortranLoggingLevel in libaster


def setFortranLoggingLevel(level):
    """Set level of logging for fortran code.

    Arguments:
        level[int]: Level of logging
    """


# built-in function resetFortranLoggingLevel in libaster


def resetFortranLoggingLevel():
    """Reset level of logging for fortran code (level = 0)."""


# class PhysicalSolutionRestitutor in libaster


class PhysicalSolutionRestitutor:
    pass

    # Method resolution order:
    #     PhysicalSolutionRestitutor
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.PhysicalSolutionRestitutor, resgen: libaster.TransientGeneralizedResult) -> None


        PhysicalSolutionRestitutor(resgen)

        Create a PhysicalSolutionRestitutor associated with a transient generalized result.

        Arguments:
            resgen (TransientGeneralizedResult):
                Pointer to the transient generalized result from which physical fields
                will be reconstructed.

        Returns:
            PhysicalSolutionRestitutor:
                A restitution object initialized with the given transient result.


        2. __init__(self: libaster.PhysicalSolutionRestitutor, resgen: libaster.TransientGeneralizedResult, nbatch: int, ar: int) -> None


        PhysicalSolutionRestitutor(resgen, nbatch, ar)

        Create a PhysicalSolutionRestitutor with control over batching and activation parameters.

        Arguments:
            resgen (TransientGeneralizedResult):
                Transient generalized result object containing modal information.
            nbatch (int):
                Number of batches used during restitution computations.
            ar (int):
                Additional parameter for enveloppe computations.

        Returns:
            PhysicalSolutionRestitutor:
                Configured restitution object for field reconstruction.
        """

    def computeMaxForFieldsOnCells(self):
        """Compute the time-maximum of all modal fields defined on cells.

        Similar to `computeMaxForFieldsOnNodes`, but applied to cell-based fields.
        Each entry in the returned dictionary corresponds to a field name and its
        cell field containing maximum values over the transient duration.

        Returns:
            dict[str, FieldOnCellsReal]:
                Mapping between field names and their corresponding
                cell fields containing the maximum values over time.
        """

    def computeMaxForFieldsOnNodes(self):
        """Compute the time-maximum of all modal fields defined on nodes.

        This function processes all nodal fields associated with the transient result
        and returns the maximum (component-wise or field-wise) observed over time
        for each field.

        Returns:
            dict[str, FieldOnNodesReal]:
                Mapping between field names and their corresponding
                nodal fields containing the maximum values over time.
        """

    def get_acceleration_coeffs(self):
        """Return the generalized acceleration coefficients.

        These coefficients represent the modal acceleration amplitudes
        used to reconstruct physical acceleration fields.

        Returns:
            list[float]:
                Reference to the vector of acceleration coefficients (no copy).
        """

    def get_displacement_coeffs(self):
        """Return the generalized displacement coefficients.

        These coefficients represent the modal displacement amplitudes
        used during physical restitution.

        Returns:
            list[float]:
                Reference to the vector of displacement coefficients (no copy).
        """

    def get_velocity_coeffs(self):
        """Return the generalized velocity coefficients.

        These coefficients represent the modal velocity amplitudes
        used to reconstruct physical velocity fields.

        Returns:
            list[float]:
                Reference to the vector of velocity coefficients (no copy).
        """


# class PostProcessing in libaster


class PostProcessing:
    pass

    # Method resolution order:
    #     PostProcessing
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0):
        pass

    def computeAnnealing(self, internVar, time_prev, time_curr, externVarPrev, externVarCurr):
        """Modification of internal state variables for annealing

        Arguments:
            internVar (FieldOnNodesReal): internal state variables before annealing
            time_prev (float): time at begin of the step
            time_curr (float): time at end of the step
            externVarPrev (FieldOnCellsReal): external state variables at previous time
            externVarCurr (FieldOnCellsReal): external state variables at current time

        Returns:
            FieldOnCellReals: internal state variables after annealing
        """

    def computeHydration(self, temp_prev, temp_curr, time_prev, time_curr, hydr_prev):
        """Compute hydration at quadrature points (HYDR_ELGA)

        Arguments:
            temp_prev (FieldOnNodesReal): temperature field at begin of current time step
            temp_curr (FieldOnNodesReal): temperature field at end of current time step
            time_prev (float): time at begin of the step
            time_curr (float): time at end of the step
            hydr_prev (FieldOnCellReals): hydration field at begin of current time step

        Returns:
            FieldOnCellReals: hydration field at end of current time step
        """

    def computeMaxResultantForPipe(self, result, field_name):
        """Computes the maximum of the EFGE_ELNO or EGRU_ELNO field in absolute value,
        based on the maximal values of the equivalent moment at each element.

        Arguments:
         result (Result) : ResultPtr
            The result object containing the fields
         field_name (str) : It should be 'EFGE_ELNO' or 'EGRU_ELNO'

        Returns:
         FieldOnCellReals: The maximal value of the field
        """

    def computeStress(self, displ, time=0.0, externVar=None, strx_elga=None):
        """Compute stress SIEF_ELGA

        Arguments:
            displ (FieldOnNodesReal): displacement
            time (float): time
            externVar (FieldOnCellsReal): external state variables
            strx_elga (FieldOnCellsReal): STRX_ELGA field

        Returns:
            FieldOnCellReals: stress SIEF_ELGA field
        """

    def computeStructuralStress(self, displ, time=0.0, externVar=None):
        """Compute stress STRX_ELGA

        Arguments:
            displ (FieldOnNodesReal): displacement
            time (float): time
            externVar (FieldOnCellsReal): external state variables

        Returns:
            FieldOnCellReals: stress STRX_ELGA field
        """


# class HHO in libaster


class HHO:
    pass

    # Method resolution order:
    #     HHO
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.HHO, arg0: libaster.PhysicalProblem) -> None

        2. __init__(self: libaster.HHO, arg0: libaster.Model) -> None
        """

    def __setstate__(self, arg0):
        pass

    def evaluateAtQuadraturePoints(self, hho_field):
        """Evaluate HHO-field at quadrature points

        Arguments:
              hho_field (FieldOnNodesReal): hho field like displacement or thermic

        Returns:
              FieldOnCellsReal: HHO field evaluated at quadrature points (ELGA)
        """

    def getModel(self):
        """Get Model.

        Returns:
              Model: model used for HHO.
        """

    def projectOnHHOCellSpace(self, *args, **kwargs):
        """Overloaded function.

        1. projectOnHHOCellSpace(self: libaster.HHO, field_elga: libaster.FieldOnCellsReal) -> libaster.FieldOnNodesReal


              Project field defined at the quadrature poitns to HHO-cell_space
              Cell space is the restriction of HHO-space to cells only
              Face values are setted to zero

              Arguments:
                    field_elga (FieldOnNodesReal): values of the field at the quadrature poitns

              Returns:
                    FieldOnNodesReal: HHO field


        2. projectOnHHOCellSpace(self: libaster.HHO, func: libaster.GenericFunction, time: float = 0.0) -> libaster.FieldOnNodesReal


              Project real function to HHO Cell-space
              Cell space is the restriction of HHO-space to cells only
              Face values are setted to zero

              Arguments:
                    func (Function): real function to project
                    time (float): time value to evaluate function (default=0.0)

              Returns:
                    FieldOnNodesReal: HHO field


        3. projectOnHHOCellSpace(self: libaster.HHO, func: list[libaster.GenericFunction], time: float = 0.0) -> libaster.FieldOnNodesReal


              Project real function to HHO Cell-space
              Cell space is the restriction of HHO-space to cells only
              Face values are setted to zero

              Arguments:
                    func (Function): real function to project
                    time (float): time value to evaluate function (default=0.0)

              Returns:
                    FieldOnNodesReal: HHO field


        4. projectOnHHOCellSpace(self: libaster.HHO, value: float) -> libaster.FieldOnNodesReal


              Project real value to HHO Cell-space
              Cell space is the restriction of HHO-space to cells only
              Face values are setted to zero

              Arguments:
                    value (float): value to project

              Returns:
                    FieldOnNodesReal: HHO field


        5. projectOnHHOCellSpace(self: libaster.HHO, value: list[float]) -> libaster.FieldOnNodesReal


              Project real value to HHO Cell-space
              Cell space is the restriction of HHO-space to cells only
              Face values are setted to zero

              Arguments:
                    value (float): value to project

              Returns:
                    FieldOnNodesReal: HHO field
        """

    def projectOnHHOSpace(self, *args, **kwargs):
        """Overloaded function.

        1. projectOnHHOSpace(self: libaster.HHO, H1_field: libaster.FieldOnNodesReal) -> libaster.FieldOnNodesReal


              Project field from Lagrange-space to HHO-space

              Arguments:
                    H1_field (FieldOnNodesReal): Lagrange field

              Returns:
                    FieldOnNodesReal: HHO field


        2. projectOnHHOSpace(self: libaster.HHO, func: libaster.GenericFunction, time: float = 0.0) -> libaster.FieldOnNodesReal


              Project real function to HHO-space

              Arguments:
                    func (Function): real function to project
                    time (float): time value to evaluate function (default=0.0)

              Returns:
                    FieldOnNodesReal: HHO field


        3. projectOnHHOSpace(self: libaster.HHO, func: list[libaster.GenericFunction], time: float = 0.0) -> libaster.FieldOnNodesReal


              Project real function to HHO-space

              Arguments:
                    func (Function): real function to project
                    time (float): time value to evaluate function (default=0.0)

              Returns:
                    FieldOnNodesReal: HHO field


        4. projectOnHHOSpace(self: libaster.HHO, value: float) -> libaster.FieldOnNodesReal


              Project real value to HHO-space

              Arguments:
                    value (float): value to project

              Returns:
                    FieldOnNodesReal: HHO field


        5. projectOnHHOSpace(self: libaster.HHO, value: list[float]) -> libaster.FieldOnNodesReal


              Project real value to HHO-space

              Arguments:
                    value (float): value to project

              Returns:
                    FieldOnNodesReal: HHO field
        """

    def projectOnLagrangeSpace(self, hho_field, groupOfCells=[], average=True):
        """Project field from HHO-space to Lagrange-space

        Arguments:
              hho_field (FieldOnNodesReal): hho field like displacement or thermic
              groupOfCells (list[str]): groups where to compute
              average (bool): average or not the field at nodes.

        Returns:
              FieldOnNodesReal: HHO field project on Lagrange space
        """

    def static_condensation(self, matr_elem, vect_elem):
        """Performs static condensation.

        Arguments:
              matr_elem (ElementaryMatrixDisplacementReal): elementary (symetric) matrix.
              vect_elem (ElementaryVectorDisplacementReal): elementary vector.

        Returns:
              [
              [AssemblyMatrixDisplacementReal, FieldOnNodesReal],
              [AssemblyMatrixDisplacementReal, FieldOnNodesReal]
              ]: return two pairs of a matrix and a rhs. First pair is the condensated system
              to solve. Second pair is used for static decondensation.
        """

    def static_decondensation(self, mD, lD, uF):
        """Performs static decondensation. Update cell DoFs.

        Arguments:
              mD (AssemblyMatrixDisplacementReal): matrix of decondensation.
              lD (FieldOnNodesReal): rhs of decondensation.
              uF (FieldOnNodesReal): solution computed after condensation.

        Returns:
              FieldOnNodesReal: solution after decondensation.
        """


# class CommGraph in libaster


class CommGraph:
    pass

    # Method resolution order:
    #     CommGraph
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def addCommunication(self, rank):
        """Add a communication with a process

        Arguments:
            rank: rank of opposite process
        """

    def getMatchings(self):
        """Get matchings of communication graph

        Returns:
            list[int]: list of process to communicate with
        """

    def synchronizeOverProcesses(self):
        """Synchronise graph over processes"""


# class ObjectBalancer in libaster


class ObjectBalancer:
    pass

    # Method resolution order:
    #     ObjectBalancer
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def addElementarySend(self, rank, elemList):
        """Add an elementary send (part of a vector to send to given process)

        Arguments:
            rank: rank of process
            elemList: list of elements to send to the process
        """

    def balanceMedVectorOverProcessesWithRenumbering(self, vector):
        """Balance a med vector of reals over processes

        Arguments:
            vector: list of reals to balance

        Returns:
            MedVector[real]: balanced med vector
        """

    def balanceVectorOverProcesses(self, *args, **kwargs):
        """Overloaded function.

        1. balanceVectorOverProcesses(self: libaster.ObjectBalancer, vector: list[float]) -> list[float]


        Balance a vector of reals over processes

        Arguments:
            vector: list of reals to balance

        Returns:
            list[real]: balanced vector


        2. balanceVectorOverProcesses(self: libaster.ObjectBalancer, vector: list[int]) -> list[int]


        Balance a vector of integers over processes

        Arguments:
            vector: list of integers to balance

        Returns:
            list[int]: balanced vector
        """

    def endElementarySendDefinition(self):
        """End the definition of sends"""

    def getRenumbering(self):
        """Get element renumbering (if necessary)"""

    def prepareCommunications(self):
        """Prepare the communications between processes"""

    def setElementsToKeep(self, elemList):
        """Add a list of elements to keep on local process

        Arguments:
            elemList: list of elements to keep
        """


# class MeshBalancer in libaster


class MeshBalancer:
    pass

    # Method resolution order:
    #     MeshBalancer
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def applyBalancingStrategy(self, vector, out_mesh=None, ghost_layer=1):
        """Apply balancing strategy to given mesh. User must give nodes that local process
        will own (without ghost nodes).
        This function returns a ParallelMesh with joints, ghosts and so on.

        Arguments:
            vector: list of nodes to get on local process
            outMesh: out mesh (optional)
            ghost_layer: ghost layer size (optional)

        Returns:
            mesh: ParallelMesh
        """

    def buildFromBaseMesh(self, mesh):
        """Build balancer on an IncompleteMesh or a Mesh

        Arguments:
            mesh: mesh to balance
        """

    def getCellObjectBalancer(self):
        """Get on cells object balancer

        Returns:
            balancer: object balancer
        """

    def getNodeObjectBalancer(self):
        """Get on nodes object balancer

        Returns:
            balancer: object balancer
        """


# class IncompleteMesh in libaster


class IncompleteMesh(Mesh):
    pass

    # Method resolution order:
    #     IncompleteMesh
    #     Mesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.IncompleteMesh) -> None

        2. __init__(self: libaster.IncompleteMesh, arg0: str) -> None
        """

    def debugCheckFromBaseMesh(self, arg0):
        pass

    def getNodesFromGroup(self, grpName):
        """Get node ids (local numbering) of nodes in a group

        Arguments:
            grpName (str) : node group name

        Returns:
             list: list of ids in local numbering
        """


# class PtScotchPartitioner in libaster


class PtScotchPartitioner:
    pass

    # Method resolution order:
    #     PtScotchPartitioner
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def buildGraph(self, *args, **kwargs):
        """Overloaded function.

        1. buildGraph(self: libaster.PtScotchPartitioner, vertices: list[int], edges: list[int], weights: list[int] = []) -> int


          Build the PtScotch graph from 2 integer vectors (PtScotch format)

          Arguments:
              vertices (list[int]): Gives the position of starts of each vertex connections in edgeloctab
              edges (list[int]): Describes vertex connections (at which vertices each vertex is connected)
              weights (list[int], optional): Vertex weights


        2. buildGraph(self: libaster.PtScotchPartitioner, meshConnectionGraph: MeshConnectionGraph, nodesToGather: list[list[int]] = []) -> int


        Build the PtScotch graph from a MeshConnectionGraph

        Arguments:
            meshConnectionGraph: MeshConnectionGraph
            nodesToGather: list of list of nodes to be gather on same MPI processor
        """

    def checkGraph(self):
        """Ask PtScotch to check the graph"""

    def partitionGraph(self, deterministic=False):
        """Call PtScotch partitioning

        Arguments:
            deterministic (bool=false) : argument to force PtScotch to have a deterministic behaviour

        Returns:
            list[int]: Owner for each nodes
        """

    def writeGraph(self, path):
        """Ask PtScotch to write the graph

        Arguments:
            path: path to output file
        """


# class ParMetisPartitioner in libaster


class ParMetisPartitioner:
    pass

    # Method resolution order:
    #     ParMetisPartitioner
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def buildGraph(self, meshConnectionGraph):
        """Build the ParMetis graph from a MeshConnectionGraph

        Arguments:
            meshConnectionGraph: MeshConnectionGraph
        """

    def partitionGraph(self):
        """Call ParMetis partitioning

        Returns:
            list[int]: Owner for each nodes
        """


# class MeshConnectionGraph in libaster


class MeshConnectionGraph:
    pass

    # Method resolution order:
    #     MeshConnectionGraph
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def buildFromIncompleteMesh(self, mesh):
        """Create the graph corresponding to given IncompleteMesh to be used by PtScotchPartitioner

        Arguments:
            mesh: IncompleteMesh.
        """

    def buildFromIncompleteMeshWithVertexWeights(self, mesh, weights):
        """Create the graph corresponding to given IncompleteMesh to be used by PtScotchPartitioner

        Arguments:
            mesh: IncompleteMesh.
            weights: vertex weights.
        """

    def debugCheck(self):
        """Check graph"""


# built-in function applyBalancingStrategy in libaster


def applyBalancingStrategy(*args, **kwargs):
    """Overloaded function.

    1. applyBalancingStrategy(result: libaster.ElasticResult, vector: list[int]) -> libaster.ElasticResult


    Apply balancing strategy to given result. User must give nodes that local process
    will own (without ghost nodes).
    This function returns a PhysicalProblem with joints, ghosts and so on.

    Arguments:
        result: result to balance
        vector: list of nodes to get on local process

    Returns:
        mesh: PhysicalProblem


    2. applyBalancingStrategy(result: libaster.NonLinearResult, vector: list[int]) -> libaster.NonLinearResult


    Apply balancing strategy to given result. User must give nodes that local process
    will own (without ghost nodes).
    This function returns a PhysicalProblem with joints, ghosts and so on.

    Arguments:
        result: result to balance
        vector: list of nodes to get on local process

    Returns:
        mesh: PhysicalProblem


    3. applyBalancingStrategy(result: libaster.ThermalResult, vector: list[int]) -> libaster.ThermalResult


    Apply balancing strategy to given result. User must give nodes that local process
    will own (without ghost nodes).
    This function returns a PhysicalProblem with joints, ghosts and so on.

    Arguments:
        result: result to balance
        vector: list of nodes to get on local process

    Returns:
        mesh: PhysicalProblem
    """


# built-in function redistributePetscMat in libaster


def redistributePetscMat(pMat, subCommSize):
    """Given a distributed matrix on an MPI communicator,
    this function returns a redistributed matrix on a subcommunicator.

    Arguments:
       pMat: the petsc4py matrix to redistribute.
       subCommSize: the size of the subcommunicator
    Outputs:
       outMat: the redistributed petsc4py matrix on a subcommunicator of size
                subCommSize.
    """


# class MedFileAccessType in libaster


class MedFileAccessType:
    """Enumeration med access type."""

    # Method resolution order:
    #     MedFileAccessType
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __eq__(self, other):
        pass

    def __getstate__(self):
        pass

    def __hash__(self):
        pass

    def __index__(self):
        pass

    def __init__(self, value):
        pass

    def __int__(self):
        pass

    def __ne__(self, other):
        pass

    def __repr__(self):
        pass

    def __setstate__(self, state):
        pass

    def __str__(self):
        pass

    # ----------------------------------------------------------------------
    # Readonly properties defined here:

    @property
    def __members__(self):
        pass

    @property
    def name(self):
        """name(self: object) -> str"""

    @property
    def value(self):
        pass

    # ----------------------------------------------------------------------
    # Data and other attributes defined here:

    MedCreate = 2

    MedReadOnly = 0

    MedReadWrite = 1


# class MedFileReader in libaster


class MedFileReader:
    pass

    # Method resolution order:
    #     MedFileReader
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def close(self):
        """Close med file"""

    def createMesh(self, name, dim, desc):
        """Create new mesh in file

        Arguments:
            name (str): mesh name (length: 64)
            dim (int): mesh dimension
            desc (str): mesh description (length: 200)

        Returns:
            MedMesh: return new med mesh object
        """

    def getField(self, *args, **kwargs):
        """Overloaded function.

        1. getField(self: libaster.MedFileReader, name: str) -> MedField


        Get field from name

        Arguments:
            name (str): field name

        Returns:
            MedField: med field of name name


        2. getField(self: libaster.MedFileReader, iterator: int) -> MedField


        Get field from iterator

        Arguments:
            iterator (int): field iterator

        Returns:
            MedField: med field
        """

    def getFieldNames(self):
        """Get all field names

        Returns:
            list: list of field names
        """

    def getFieldNumber(self):
        """Get field number in field

        Returns:
            int: field number
        """

    def getMesh(self, iterator):
        """Get mesh from iterator

        Arguments:
            iterator (int): iterator on mesh

        Returns:
            MedMesh: med mesh
        """

    def getMeshNumber(self):
        """Get mesh number

        Returns:
            int: mesh number
        """

    def getProfileNumber(self):
        """Get profile number

        Returns:
            int: profile number
        """

    def open(self, path, accessType, version=[0, 0, 0]):
        """Open med file

        Arguments:
            path (Path|str): path to med file
            accessType (MedFileAccessType): med access type
            version (list): list of size 3 ([major, minor, release])

        Returns:
            int: return code (0 if open is ok)
        """

    def openParallel(self, path, accessType):
        """Open med file in parallel

        Arguments:
            path (Path|str): path to med file
            accessType (MedFileAccessType): med access type

        Returns:
            int: return code (0 if open is ok)
        """


# class MedField in libaster


class MedField:
    pass

    # Method resolution order:
    #     MedField
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __pickling_disabled__(self):
        pass

    def getAllSupportEntitiesAtSequence(self, numdt, numit):
        """Get list of all entity type and geometric type in calculation sequence

        Arguments:
            numdt (int): time step id
            numit (int): iteration id

        Returns:
            list: list of pair of entity type and geometry type
        """

    def getComponentName(self):
        """Get field component name"""

    def getComponentNumber(self):
        """Get field component number"""

    def getName(self):
        """Get field name"""

    def getProfileNumberAtSequenceOnEntity(self, arg0, arg1, arg2, arg3):
        """Get profile number in calculation sequence for a given entity and geometric type"""

    def getSequence(self, arg0):
        """Get time step id and iteration id for a given sequence id

        Returns:
            list: time step id and iteration id
        """

    def getSequenceNumber(self):
        """Get calculation sequence number"""

    def getTime(self, arg0):
        """Get time for a given sequence id

        Returns:
            float: time
        """

    def getValuesAtSequenceOnCellTypesList(self, numdt, numit, geomtyp):
        """Get cell field values at calculation sequence from geometric type list

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            geomtyp (list): list of geomtric types

        Returns:
            list: values on cells (same sort as list of geomtric types)
        """

    def getValuesAtSequenceOnEntityAndProfile(self, numdt, numit, entity, geometry, iterator):
        """Get field values

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            entity (int): entity type
            geometry (int): geometric type
            iterator (int): iterator on profile

        Returns:
            list: values
        """

    def getValuesAtSequenceOnNodes(self, numdt, numit):
        """Get node field values at calculation sequence

        Arguments:
            numdt (int): time step id
            numit (int): iteration id

        Returns:
            list: values on nodes
        """


# class MedMesh in libaster


class MedMesh:
    pass

    # Method resolution order:
    #     MedMesh
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __pickling_disabled__(self):
        pass

    def addFamily(self, name, num, grps):
        """Add family to mesh

        Arguments:
            name (str): family name
            num (int): family id
            grps (list): group list
        """

    def getCellFamilyAtSequence(self, numdt, numit, type_iterator):
        """Get cell family in calculation sequence for given profile

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            profile_iterator (int): iterator on profile

        Returns:
            list: family id for cells
        """

    def getCellFamilyForGeometricTypeAtSequence(self, numdt, numit, geom_type):
        """Get cell family for calculation sequence and geometric type

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            geom_type (int): geomtric type

        Returns:
            list: family id for cells
        """

    def getCellNumberAtSequence(self, numdt, numit, geomtype_iterator):
        """Get cell number for calculation sequence and geometric type iterator

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            geomtype_iterator (int): iterator on geometric type

        Returns:
            int: cell number
        """

    def getCellNumberForGeometricTypeAtSequence(self, numdt, numit, geomtype):
        """Get cell number for calculation sequence and geometric type

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            geomtype (int): geometric type

        Returns:
            int: cell number
        """

    def getCellTypeAtSequence(self, numdt, numit, geomtype_iterator):
        """Get cell geometric type for calculation sequence and geomtype_iterator

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            geomtype_iterator (int): iterator on geometric type

        Returns:
            int: geometric type
        """

    def getCellTypeNumberAtSequence(self, numdt, numit):
        """Get cell type number for calculation sequence

        Arguments:
            numdt (int): time step id
            numit (int): iteration id

        Returns:
            int: cell type number
        """

    def getConnectivityAtSequence(self, numdt, numit, geomtype_iterator):
        """Get cell connectivity for calculation sequence and geometric type iterator

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            geomtype_iterator (int): iterator on geometric type

        Returns:
            list: cell connectivity
        """

    def getConnectivityForGeometricTypeAtSequence(self, numdt, numit, geomtype):
        """Get cell connectivity for calculation sequence and geometric type

        Arguments:
            numdt (int): time step id
            numit (int): iteration id
            geomtype (int): geometric type

        Returns:
            list: cell connectivity
        """

    def getDimension(self):
        """Get mesh dimension"""

    def getFamilies(self):
        """Get family list

        Returns:
            list: MedFamily list
        """

    def getGeometricTypesAtSequence(self, numdt, numit):
        """Get all cell geometric types

        Arguments:
            numdt (int): time step id
            numit (int): iteration id

        Returns:
            list: cell geometric type list
        """

    def getName(self):
        """Get mesh name"""

    def getNodeFamilyAtSequence(self, numdt, numit):
        """Get node families for calculation sequence

        Arguments:
            numdt (int): time step id
            numit (int): iteration id

        Returns:
            list: node families
        """

    def getNodeNumberAtSequence(self, numdt, numit):
        """Get node number for calculation sequence

        Arguments:
            numdt (int): time step id
            numit (int): iteration id

        Returns:
            int: node number
        """

    def getNodeNumberForGeometricType(self, geotype):
        """Get node number from a geometric type

        Arguments:
            geotype (int): geometric type

        Returns:
            int: node number
        """

    def getSequence(self, seq_iterator):
        """Get calculation sequence

        Arguments:
            seq_iterator (int): iterator on sequence

        Returns:
            list: pair time step id/iterator id
        """

    def getSequenceNumber(self):
        """Get calculation sequence number"""

    def readCoordinates(self, numdt, numit):
        """Get coordinates for calculation sequence

        Arguments:
            numdt (int): time step id
            numit (int): iteration id

        Returns:
            list: coordinates list
        """


# class MedFamily in libaster


class MedFamily:
    pass

    # Method resolution order:
    #     MedFamily
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __pickling_disabled__(self):
        pass

    def getGroups(self):
        """Get list of groups in family"""

    def getId(self):
        """Get family med id"""

    def getName(self):
        """Get family name"""


# class MedVector in libaster


class MedVector:
    pass

    # Method resolution order:
    #     MedVector
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __getstate__(self):
        pass

    def __init__(self, *args, **kwargs):
        """Overloaded function.

        1. __init__(self: libaster.MedVector) -> None

        2. __init__(self: libaster.MedVector, arg0: int) -> None
        """

    def __setstate__(self, arg0):
        pass

    def getComponentName(self):
        """Get component name"""

    def getComponentNumber(self):
        """Get component name"""

    def getComponentVector(self):
        """Get component on element vector"""

    def getCumulatedSizesVector(self):
        """Get cumulated sizes vector

        Returns:
            list: Cumulated sizes for each element
        """

    def getValues(self):
        """Get vector values (WARNING values are owned by MedVector: no copy)

        Returns:
            numpy array: all field values
        """

    def setComponentName(self, arg0):
        """Set component name"""

    def setComponentNumber(self, arg0):
        """Set component number"""

    def setComponentVector(self, arg0):
        """Set component on element vector"""

    def setCumulatedSizesVector(self, arg0):
        """Set cumulated sizes vector"""

    def setValues(self, arg0):
        """Set vector values (WARNING values are owned by MedVector: no copy)"""

    def size(self):
        """Get vector size, ie: number of elements (cells or nodes)"""


# class MedToAsterReader in libaster


class MedToAsterReader:
    pass

    # Method resolution order:
    #     MedToAsterReader
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def readIncompleteMeshFromMedFile(self, mesh, path, mesh_name="", verbosity=0):
        """Open med file

        Arguments:
            IncompleteMesh: return mesh to fill
            path (Path|str): path to med file
            mesh_name (str): mesh name (optional)
            verbosity (int): verbosity (optional)
        """

    def readMeshFromMedFile(self, mesh, path, mesh_name="", verbosity=0):
        """Open med file

        Arguments:
            Mesh: return mesh to fill
            path (Path|str): path to med file
            mesh_name (str): mesh name (optional)
            verbosity (int): verbosity (optional)
        """

    def readParallelMeshFromMedFile(self, mesh, path, mesh_name="", verbosity=0):
        """Open med file

        Arguments:
            ParallelMesh: return mesh to fill
            path (Path|str): path to med file
            mesh_name (str): mesh name (optional)
            verbosity (int): verbosity (optional)
        """


# class AsterToMedWriter in libaster


class AsterToMedWriter:
    pass

    # Method resolution order:
    #     AsterToMedWriter
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self):
        pass

    def __pickling_disabled__(self):
        pass

    def printMesh(self, *args, **kwargs):
        """Overloaded function.

        1. printMesh(self: libaster.AsterToMedWriter, mesh: libaster.Mesh, path: os.PathLike, parallelPrint: bool = False, mesh_name: str = '', mode: str = 'a') -> bool


        Print mesh to med file

        Arguments:
            Mesh: mesh to print
            path (Path|str): path to med file
            parallelPrint (bool): false by default. If true print in one parallel file (optional)
            mesh_name (str): mesh name (optional)
            mode (str): 'a' or 'w' ('w' to overwrite file) (optional)


        2. printMesh(self: libaster.AsterToMedWriter, mesh: libaster.ParallelMesh, path: os.PathLike, parallelPrint: bool = False, mesh_name: str = '', mode: str = 'a') -> bool


        Print mesh to med file

        Arguments:
            Mesh: mesh to print
            path (Path|str): path to med file
            parallelPrint (bool): false by default. If true print in one parallel file (optional)
            mesh_name (str): mesh name (optional)
            mode (str): 'a' or 'w' ('w' to overwrite file) (optional)


        3. printMesh(self: libaster.AsterToMedWriter, mesh: libaster.ConnectionMesh, path: os.PathLike, parallelPrint: bool = False, mesh_name: str = '', mode: str = 'a') -> bool


        Print mesh to med file

        Arguments:
            Mesh: mesh to print
            path (Path|str): path to med file
            parallelPrint (bool): false by default. If true print in one parallel file (optional)
            mesh_name (str): mesh name (optional)
            mode (str): 'a' or 'w' ('w' to overwrite file) (optional)
        """

    def printResult(self, result, path, parallelPrint=False, mode="a"):
        """Print result to med file

        Arguments:
            result: result to print
            path (Path|str): path to med file
            parallelPrint (bool): false by default. If true print in one parallel file (optional)
            mode (str): 'a' or 'w' ('w' to overwrite file) (optional)
        """


# class FieldCharacteristics in libaster


class FieldCharacteristics:
    pass

    # Method resolution order:
    #     FieldCharacteristics
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0):
        pass

    def __pickling_disabled__(self):
        pass

    def getLocalization(self):
        pass

    def getName(self):
        pass

    def getOption(self):
        pass

    def getParameter(self):
        pass

    def getQuantity(self):
        pass


# built-in function getModelings in libaster


def getModelings():
    pass


# class SyntaxSaver in libaster


class SyntaxSaver:
    pass

    # Method resolution order:
    #     SyntaxSaver
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, arg0, arg1, arg2):
        pass


# built-in function projectionAlongDirection in libaster


def projectionAlongDirection(
    type_elem, nb_node, nb_dim, elem_coor, pt_coor, iter_maxi, tole_maxi, proj_dire
):
    """Do the intersection of a node  with a given element along a given direction

    Arguments:
        type_elem (str)         : type of the element
        nb_node (int)           : number of nodes on element
        nb_dim (int)            : dimension of the problem(2 or 3)
        elem_coor (list[float]) : coordinates of nodes of element
        pt_coor (list[flot])    : coordinates of point to project
        iter_maxi (int)         : maximum number of ierations of the Newton algorithm
        tole_maxi (float)       : tolerance of newton algorithm
        proj_dire (list[float]) : direction of projection
        tang_1 (list[float])    : first tangent of local basis for the projection of point on element
        tang_2  (list[float])   : second tangent of local basis for the projection of point on element
    Returns:
        int : 1 if error detected
        float : scalar value that multiply proj_dire to obtain the projected position
        float : first parametric coordinate of projection of point on element
        float : second parametric coordinate of projection of point on element
    """


# class ResultManager in libaster


class ResultManager:
    pass

    # Method resolution order:
    #     ResultManager
    #     pybind11_builtins.pybind11_object
    #     builtins.object

    # Methods defined here:

    def __init__(self, /, *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature."""

    def releaseCurrentResult(self):
        """Release result"""

    def setCurrentResult(self, result):
        """Set result which will be enrich during calculation

        Argument:
            result (Result): result to enrich
        """

    # ----------------------------------------------------------------------
    # Static methods defined here:
