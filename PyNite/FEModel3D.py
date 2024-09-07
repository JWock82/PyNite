# %%
from os import rename
import warnings
from math import isclose


from PyNite.Solvers import solve
from numpy import array, zeros, matmul, divide, subtract, atleast_2d, all

from PyNite.Node3D import Node3D
from PyNite.Material import Material
from PyNite.Section import Section
from PyNite.PhysMember import PhysMember
from PyNite.Spring3D import Spring3D
from PyNite.Member3D import Member3D
from PyNite.Quad3D import Quad3D
from PyNite.Plate3D import Plate3D
from PyNite.LoadCombo import LoadCombo
from PyNite.Mesh import Mesh, RectangleMesh, AnnulusMesh, FrustrumMesh, CylinderMesh
from PyNite import Analysis

# %%
class FEModel3D():
    """A 3D finite element model object. This object has methods and dictionaries to create, store,
       and retrieve results from a finite element model.
    """

    def __init__(self):
        """Creates a new 3D finite element model.
        """
        
        # Initialize the model's various dictionaries. The dictionaries will be prepopulated with
        # the data types they store, and then those types will be removed. This will give us the
        # ability to get type-based hints when using the dictionaries.

        self.nodes = {str:Node3D}          # A dictionary of the model's nodes
        self.nodes.pop(str)
        self.aux_nodes = {str:Node3D}       # A dictionary of the model's auxiliary nodes
        self.aux_nodes.pop(str)
        self.materials = {str:Material}    # A dictionary of the model's materials
        self.materials.pop(str)
        self.sections = {str:Section}      # A dictonary of the model's cross-sections
        self.sections.pop(str)
        self.springs = {str:Spring3D}      # A dictionary of the model's springs
        self.springs.pop(str)
        self.members = {str:PhysMember}    # A dictionary of the model's physical members
        self.members.pop(str)
        self.quads = {str:Quad3D}          # A dictionary of the model's quadiralterals
        self.quads.pop(str)
        self.plates = {str:Plate3D}        # A dictionary of the model's rectangular plates
        self.plates.pop(str)
        self.meshes = {str:Mesh}           # A dictionary of the model's meshes
        self.meshes.pop(str)         
        self.load_combos = {str:LoadCombo}  # A dictionary of the model's load combinations
        self.load_combos.pop(str)
        self._D = {str:[]}                 # A dictionary of the model's nodal displacements by load combination
        self._D.pop(str)
        
        self.solution = None  # Indicates the solution type for the latest run of the model

    @property
    def load_cases(self):
        """Returns a list of all the load cases in the model (in alphabetical order).
        """
        
        # Create an empty list of load cases
        cases = []

        # Step through each node
        for node in self.nodes.values():
            # Step through each nodal load
            for load in node.NodeLoads:
                # Get the load case for each nodal laod
                cases.append(load[2])
        
        # Step through each member
        for member in self.members.values():
            # Step through each member point load
            for load in member.PtLoads:
                # Get the load case for each member point load
                cases.append(load[3])
            # Step through each member distributed load
            for load in member.DistLoads:
                # Get the load case for each member distributed load
                cases.append(load[5])
        
        # Step through each plate/quad
        for plate in list(self.plates.values()) + list(self.quads.values()):
            # Step through each surface load
            for load in plate.pressures:
                # Get the load case for each plate/quad pressure
                cases.append(load[1])

        # Remove duplicates and return the list (sorted ascending)
        return sorted(list(dict.fromkeys(cases)))

    def add_node(self, name, X, Y, Z):
        """Adds a new node to the model.

        :param name: A unique user-defined name for the node. If set to None or "" a name will be
                     automatically assigned.
        :type name: str
        :param X: The node's global X-coordinate.
        :type X: number
        :param Y: The node's global Y-coordinate.
        :type Y: number
        :param Z: The node's global Z-coordinate.
        :type Z: number
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the node added to the model.
        :rtype: str
        """
        
        # Name the node or check it doesn't already exist
        if name:
            if name in self.nodes:
                raise NameError(f"Node name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "N" + str(len(self.nodes))
            count = 1
            while name in self.nodes: 
                name = "N" + str(len(self.nodes) + count)
                count += 1
        
        # Create a new node
        new_node = Node3D(name, X, Y, Z)
        
        # Add the new node to the list
        self.nodes[name] = new_node
        
        # Flag the model as unsolved
        self.solution = None

        #Return the node name
        return name

    def add_auxnode(self, name, X, Y, Z):
        """Adds a new auxiliary node to the model. Together with a member's `i` and `j` nodes, an
        auxiliary node defines the plane in which the member's local z-axis lies, and the side of
        the member the z-axis points toward. If no auxiliary node is specified for a member, PyNite
        uses its own default configuration.

        :param name: A unique user-defined name for the node. If None or "", a name will be
                     automatically assigned.
        :type name: str
        :param X: The global X-coordinate of the node.
        :type X: number
        :param Y: The global Y-coordinate of the node.
        :type Y: number
        :param Z: The global Z-coordinate of the node.
        :type Z: number
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the auxiliary node that was added to the model.
        :rtype: str
        """ 
        
        # Name the node or check it doesn't already exist
        if name:
            if name in self.aux_nodes:
                raise NameError(f"Auxnode name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "AN" + str(len(self.aux_nodes))
            count = 1
            while name in self.aux_nodes: 
                name = "AN" + str(len(self.aux_nodes)+count)
                count += 1
                
        # Create a new node
        new_node = Node3D(name, X, Y, Z)
        
        # Add the new node to the list
        self.aux_nodes[name] = new_node

        # Flag the model as unsolved
        self.solution = None
        
        #Return the node name
        return name 

    def add_material(self, name, E, G, nu, rho, fy=None):
        """Adds a new material to the model.

        :param name: A unique user-defined name for the material.
        :type name: str
        :param E: The modulus of elasticity of the material.
        :type E: number
        :param G: The shear modulus of elasticity of the material.
        :type G: number
        :param nu: Poisson's ratio of the material.
        :type nu: number
        :param rho: The density of the material
        :type rho: number
        :raises NameError: Occurs when the specified name already exists in the model.
        """

        # Name the material or check it doesn't already exist
        if name:
            if name in self.materials:
                raise NameError(f"Material name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "M" + str(len(self.materials))
            count = 1
            while name in self.materials: 
                name = "M" + str(len(self.materials) + count)
                count += 1
                
        # Create a new material
        new_material = Material(name, E, G, nu, rho, fy)
        
        # Add the new material to the list
        self.materials[name] = new_material
        
        # Flag the model as unsolved
        self.solution = None

    def add_section(self, name, section):
        """Adds a cross-section to the model.

        :param name: A unique name for the cross-section.
        :type name: string
        :param section: A 'PyNite' `Section` object.
        :type section: Section
        """

        # Check if the section name has already been used
        if name not in self.sections.keys():
            # Store the section in the `Sections` dictionary
            self.sections[name] = section
        else:
            # Stop execution and notify the user that the section name is already being used
            raise Exception('Cross-section name ' + name + ' already exists in the model.')

    def add_spring(self, name, i_node, j_node, ks, tension_only=False, comp_only=False):
        """Adds a new spring to the model.

        :param name: A unique user-defined name for the member. If None or "", a name will be
                    automatically assigned
        :type name: str
        :param i_node: The name of the i-node (start node).
        :type i_node: str
        :param j_node: The name of the j-node (end node).
        :type j_node: str
        :param ks: The spring constant (force/displacement).
        :type ks: number
        :param tension_only: Indicates if the member is tension-only, defaults to False
        :type tension_only: bool, optional
        :param comp_only: Indicates if the member is compression-only, defaults to False
        :type comp_only: bool, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the spring that was added to the model.
        :rtype: str
        """ 
        
        # Name the spring or check it doesn't already exist
        if name:
            if name in self.springs:
                raise NameError(f"Spring name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "S" + str(len(self.springs))
            count = 1
            while name in self.springs: 
                name = "S" + str(len(self.springs) + count)
                count += 1
                
        # Create a new spring
        new_spring = Spring3D(name, self.nodes[i_node], self.nodes[j_node],
                              ks, self.load_combos, tension_only=tension_only,
                              comp_only=comp_only)
        
        # Add the new spring to the list
        self.springs[name] = new_spring
        
        # Flag the model as unsolved
        self.solution = None

        # Return the spring name
        return name

    def add_member(self, name, i_node, j_node, material_name, Iy=None, Iz=None, J=None, A=None, aux_node=None, tension_only=False, comp_only=False, section_name=None):
        """Adds a new physical member to the model.

        :param name: A unique user-defined name for the member. If ``None`` or ``""``, a name will be automatically assigned
        :type name: str
        :param i_node: The name of the i-node (start node).
        :type i_node: str
        :param j_node: The name of the j-node (end node).
        :type j_node: str
        :param material_name: The name of the material of the member.
        :type material: str
        :param Iy: The moment of inertia of the member about its local y-axis.
        :type Iy: number
        :param Iz: The moment of inertia of the member about its local z-axis.
        :type Iz: number
        :param J: The polar moment of inertia of the member.
        :type J: number
        :param A: The cross-sectional area of the member.
        :type A: number
        :param auxNode: The name of the auxiliary node used to define the local z-axis. The default is ``None``, in which case the program defines the axis instead of using an auxiliary node.
        :type aux_node: str, optional
        :param tension_only: Indicates if the member is tension-only, defaults to False
        :type tension_only: bool, optional
        :param comp_only: Indicates if the member is compression-only, defaults to False
        :type comp_only: bool, optional
        :raises NameError: Occurs if the specified name already exists.
        :param section_name: The name of the cross section to use for section properties. If section is `None` the user must provide `Iy`, `Iz`, `A`, and `J`. If section is not `None` any values of `Iy`, `Iz`, `A`, and `J` will be ignored and the cross-section's values will be used instead.
        :type section: string
        :return: The name of the member added to the model.
        :rtype: str
        """

        # Name the member or check it doesn't already exist
        if name:
            if name in self.members: raise NameError(f"Member name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "M" + str(len(self.members))
            count = 1
            while name in self.members: 
                name = "M" + str(len(self.members)+count)
                count += 1
                
        # Create a new member
        if aux_node == None:
            new_member = PhysMember(name, self.nodes[i_node], self.nodes[j_node], material_name, self, Iy, Iz, J, A, tension_only=tension_only, comp_only=comp_only, section_name=section_name)
        else:
            new_member = PhysMember(name, self.nodes[i_node], self.nodes[j_node], material_name, self, Iy, Iz, J, A, aux_node=self.aux_nodes[aux_node], tension_only=tension_only, comp_only=comp_only, section_name=section_name)
        
        # Add the new member to the list
        self.members[name] = new_member
        
        # Flag the model as unsolved
        self.solution = None

        # Return the member name
        return name

    def add_plate(self, name, i_node, j_node, m_node, n_node, t, material_name, kx_mod=1.0, ky_mod=1.0):
        """Adds a new rectangular plate to the model. The plate formulation for in-plane (membrane)
        stiffness is based on an isoparametric formulation. For bending, it is based on a 12-term
        polynomial formulation. This element must be rectangular, and must not be used where a
        thick plate formulation is needed. For a more versatile plate element that can handle
        distortion and thick plate conditions, consider using the `add_quad` method instead.

        :param name: A unique user-defined name for the plate. If None or "", a name will be
                     automatically assigned.
        :type name: str
        :param i_node: The name of the i-node.
        :type i_node: str
        :param j_node: The name of the j-node.
        :type j_node: str
        :param m_node: The name of the m-node.
        :type m_node: str
        :param n_node: The name of the n-node.
        :type n_node: str
        :param t: The thickness of the element.
        :type t: number
        :param material_name: The name of the material for the element.
        :type material_name: str
        :param kx_mod: Stiffness modification factor for in-plane stiffness in the element's local
                       x-direction, defaults to 1 (no modification).
        :type kx_mod: number, optional
        :param ky_mod: Stiffness modification factor for in-plane stiffness in the element's local
                       y-direction, defaults to 1 (no modification).
        :type ky_mod: number, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the element added to the model.
        :rtype: str
        """
        
        # Name the plate or check it doesn't already exist
        if name:
            if name in self.plates: raise NameError(f"Plate name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "P" + str(len(self.plates))
            count = 1
            while name in self.plates: 
                name = "P" + str(len(self.plates)+count)
                count += 1
        
        # Create a new plate
        new_plate = Plate3D(name, self.nodes[i_node], self.nodes[j_node], self.nodes[m_node],
                           self.nodes[n_node], t, material_name, self, kx_mod, ky_mod)
        
        # Add the new plate to the list
        self.plates[name] = new_plate

        # Flag the model as unsolved
        self.solution = None
        
        # Return the plate name
        return name

    def add_quad(self, name, i_node, j_node, m_node, n_node, t, material_name, kx_mod=1.0, ky_mod=1.0):
        """Adds a new quadrilateral to the model. The quad formulation for in-plane (membrane)
        stiffness is based on an isoparametric formulation. For bending, it is based on an MITC4
        formulation. This element handles distortion relatively well, and is appropriate for thick
        and thin plates. One limitation with this element is that it does a poor job of reporting
        corner stresses. Corner forces, however are very accurate. Center stresses are very
        accurate as well. For cases where corner stress results are important, consider using the
        `add_plate` method instead.

        :param name: A unique user-defined name for the quadrilateral. If None or "", a name will
                     be automatically assigned.
        :type name: str
        :param i_node: The name of the i-node.
        :type i_node: str
        :param j_node: The name of the j-node.
        :type j_node: str
        :param m_node: The name of the m-node.
        :type m_node: str
        :param n_node: The name of the n-node.
        :type n_node: str
        :param t: The thickness of the element.
        :type t: number
        :param material_name: The name of the material for the element.
        :type material_name: str
        :param kx_mod: Stiffness modification factor for in-plane stiffness in the element's local
            x-direction, defaults to 1 (no modification).
        :type kx_mod: number, optional
        :param ky_mod: Stiffness modification factor for in-plane stiffness in the element's local
            y-direction, defaults to 1 (no modification).
        :type ky_mod: number, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the element added to the model.
        :rtype: str
        """
        
        # Name the quad or check it doesn't already exist
        if name:
            if name in self.quads: raise NameError(f"Quad name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "Q" + str(len(self.quads))
            count = 1
            while name in self.quads: 
                name = "Q" + str(len(self.quads) + count)
                count += 1
        
        # Create a new member
        new_quad = Quad3D(name, self.nodes[i_node], self.nodes[j_node], self.nodes[m_node],
                         self.nodes[n_node], t, material_name, self, kx_mod, ky_mod)
        
        # Add the new member to the list
        self.quads[name] = new_quad

        # Flag the model as unsolved
        self.solution = None
        
        #Return the quad name
        return name

    def add_rectangle_mesh(self, name, mesh_size, width, height, thickness, material_name, kx_mod=1.0, ky_mod=1.0, origin=[0, 0, 0], plane='XY', x_control=None, y_control=None, start_node=None, start_element = None, element_type='Quad'):
        """Adds a rectangular mesh of elements to the model.

        :param name: A unique name for the mesh.
        :type name: str
        :param mesh_size: The desired mesh size.
        :type mesh_size: number
        :param width: The overall width of the rectangular mesh measured along its local x-axis.
        :type width: number
        :param height: The overall height of the rectangular mesh measured along its local y-axis.
        :type height: number
        :param thickness: The thickness of each element in the mesh.
        :type thickness: number
        :param material_name: The name of the material for elements in the mesh.
        :type material_name: str
        :param kx_mod: Stiffness modification factor for in-plane stiffness in the element's local x-direction. Defaults to 1.0 (no modification).
        :type kx_mod: float, optional
        :param ky_mod: Stiffness modification factor for in-plane stiffness in the element's local y-direction. Defaults to 1.0 (no modification).
        :type ky_mod: float, optional
        :param origin: The origin of the regtangular mesh's local coordinate system. Defaults to [0, 0, 0]
        :type origin: list, optional
        :param plane: The plane the mesh will be parallel to. Options are 'XY', 'YZ', and 'XZ'. Defaults to 'XY'.
        :type plane: str, optional
        :param x_control: A list of control points along the mesh's local x-axis to work into the mesh. Defaults to `None`.
        :type x_control: list, optional
        :param y_control: A list of control points along the mesh's local y-axis to work into the mesh. Defaults to None.
        :type y_control: list, optional
        :param start_node: The name of the first node in the mesh. If set to `None` the program will use the next available node name. Default is `None`.
        :type start_node: str, optional
        :param start_element: The name of the first element in the mesh. If set to `None` the program will use the next available element name. Default is `None`.
        :type start_element: str, optional
        :param element_type: They type of element to make the mesh out of. Either 'Quad' or 'Rect'. Defaults to 'Quad'.
        :type element_type: str, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the mesh added to the model.
        :rtype: str
        """
        
        # Check if a mesh name has been provided
        if name:
            # Check that the mesh name isn't already being used
            if name in self.meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Rename the mesh if necessary
        else:
            name = self.unique_name(self.meshes, 'MSH')
        
        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.nodes, 'N')
        if element_type == 'Rect' and start_element is None:
            start_element = self.unique_name(self.plates, 'R')
        elif element_type == 'Quad' and start_element is None:
            start_element = self.unique_name(self.quads, 'Q')
        
        # Create the mesh
        new_mesh = RectangleMesh(mesh_size, width, height, thickness, material_name, self, kx_mod,
                                 ky_mod, origin, plane, x_control, y_control, start_node,
                                 start_element, element_type=element_type)

        # Add the new mesh to the `Meshes` dictionary
        self.meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name
    
    def add_annulus_mesh(self, name, mesh_size, outer_radius, inner_radius, thickness, material_name, kx_mod=1.0, 
            ky_mod=1.0, origin=[0, 0, 0], axis='Y', start_node=None, start_element=None):
        """Adds a mesh of quadrilaterals forming an annulus (a donut).

        :param name: A unique name for the mesh.
        :type name: str
        :param mesh_size: The target mesh size.
        :type mesh_size: float
        :param outer_radius: The radius to the outside of the annulus.
        :type outer_radius: float
        :param inner_radius: The radius to the inside of the annulus.
        :type inner_radius: float
        :param thickness: Element thickness.
        :type thickness: float
        :param material_name: The name of the element material.
        :type material_name: str
        :param kx_mod: Stiffness modification factor for radial stiffness in the element's local
                       x-direction. Default is 1.0 (no modification).
        :type kx_mod: float, optional
        :param ky_mod: Stiffness modification factor for meridional stiffness in the element's
                       local y-direction. Default is 1.0 (no modification).
        :type ky_mod: float, optional
        :param origin: The origin of the mesh. The default is [0, 0, 0].
        :type origin: list, optional
        :param axis: The global axis about which the mesh will be generated. The default is 'Y'.
        :type axis: str, optional
        :param start_node: The name of the first node in the mesh. If set to `None` the program
                           will use the next available node name. Default is `None`.
        :type start_node: str, optional
        :param start_element: The name of the first element in the mesh. If set to `None` the
                              program will use the next available element name. Default is `None`.
        :type start_element: str, optional
        :raises NameError: Occurs if the specified name already exists in the model.
        :return: The name of the mesh added to the model.
        :rtype: str
        """
        
        # Check if a mesh name has been provided
        if name:
            # Check that the mesh name doesn't already exist
            if name in self.meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Give the mesh a new name if necessary
        else:
            name = self.unique_name(self.meshes, 'MSH')

        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.nodes, 'N')
        if start_element is None:
            start_element = self.unique_name(self.quads, 'Q')
        
        # Create a new mesh
        new_mesh = AnnulusMesh(mesh_size, outer_radius, inner_radius, thickness, material_name, self,
                               kx_mod, ky_mod, origin, axis, start_node, start_element)

        # Add the new mesh to the `Meshes` dictionary
        self.meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name

    def add_frustrum_mesh(self, name, mesh_size, large_radius, small_radius, height, thickness,material_name, kx_mod=1.0, ky_mod=1.0, origin=[0, 0, 0], axis='Y', start_node=None, start_element=None):
        """Adds a mesh of quadrilaterals forming a frustrum (a cone intersected by a horizontal plane).

        :param name: A unique name for the mesh.
        :type name: str
        :param mesh_size: The target mesh size
        :type mesh_size: number
        :param large_radius: The larger of the two end radii.
        :type large_radius: number
        :param small_radius: The smaller of the two end radii.
        :type small_radius: number
        :param height: The height of the frustrum.
        :type height: number
        :param thickness: The thickness of the elements.
        :type thickness: number
        :param material_name: The name of the element material.
        :type material_name: str
        :param kx_mod: Stiffness modification factor for radial stiffness in each element's local x-direction, defaults to 1 (no modification).
        :type kx_mod: number, optional
        :param ky_mod: Stiffness modification factor for meridional stiffness in each element's local y-direction, defaults to 1 (no modification).
        :type ky_mod: number, optional
        :param origin: The origin of the mesh, defaults to [0, 0, 0].
        :type origin: list, optional
        :param axis: The global axis about which the mesh will be generated, defaults to 'Y'.
        :type axis: str, optional
        :param start_node: The name of the first node in the mesh. If set to None the program will use the next available node name, defaults to None.
        :type start_node: str, optional
        :param start_element: The name of the first element in the mesh. If set to `None` the
                              program will use the next available element name, defaults to None
        :type start_element: str, optional
        :raises NameError: Occurs if the specified name already exists.
        :return: The name of the mesh added to the model.
        :rtype: str
        """

        # Check if a name has been provided
        if name:
            # Check that the mesh name doesn't already exist
            if name in self.meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Give the mesh a new name if necessary
        else:
            name = self.unique_name(self.meshes, 'MSH')

        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.nodes, 'N')
        if start_element is None:
            start_element = self.unique_name(self.quads, 'Q')
        
        # Create a new mesh
        new_mesh = FrustrumMesh(mesh_size, large_radius, small_radius, height, thickness, material_name,
                                self, kx_mod, ky_mod, origin, axis, start_node, start_element)

        # Add the new mesh to the `Meshes` dictionary
        self.meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name
    
    def add_cylinder_mesh(self, name, mesh_size, radius, height, thickness, material_name, kx_mod=1,
                          ky_mod=1, origin=[0, 0, 0], axis='Y', num_elements=None, start_node=None,
                          start_element=None, element_type='Quad'):
        """Adds a mesh of elements forming a cylinder.

        :param name: A unique name for the mesh.
        :type name: str
        :param mesh_size: The target mesh size.
        :type mesh_size: float
        :param radius: The radius of the cylinder.
        :type radius: float
        :param height: The height of the cylinder.
        :type height: float
        :param thickness: Element thickness.
        :type thickness: float
        :param material_name: The name of the element material.
        :type material_name: str
        :param kx_mod: Stiffness modification factor for hoop stiffness in each element's local
                       x-direction. Defaults to 1.0 (no modification).
        :type kx_mod: int, optional
        :param ky_mod: Stiffness modification factor for meridional stiffness in each element's
                       local y-direction. Defaults to 1.0 (no modification).
        :type ky_mod: int, optional
        :param origin: The origin [X, Y, Z] of the mesh. Defaults to [0, 0, 0].
        :type origin: list, optional
        :param axis: The global axis about which the mesh will be generated. Defaults to 'Y'.
        :type axis: str, optional
        :param num_elements: The number of elements to use to form each course of elements. This
                             is typically only used if you are trying to match the nodes to another
                             mesh's nodes. If set to `None` the program will automatically
                             calculate the number of elements to use based on the mesh size.
                             Defaults to None.
        :type num_elements: int, optional
        :param start_node: The name of the first node in the mesh. If set to `None` the program
                           will use the next available node name. Defaults to `None`.
        :type start_node: str, optional
        :param start_element: The name of the first element in the mesh. If set to `None` the
                              program will use the next available element name. Defaults to `None`.
        :type start_element: str, optional
        :param element_type: The type of element to make the mesh out of. Either 'Quad' or 'Rect'.
                             Defaults to 'Quad'.
        :type element_type: str, optional
        :raises NameError: Occurs when the specified mesh name is already being used in the model.
        :return: The name of the mesh added to the model
        :rtype: str
        """            
        
        # Check if a name has been provided
        if name:
            # Check that the mesh name doesn't already exist
            if name in self.meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Give the mesh a new name if necessary
        else:
            name = self.unique_name(self.meshes, 'MSH')

        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.nodes, 'N')
        if element_type == 'Rect' and start_element is None:
            start_element = self.unique_name(self.plates, 'R')
        elif element_type == 'Quad' and start_element is None:
            start_element = self.unique_name(self.quads, 'Q')
        
        # Create a new mesh
        new_mesh = CylinderMesh(mesh_size, radius, height, thickness, material_name, self,
                               kx_mod, ky_mod, origin, axis, start_node, start_element,
                               num_elements, element_type)

        # Add the new mesh to the `Meshes` dictionary
        self.meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name

    def merge_duplicate_nodes(self, tolerance=0.001):
        """Removes duplicate nodes from the model and returns a list of the removed node names.

        :param tolerance: The maximum distance between two nodes in order to consider them
                          duplicates. Defaults to 0.001.
        :type tolerance: float, optional
        """

        # Initialize a dictionary marking where each node is used
        node_lookup = {node_name: [] for node_name in self.nodes.keys()}
        element_dicts = ('springs', 'members', 'plates', 'quads')
        node_types = ('i_node', 'j_node', 'm_node', 'n_node')

        # Step through each dictionary of elements in the model (springs, members, plates, quads)
        for element_dict in element_dicts:

            # Step through each element in the dictionary
            for element in getattr(self, element_dict).values():

                # Step through each possible node type in the element (i-node, j-node, m-node, n-node)
                for node_type in node_types:

                    # Get the current element's node having the current type
                    # Return `None` if the element doesn't have this node type
                    node = getattr(element, node_type, None)

                    # Determine if the node exists on the element
                    if node is not None:
                        # Add the element to the list of elements attached to the node
                        node_lookup[node.name].append((element, node_type))

        # Make a list of the names of each node in the model
        node_names = list(self.nodes.keys())

        # Make a list of nodes to be removed from the model
        remove_list = []

        # Step through each node in the copy of the `Nodes` dictionary
        for i, node_1_name in enumerate(node_names):

            # Skip iteration if `node_1` has already been removed
            if node_lookup[node_1_name] is None:
                continue

            # There is no need to check `node_1` against itself
            for node_2_name in node_names[i + 1:]:

                # Skip iteration if node_2 has already been removed
                if node_lookup[node_2_name] is None:
                    continue

                # Calculate the distance between nodes
                if self.nodes[node_1_name].distance(self.nodes[node_2_name]) > tolerance:
                    continue

                # Replace references to `node_2` in each element with references to `node_1`
                for element, node_type in node_lookup[node_2_name]:
                    setattr(element, node_type, self.nodes[node_1_name])

                # Flag `node_2` as no longer used
                node_lookup[node_2_name] = None

                # Merge any boundary conditions
                support_cond = ('support_DX', 'support_DY', 'support_DZ', 'support_RX', 'support_RY', 'support_RZ')
                for dof in support_cond:
                    if getattr(self.nodes[node_2_name], dof) == True:
                        setattr(self.nodes[node_1_name], dof, True)
                
                # Merge any spring supports
                spring_cond = ('spring_DX', 'spring_DY', 'spring_DZ', 'spring_RX', 'spring_RY', 'spring_RZ')
                for dof in spring_cond:
                    value = getattr(self.nodes[node_2_name], dof)
                    if value != [None, None, None]:
                        setattr(self.nodes[node_1_name], dof, value)
                
                # Fix the mesh labels
                for mesh in self.meshes.values():

                    # Fix the nodes in the mesh
                    if node_2_name in mesh.nodes.keys():

                        # Attach the correct node to the mesh
                        mesh.nodes[node_2_name] = self.nodes[node_1_name]

                        # Fix the dictionary key
                        mesh.nodes[node_1_name] = mesh.nodes.pop(node_2_name)

                    # Fix the elements in the mesh
                    for element in mesh.elements.values():
                        if node_2_name == element.i_node.name: element.i_node = self.nodes[node_1_name]
                        if node_2_name == element.j_node.name: element.j_node = self.nodes[node_1_name]
                        if node_2_name == element.m_node.name: element.m_node = self.nodes[node_1_name]
                        if node_2_name == element.n_node.name: element.n_node = self.nodes[node_1_name]
                    
                # Add the node to the `remove` list
                remove_list.append(node_2_name)

        # Remove `node_2` from the model's `Nodes` dictionary
        for node_name in remove_list:
            self.nodes.pop(node_name)
        
        # Flag the model as unsolved
        self.solution = None

    def delete_node(self, node_name):
        """Removes a node from the model. All nodal loads associated with the node and elements attached to the node will also be removed.

        :param node_name: The name of the node to be removed.
        :type node_name: str
        """
            
        # Remove the node. Nodal loads are stored within the node, so they
        # will be deleted automatically when the node is deleted.
        self.nodes.pop(node_name)
        
        # Find any elements attached to the node and remove them
        self.members = {name: member for name, member in self.members.items() if member.i_node.name != node_name and member.j_node.name != node_name}
        self.plates = {name: plate for name, plate in self.plates.items() if plate.i_node.name != node_name and plate.j_node.name != node_name and plate.m_node.name != node_name and plate.n_node.name != node_name}
        self.quads = {name: quad for name, quad in self.quads.items() if quad.i_node.name != node_name and quad.j_node.name != node_name and quad.m_node.name != node_name and quad.n_node.name != node_name}

        # Flag the model as unsolved
        self.solution = None

    def delete_auxnode(self, auxnode_name):
        """Removes an auxiliary node from the model.

        :param auxnode_name: The name of the auxiliary node to be removed.
        :type auxnode_name: str
        """

        # Remove the auxiliary node
        self.aux_nodes.pop(auxnode_name)

        # Remove the auxiliary node from any members that were using it
        for member in self.members.values():
            if member.auxNode == auxnode_name:
                member.auxNode = None
        
        # Flag the model as unsolved
        self.solution = None

    def delete_spring(self, spring_name):
        """Removes a spring from the model.

        :param spring_name: The name of the spring to be removed.
        :type spring_name: str
        """
        
        # Remove the spring
        self.springs.pop(spring_name)

        # Flag the model as unsolved
        self.solution = None

    def delete_member(self, member_name):
        """Removes a member from the model. All member loads associated with the member will also
           be removed.

        :param member_name: The name of the member to be removed.
        :type member_name: str
        """
        
        # Remove the member. Member loads are stored within the member, so they
        # will be deleted automatically when the member is deleted.
        self.members.pop(member_name)

        # Flag the model as unsolved
        self.solution = None
        
    def def_support(self, node_name, support_DX=False, support_DY=False, support_DZ=False, support_RX=False, support_RY=False, support_RZ=False):
        """Defines the support conditions at a node. Nodes will default to fully unsupported
           unless specified otherwise.

        :param node_name: The name of the node where the support is being defined.
        :type node_name: str
        :param support_DX: Indicates whether the node is supported against translation in the
                           global X-direction. Defaults to False.
        :type support_DX: bool, optional
        :param support_DY: Indicates whether the node is supported against translation in the
                           global Y-direction. Defaults to False.
        :type support_DY: bool, optional
        :param support_DZ: Indicates whether the node is supported against translation in the
                           global Z-direction. Defaults to False.
        :type support_DZ: bool, optional
        :param support_RX: Indicates whether the node is supported against rotation about the
                           global X-axis. Defaults to False.
        :type support_RX: bool, optional
        :param support_RY: Indicates whether the node is supported against rotation about the
                           global Y-axis. Defaults to False.
        :type support_RY: bool, optional
        :param support_RZ: Indicates whether the node is supported against rotation about the
                           global Z-axis. Defaults to False.
        :type support_RZ: bool, optional
        """            
        
        # Get the node to be supported
        node = self.nodes[node_name]
                
        # Set the node's support conditions
        node.support_DX = support_DX
        node.support_DY = support_DY
        node.support_DZ = support_DZ
        node.support_RX = support_RX
        node.support_RY = support_RY
        node.support_RZ = support_RZ

        # Flag the model as unsolved
        self.solution = None

    def def_support_spring(self, node_name, dof, stiffness, direction=None):
        """Defines a spring support at a node.

        :param node_name: The name of the node to apply the spring support to.
        :type node_name: str
        :param dof: The degree of freedom to apply the spring support to.
        :type dof: str ('DX', 'DY', 'DZ', 'RX', 'RY', or 'RZ')
        :param stiffness: The translational or rotational stiffness of the spring support.
        :type stiffness: float
        :param direction: The direction in which the spring can act. '+' allows the spring to resist positive displacements. '-' allows the spring to resist negative displacements. None allows the spring to act in both directions. Default is None.
        :type direction: str or None ('+', '-', None), optional
        :raises ValueError: Occurs when an invalid support spring direction has been specified.
        :raises ValueError: Occurs when an invalid support spring degree of freedom has been specified.
        """        
        
        if dof in ('DX', 'DY', 'DZ', 'RX', 'RY', 'RZ'):
            if direction in ('+', '-', None):
                if dof == 'DX':
                    self.nodes[node_name].spring_DX = [stiffness, direction, True]
                elif dof == 'DY':
                    self.nodes[node_name].spring_DY = [stiffness, direction, True]
                elif dof == 'DZ':
                    self.nodes[node_name].spring_DZ = [stiffness, direction, True]
                elif dof == 'RX':
                    self.nodes[node_name].spring_RX = [stiffness, direction, True]
                elif dof == 'RY':
                    self.nodes[node_name].spring_RY = [stiffness, direction, True]
                elif dof == 'RZ':
                    self.nodes[node_name].spring_RZ = [stiffness, direction, True]
            else:
                raise ValueError('Invalid support spring direction. Specify \'+\', \'-\', or None.')
        else:
            raise ValueError('Invalid support spring degree of freedom. Specify \'DX\', \'DY\', \'DZ\', \'RX\', \'RY\', or \'RZ\'')
        
        # Flag the model as unsolved
        self.solution = None

    def def_node_disp(self, node_name, direction, magnitude): 
        """Defines a nodal displacement at a node.

        :param node_name: The name of the node where the nodal displacement is being applied.
        :type node_name: str
        :param direction: The global direction the nodal displacement is being applied in. Displacements are 'DX', 'DY', and 'DZ'. Rotations are 'RX', 'RY', and 'RZ'.
        :type direction: str
        :param magnitude: The magnitude of the displacement.
        :type magnitude: float
        :raises ValueError: _description_
        """
            
        # Validate the value of direction
        if direction not in ('DX', 'DY', 'DZ', 'RX', 'RY', 'RZ'):
            raise ValueError(f"direction must be 'DX', 'DY', 'DZ', 'RX', 'RY', or 'RZ'. {direction} was given.")
        # Get the node
        node = self.nodes[node_name]

        if direction == 'DX':
            node.EnforcedDX = magnitude
        if direction == 'DY':
            node.EnforcedDY = magnitude
        if direction == 'DZ':
            node.EnforcedDZ = magnitude
        if direction == 'RX':
            node.EnforcedRX = magnitude
        if direction == 'RY':
            node.EnforcedRY = magnitude
        if direction == 'RZ':
            node.EnforcedRZ = magnitude
        
        # Flag the model as unsolved
        self.solution = None

    def def_releases(self, member_name, Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=False, Rzi=False, Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=False, Rzj=False):
        """Defines member end realeses for a member. All member end releases will default to unreleased unless specified otherwise.

        :param member_name: The name of the member to have its releases modified.
        :type member_name: str
        :param Dxi: Indicates whether the member is released axially at its start. Defaults to False.
        :type Dxi: bool, optional
        :param Dyi: Indicates whether the member is released for shear in the local y-axis at its start. Defaults to False.
        :type Dyi: bool, optional
        :param Dzi: Indicates whether the member is released for shear in the local z-axis at its start. Defaults to False.
        :type Dzi: bool, optional
        :param Rxi: Indicates whether the member is released for torsion at its start. Defaults to False.
        :type Rxi: bool, optional
        :param Ryi: Indicates whether the member is released for moment about the local y-axis at its start. Defaults to False.
        :type Ryi: bool, optional
        :param Rzi: Indicates whether the member is released for moment about the local z-axis at its start. Defaults to False.
        :type Rzi: bool, optional
        :param Dxj: Indicates whether the member is released axially at its end. Defaults to False.
        :type Dxj: bool, optional
        :param Dyj: Indicates whether the member is released for shear in the local y-axis at its end. Defaults to False.
        :type Dyj: bool, optional
        :param Dzj: Indicates whether the member is released for shear in the local z-axis. Defaults to False.
        :type Dzj: bool, optional
        :param Rxj: Indicates whether the member is released for torsion at its end. Defaults to False.
        :type Rxj: bool, optional
        :param Ryj: Indicates whether the member is released for moment about the local y-axis at its end. Defaults to False.
        :type Ryj: bool, optional
        :param Rzj: Indicates whether the member is released for moment about the local z-axis at its end. Defaults to False.
        :type Rzj: bool, optional
        """
        
        # Apply the end releases to the member
        self.members[member_name].Releases = [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]     

        # Flag the model as unsolved
        self.solution = None

    def add_load_combo(self, name, factors, combo_tags=None):
        """Adds a load combination to the model.

        :param name: A unique name for the load combination (e.g. '1.2D+1.6L+0.5S' or 'Gravity Combo').
        :type name: str
        :param factors: A dictionary containing load cases and their corresponding factors (e.g. {'D':1.2, 'L':1.6, 'S':0.5}).
        :type factors: dict
        :param combo_tags: A list of tags used to categorize load combinations. Default is `None`. This can be useful for filtering results later on, or for limiting analysis to only those combinations with certain tags. This feature is provided for convenience. It is not necessary to use tags.
        :type combo_tags: list, optional
        """            

        # Create a new load combination object
        new_combo = LoadCombo(name, combo_tags, factors)

        # Add the load combination to the dictionary of load combinations
        self.load_combos[name] = new_combo

        # Flag the model as solved
        self.solution = None

    def add_node_load(self, node_name, direction, P, case='Case 1'):
        """Adds a nodal load to the model.

        :param node_name: The name of the node where the load is being applied.
        :type node_name: str
        :param direction: The global direction the load is being applied in. Forces are `'FX'`,
                          `'FY'`, and `'FZ'`. Moments are `'MX'`, `'MY'`, and `'MZ'`.
        :type direction: str
        :param P: The numeric value (magnitude) of the load.
        :type P: float
        :param case: The name of the load case the load belongs to. Defaults to 'Case 1'.
        :type case: str, optional
        :raises ValueError: Occurs when an invalid load direction was specified.
        """
        
        # Validate the value of direction
        if direction not in ('FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'):
            raise ValueError(f"direction must be 'FX', 'FY', 'FZ', 'MX', 'MY', or 'MZ'. {direction} was given.")
        # Add the node load to the model
        self.nodes[node_name].NodeLoads.append((direction, P, case))

        # Flag the model as unsolved
        self.solution = None

    def add_member_pt_load(self, member_name, direction, P, x, case='Case 1'):
        """Adds a member point load to the model.

        :param member_name: The name of the member the load is being applied to.
        :type member_name: str
        :param direction: The direction in which the load is to be applied. Valid values are `'Fx'`,
                          `'Fy'`, `'Fz'`, `'Mx'`, `'My'`, `'Mz'`, `'FX'`, `'FY'`, `'FZ'`, `'MX'`, `'MY'`, or `'MZ'`.
                          Note that lower-case notation indicates use of the beam's local
                          coordinate system, while upper-case indicates use of the model's globl
                          coordinate system.
        :type direction: str
        :param P: The numeric value (magnitude) of the load.
        :type P: float
        :param x: The load's location along the member's local x-axis.
        :type x: float
        :param case: The load case to categorize the load under. Defaults to 'Case 1'.
        :type case: str, optional
        :raises ValueError: Occurs when an invalid load direction has been specified.
        """            

        # Validate the value of direction
        if direction not in ('Fx', 'Fy', 'Fz', 'FX', 'FY', 'FZ', 'Mx', 'My', 'Mz', 'MX', 'MY', 'MZ'):
            raise ValueError(f"direction must be 'Fx', 'Fy', 'Fz', 'FX', 'FY', FZ', 'Mx', 'My', 'Mz', 'MX', 'MY', or 'MZ'. {direction} was given.")
        
        # Add the point load to the member
        self.members[member_name].PtLoads.append((direction, P, x, case))
        
        # Flag the model as unsolved
        self.solution = None

    def add_member_dist_load(self, member_name, direction, w1, w2, x1=None, x2=None, case='Case 1'):
        """Adds a member distributed load to the model.

        :param member_name: The name of the member the load is being appied to.
        :type member_name: str
        :param direction: The direction in which the load is to be applied. Valid values are `'Fx'`,
                          `'Fy'`, `'Fz'`, `'FX'`, `'FY'`, or `'FZ'`.
                          Note that lower-case notation indicates use of the beam's local
                          coordinate system, while upper-case indicates use of the model's globl
                          coordinate system.
        :type direction: str
        :param w1: The starting value (magnitude) of the load.
        :type w1: float
        :param w2: The ending value (magnitude) of the load.
        :type w2: float
        :param x1: The load's start location along the member's local x-axis. If this argument is
                   not specified, the start of the member will be used. Defaults to `None`
        :type x1: float, optional
        :param x2: The load's end location along the member's local x-axis. If this argument is not
                   specified, the end of the member will be used. Defaults to `None`.
        :type x2: float, optional
        :param case: _description_, defaults to 'Case 1'
        :type case: str, optional
        :raises ValueError: Occurs when an invalid load direction has been specified.
        """
       
        # Validate the value of direction
        if direction not in ('Fx', 'Fy', 'Fz', 'FX', 'FY', 'FZ'):
            raise ValueError(f"direction must be 'Fx', 'Fy', 'Fz', 'FX', 'FY', or 'FZ'. {direction} was given.")
        # Determine if a starting and ending points for the load have been specified.
        # If not, use the member start and end as defaults
        if x1 == None:
            start = 0
        else:
            start = x1
        
        if x2 == None:
            end = self.members[member_name].L()
        else:
            end = x2

        # Add the distributed load to the member
        self.members[member_name].DistLoads.append((direction, w1, w2, start, end, case))
        
        # Flag the model as unsolved
        self.solution = None

    def add_member_self_weight(self, global_direction, factor, case='Case 1'):
        """Adds self weight to all members in the model. Note that this only works for members. Plate and Quad elements will be ignored by this command.

        :param global_direction: The global direction to apply the member load in: 'FX', 'FY', or 'FZ'.
        :type global_direction: string
        :param factor: A factor to apply to the member self-weight. Can be used to account for items like connections, or to switch the direction of the self-weight load.
        :type factor: float
        :param case: The load case to apply the self-weight to. Defaults to 'Case 1'
        :type case: str, optional
        """

        # Step through each member in the model
        for member in self.members.values():

            # Calculate the self weight of the member
            self_weight = factor*self.materials[member.material_name].rho*member.A

            # Add the self-weight load to the member
            self.add_member_dist_load(member.name, global_direction, self_weight, self_weight, case=case)
        
        # No need to flag the model as unsolved. That has already been taken care of by our call to `add_member_dist_load`

    def add_plate_surface_pressure(self, plate_name, pressure, case='Case 1'):
        """Adds a surface pressure to the rectangular plate element.
        

        :param plate_name: The name for the rectangular plate to add the surface pressure to.
        :type plate_name: str
        :param pressure: The value (magnitude) for the surface pressure.
        :type pressure: float
        :param case: The load case to add the surface pressure to. Defaults to 'Case 1'.
        :type case: str, optional
        :raises Exception: Occurs when an invalid plate name has been specified.
        """   

        # Add the surface pressure to the rectangle
        if plate_name in self.plates.keys():
            self.plates[plate_name].pressures.append([pressure, case])
        else:
            raise Exception('Invalid plate name specified for plate surface pressure.')
        
        # Flag the model as unsolved
        self.solution = None

    def add_quad_surface_pressure(self, quad_name, pressure, case='Case 1'):
        """Adds a surface pressure to the quadrilateral element.

        :param quad_name: The name for the quad to add the surface pressure to.
        :type quad_name: str
        :param pressure: The value (magnitude) for the surface pressure.
        :type pressure: float
        :param case: The load case to add the surface pressure to. Defaults to 'Case 1'.
        :type case: str, optional
        :raises Exception: Occurs when an invalid quad name has been specified.
        """

        # Add the surface pressure to the quadrilateral
        if quad_name in self.quads.keys():
            self.quads[quad_name].pressures.append([pressure, case])
        else:
            raise Exception('Invalid quad name specified for quad surface pressure.')
        
        # Flag the model as unsolved
        self.solution = None

    def delete_loads(self):
        """Deletes all loads from the model along with any results based on the loads.
        """

        # Delete the member loads and the calculated internal forces
        for member in self.members.values():
            member.DistLoads = []
            member.PtLoads = []
            member.SegmentsZ = []
            member.SegmentsY = []
            member.SegmentsX = []
        
        # Delete the plate loads
        for plate in self.plates.values():
            plate.pressures = []
        
        # Delete the quadrilateral loads
        for quad in self.quads.values():
            quad.pressures = []
        
        # Delete the nodal loads, calculated displacements, and calculated reactions
        for node in self.nodes.values():

            node.NodeLoads = []

            node.DX = {}
            node.DY = {}
            node.DZ = {}
            node.RX = {}
            node.RY = {}
            node.RZ = {}

            node.RxnFX = {}
            node.RxnFY = {}
            node.RxnFZ = {}
            node.RxnMX = {}
            node.RxnMY = {}
            node.RxnMZ = {}
        
        # Flag the model as unsolved
        self.solution = None
               
    def K(self, combo_name='Combo 1', log=False, check_stability=True, sparse=True):
        """Returns the model's global stiffness matrix. The stiffness matrix will be returned in
           scipy's sparse lil format, which reduces memory usage and can be easily converted to
           other formats.

        :param combo_name: The load combination to get the stiffness matrix for. Defaults to 'Combo 1'.
        :type combo_name: str, optional
        :param log: Prints updates to the console if set to True. Defaults to False.
        :type log: bool, optional
        :param check_stability: Causes Pynite to check for instabilities if set to True. Defaults
                                to True. Set to False if you want the model to run faster.
        :type check_stability: bool, optional
        :param sparse: Returns a sparse matrix if set to True, and a dense matrix otherwise.
                       Defaults to True.
        :type sparse: bool, optional
        :return: The global stiffness matrix for the structure.
        :rtype: ndarray or coo_matrix
        """           
        
        # Determine if a sparse matrix has been requested
        if sparse == True:
            # The stiffness matrix will be stored as a scipy `coo_matrix`. Scipy's
            # documentation states that this type of matrix is ideal for efficient
            # construction of finite element matrices. When converted to another
            # format, the `coo_matrix` sums values at the same (i, j) index. We'll
            # build the matrix from three lists.
            row = []
            col = []
            data = []
        else:
            # Initialize a dense matrix of zeros
            K = zeros((len(self.nodes)*6, len(self.nodes)*6))
        
        # Add stiffness terms for each nodal spring in the model
        if log: print('- Adding nodal spring support stiffness terms to global stiffness matrix')
        for node in self.nodes.values():
            
            # Determine if the node has any spring supports
            if node.spring_DX[0] != None:

                # Check for an active spring support
                if node.spring_DX[2] == True:
                    m, n = node.ID*6, node.ID*6
                    if sparse == True:
                        data.append(float(node.spring_DX[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_DX[0])

            if node.spring_DY[0] != None:

                # Check for an active spring support
                if node.spring_DY[2] == True:
                    m, n = node.ID*6 + 1, node.ID*6 + 1
                    if sparse == True:
                        data.append(float(node.spring_DY[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_DY[0])

            if node.spring_DZ[0] != None:

                # Check for an active spring support
                if node.spring_DZ[2] == True:
                    m, n = node.ID*6 + 2, node.ID*6 + 2
                    if sparse == True:
                        data.append(float(node.spring_DZ[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_DZ[0])
        
            if node.spring_RX[0] != None:

                # Check for an active spring support
                if node.spring_RX[2] == True:
                    m, n = node.ID*6 + 3, node.ID*6 + 3
                    if sparse == True:
                        data.append(float(node.spring_RX[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_RX[0])

            if node.spring_RY[0] != None:

                # Check for an active spring support
                if node.spring_RY[2] == True:
                    m, n = node.ID*6 + 4, node.ID*6 + 4
                    if sparse == True:
                        data.append(float(node.spring_RY[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_RY[0])
            
            if node.spring_RZ[0] != None:

                # Check for an active spring support
                if node.spring_RZ[2] == True:
                    m, n = node.ID*6 + 5, node.ID*6 + 5
                    if sparse == True:
                        data.append(float(node.spring_RZ[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_RZ[0])

        # Add stiffness terms for each spring in the model
        if log: print('- Adding spring stiffness terms to global stiffness matrix')
        for spring in self.springs.values():
            
            if spring.active[combo_name] == True:

                # Get the spring's global stiffness matrix
                # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                spring_K = spring.K()

                # Step through each term in the spring's stiffness matrix
                # 'a' & 'b' below are row/column indices in the spring's stiffness matrix
                # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
                for a in range(12):
                
                    # Determine if index 'a' is related to the i-node or j-node
                    if a < 6:
                        # Find the corresponding index 'm' in the global stiffness matrix
                        m = spring.i_node.ID*6 + a
                    else:
                        # Find the corresponding index 'm' in the global stiffness matrix
                        m = spring.j_node.ID*6 + (a-6)
                    
                    for b in range(12):
                    
                        # Determine if index 'b' is related to the i-node or j-node
                        if b < 6:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = spring.i_node.ID*6 + b
                        else:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = spring.j_node.ID*6 + (b-6)
                    
                        # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                        if sparse == True:
                            row.append(m)
                            col.append(n)
                            data.append(spring_K[a, b])
                        else:
                            K[m, n] += spring_K[a, b]

        # Add stiffness terms for each physical member in the model
        if log: print('- Adding member stiffness terms to global stiffness matrix')
        for phys_member in self.members.values():
            
            # Check to see if the physical member is active for the given load combination
            if phys_member.active[combo_name] == True:

                # Step through each sub-member in the physical member and add terms
                for member in phys_member.sub_members.values():
                    
                    # Get the member's global stiffness matrix
                    # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                    member_K = member.K()

                    # Step through each term in the member's stiffness matrix
                    # 'a' & 'b' below are row/column indices in the member's stiffness matrix
                    # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
                    for a in range(12):
                    
                        # Determine if index 'a' is related to the i-node or j-node
                        if a < 6:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.i_node.ID*6 + a
                        else:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.j_node.ID*6 + (a-6)
                        
                        for b in range(12):
                        
                            # Determine if index 'b' is related to the i-node or j-node
                            if b < 6:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.i_node.ID*6 + b
                            else:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.j_node.ID*6 + (b-6)
                        
                            # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                            if sparse == True:
                                row.append(m)
                                col.append(n)
                                data.append(member_K[a, b])
                            else:
                                K[m, n] += member_K[a, b]
                
        # Add stiffness terms for each quadrilateral in the model
        if log: print('- Adding quadrilateral stiffness terms to global stiffness matrix')
        for quad in self.quads.values():
            
            # Get the quadrilateral's global stiffness matrix
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            quad_K = quad.K()

            # Step through each term in the quadrilateral's stiffness matrix
            # 'a' & 'b' below are row/column indices in the quadrilateral's stiffness matrix
            # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
            for a in range(24):

                # Determine which node the index 'a' is related to
                if a < 6:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.i_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.j_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.m_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.n_node.ID*6 + (a - 18)

                for b in range(24):

                    # Determine which node the index 'b' is related to
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.i_node.ID*6 + b
                    elif b < 12:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.j_node.ID*6 + (b - 6)
                    elif b < 18:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.m_node.ID*6 + (b - 12)
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.n_node.ID*6 + (b - 18)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                    if sparse == True:
                        row.append(m)
                        col.append(n)
                        data.append(quad_K[a, b])
                    else:
                        K[m, n] += quad_K[a, b]
        
        # Add stiffness terms for each plate in the model
        if log: print('- Adding plate stiffness terms to global stiffness matrix')
        for plate in self.plates.values():
            
            # Get the plate's global stiffness matrix
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            plate_K = plate.K()

            # Step through each term in the plate's stiffness matrix
            # 'a' & 'b' below are row/column indices in the plate's stiffness matrix
            # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
            for a in range(24):

                # Determine which node the index 'a' is related to
                if a < 6:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.i_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.j_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.m_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.n_node.ID*6 + (a - 18)

                for b in range(24):

                    # Determine which node the index 'b' is related to
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.i_node.ID*6 + b
                    elif b < 12:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.j_node.ID*6 + (b - 6)
                    elif b < 18:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.m_node.ID*6 + (b - 12)
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.n_node.ID*6 + (b - 18)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                    if sparse == True:
                        row.append(m)
                        col.append(n)
                        data.append(plate_K[a, b])
                    else:
                        K[m, n] += plate_K[a, b]

        if sparse:
            # The stiffness matrix will be stored as a scipy `coo_matrix`. Scipy's
            # documentation states that this type of matrix is ideal for efficient
            # construction of finite element matrices. When converted to another
            # format, the `coo_matrix` sums values at the same (i, j) index.
            from scipy.sparse import coo_matrix
            row = array(row)
            col = array(col)
            data = array(data)
            K = coo_matrix((data, (row, col)), shape=(len(self.nodes)*6, len(self.nodes)*6))

        # Check that there are no nodal instabilities
        if check_stability:
            if log: print('- Checking nodal stability')
            if sparse: Analysis._check_stability(self, K.tocsr())
            else: Analysis._check_stability(self, K)

        # Return the global stiffness matrix
        return K    
   
    def Kg(self, combo_name='Combo 1', log=False, sparse=True, first_step=True):
        """Returns the model's global geometric stiffness matrix. Geometric stiffness of plates is not considered.

        :param combo_name: The name of the load combination to derive the matrix for. Defaults to 'Combo 1'.
        :type combo_name: str, optional
        :param log: Prints updates to the console if set to `True`. Defaults to `False`.
        :type log: bool, optional
        :param sparse: Returns a sparse matrix if set to `True`, and a dense matrix otherwise. Defaults to `True`.
        :type sparse: bool, optional
        :param first_step: Used to indicate if the analysis is occuring at the first load step. Used in nonlinear analysis where the load is broken into multiple steps. Default is `True`.
        :type first_step: book, optional
        :return: The global geometric stiffness matrix for the structure.
        :rtype: ndarray or coo_matrix
        """
        
        if sparse == True:
            # Initialize a zero matrix to hold all the stiffness terms. The matrix will be stored as a scipy sparse `lil_matrix`. This matrix format has several advantages. It uses less memory if the matrix is sparse, supports slicing, and can be converted to other formats (sparse or dense) later on for mathematical operations.
            from scipy.sparse import lil_matrix
            Kg = lil_matrix((len(self.nodes)*6, len(self.nodes)*6))
        else:
            Kg = zeros(len(self.nodes)*6, len(self.nodes)*6)
        
        # Add stiffness terms for each physical member in the model
        if log: print('- Adding member geometric stiffness terms to global geometric stiffness matrix')
        for phys_member in self.members.values():
            
            # Check to see if the physical member is active for the given load combination
            if phys_member.active[combo_name] == True:

                # Step through each sub-member in the physical member and add terms
                for member in phys_member.sub_members.values():

                    # Calculate the axial force in the member
                    E = member.E
                    A = member.A
                    L = member.L()

                    # Calculate the axial force acting on the member
                    if first_step:
                        # For the first load step take P = 0
                        P = 0
                    else:
                        # Calculate the member axial force due to axial strain
                        d = member.d(combo_name)
                        P = E*A/L*(d[6, 0] - d[0, 0])

                    # Get the member's global stiffness matrix
                    # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                    member_Kg = member.Kg(P)

                    # Step through each term in the member's stiffness matrix
                    # 'a' & 'b' below are row/column indices in the member's stiffness matrix
                    # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
                    for a in range(12):
                    
                        # Determine if index 'a' is related to the i-node or j-node
                        if a < 6:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.i_node.ID*6 + a
                        else:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.j_node.ID*6 + (a-6)
                        
                        for b in range(12):
                        
                            # Determine if index 'b' is related to the i-node or j-node
                            if b < 6:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.i_node.ID*6 + b
                            else:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.j_node.ID*6 + (b-6)
                        
                            # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                            Kg[m, n] += member_Kg[(a, b)]

        # Return the global geometric stiffness matrix
        return Kg
      
    def Km(self, combo_name='Combo 1', push_combo='Push', step_num=1, log=False, sparse=True):
        """Calculates the structure's global plastic reduction matrix, which is used for nonlinear inelastic analysis.

        :param combo_name: The name of the load combination to get the plastic reduction matrix for. Defaults to 'Combo 1'.
        :type combo_name: str, optional
        :param push_combo: The name of the load combination that contains the pushover load definition. Defaults to 'Push'.
        :type push_combo: str, optional
        :param step_num: The load step used to generate the plastic reduction matrix. Defaults to 1.
        :type step_num: int, optional
        :param log: Determines whether this method writes output to the console as it runs. Defaults to False.
        :type log: bool, optional
        :param sparse: Indicates whether the sparse solver should be used. Defaults to True.
        :type sparse: bool, optional
        :return: The gloabl plastic reduction matrix.
        :rtype: array
        """
       
        # Determine if a sparse matrix has been requested
        if sparse == True:
            # The plastic reduction matrix will be stored as a scipy `coo_matrix`. Scipy's documentation states that this type of matrix is ideal for efficient construction of finite element matrices. When converted to another format, the `coo_matrix` sums values at the same (i, j) index. We'll build the matrix from three lists.
            row = []
            col = []
            data = []
        else:
            # Initialize a dense matrix of zeros
            Km = zeros((len(self.nodes)*6, len(self.nodes)*6))

        # Add stiffness terms for each physical member in the model
        if log: print('- Calculating the plastic reduction matrix')
        for phys_member in self.members.values():
            
            # Check to see if the physical member is active for the given load combination
            if phys_member.active[combo_name] == True:

                # Step through each sub-member in the physical member and add terms
                for member in phys_member.sub_members.values():
                    
                    # Get the member's global plastic reduction matrix
                    # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                    member_Km = member.Km(combo_name, push_combo, step_num)

                    # Step through each term in the member's plastic reduction matrix
                    # 'a' & 'b' below are row/column indices in the member's matrix
                    # 'm' & 'n' are corresponding row/column indices in the structure's global matrix
                    for a in range(12):
                    
                        # Determine if index 'a' is related to the i-node or j-node
                        if a < 6:
                            # Find the corresponding index 'm' in the global plastic reduction matrix
                            m = member.i_node.ID*6 + a
                        else:
                            # Find the corresponding index 'm' in the global plastic reduction matrix
                            m = member.j_node.ID*6 + (a-6)
                        
                        for b in range(12):
                        
                            # Determine if index 'b' is related to the i-node or j-node
                            if b < 6:
                                # Find the corresponding index 'n' in the global plastic reduction matrix
                                n = member.i_node.ID*6 + b
                            else:
                                # Find the corresponding index 'n' in the global plastic reduction matrix
                                n = member.j_node.ID*6 + (b-6)
                        
                            # Now that 'm' and 'n' are known, place the term in the global plastic reduction matrix
                            if sparse == True:
                                row.append(m)
                                col.append(n)
                                data.append(member_Km[a, b])
                            else:
                                Km[m, n] += member_Km[a, b]

        if sparse:
            # The plastic reduction matrix will be stored as a scipy `coo_matrix`. Scipy's documentation states that this type of matrix is ideal for efficient construction of finite element matrices. When converted to another format, the `coo_matrix` sums values at the same (i, j) index.
            from scipy.sparse import coo_matrix
            row = array(row)
            col = array(col)
            data = array(data)
            Km = coo_matrix((data, (row, col)), shape=(len(self.nodes)*6, len(self.nodes)*6))

        # Check that there are no nodal instabilities
        # if check_stability:
        #     if log: print('- Checking nodal stability')
        #     if sparse: Analysis._check_stability(self, Km.tocsr())
        #     else: Analysis._check_stability(self, Km)

        # Return the global plastic reduction matrix
        return Km 

    def FER(self, combo_name='Combo 1'):
        """Assembles and returns the global fixed end reaction vector for any given load combo.

        :param combo_name: The name of the load combination to get the fixed end reaction vector
                           for. Defaults to 'Combo 1'.
        :type combo_name: str, optional
        :return: The fixed end reaction vector
        :rtype: array
        """        
        
        # Initialize a zero vector to hold all the terms
        FER = zeros((len(self.nodes) * 6, 1))
        
        # Step through each physical member in the model
        for phys_member in self.members.values():
            
            # Step through each sub-member and add terms
            for member in phys_member.sub_members.values():

                # Get the member's global fixed end reaction vector
                # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                member_FER = member.FER(combo_name)

                # Step through each term in the member's fixed end reaction vector
                # 'a' below is the row index in the member's fixed end reaction vector
                # 'm' below is the corresponding row index in the global fixed end reaction vector
                for a in range(12):
                    
                    # Determine if index 'a' is related to the i-node or j-node
                    if a < 6:
                        # Find the corresponding index 'm' in the global fixed end reaction vector
                        m = member.i_node.ID * 6 + a
                    else:
                        # Find the corresponding index 'm' in the global fixed end reaction vector
                        m = member.j_node.ID * 6 + (a - 6)
                    
                    # Now that 'm' is known, place the term in the global fixed end reaction vector
                    FER[m, 0] += member_FER[a, 0]
        
        # Add terms for each rectangle in the model
        for plate in self.plates.values():
            
            # Get the quadrilateral's global fixed end reaction vector
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            plate_FER = plate.FER(combo_name)

            # Step through each term in the quadrilateral's fixed end reaction vector
            # 'a' below is the row index in the quadrilateral's fixed end reaction vector
            # 'm' below is the corresponding row index in the global fixed end reaction vector
            for a in range(24):
                
                # Determine if index 'a' is related to the i-node, j-node, m-node, or n-node
                if a < 6:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.i_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.j_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.m_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.n_node.ID*6 + (a - 18)
                
                # Now that 'm' is known, place the term in the global fixed end reaction vector
                FER[m, 0] += plate_FER[a, 0]

        # Add terms for each quadrilateral in the model
        for quad in self.quads.values():
            
            # Get the quadrilateral's global fixed end reaction vector
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            quad_FER = quad.FER(combo_name)

            # Step through each term in the quadrilateral's fixed end reaction vector
            # 'a' below is the row index in the quadrilateral's fixed end reaction vector
            # 'm' below is the corresponding row index in the global fixed end reaction vector
            for a in range(24):
                
                # Determine if index 'a' is related to the i-node, j-node, m-node, or n-node
                if a < 6:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.i_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.j_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.m_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.n_node.ID*6 + (a - 18)
                
                # Now that 'm' is known, place the term in the global fixed end reaction vector
                FER[m, 0] += quad_FER[a, 0]

        # Return the global fixed end reaction vector
        return FER
    
    def P(self, combo_name='Combo 1'):
        """Assembles and returns the global nodal force vector.

        :param combo_name: The name of the load combination to get the force vector for. Defaults
                           to 'Combo 1'.
        :type combo_name: str, optional
        :return: The global nodal force vector.
        :rtype: array
        """
            
        # Initialize a zero vector to hold all the terms
        P = zeros((len(self.nodes)*6, 1))
        
        # Get the load combination for the given 'combo_name'
        combo = self.load_combos[combo_name]

        # Add terms for each node in the model
        for node in self.nodes.values():
            
            # Get the node's ID
            ID = node.ID

            # Step through each load factor in the load combination
            for case, factor in combo.factors.items():

                # Add the node's loads to the global nodal load vector
                for load in node.NodeLoads:

                    if load[2] == case:

                        if load[0] == 'FX':
                            P[ID*6 + 0, 0] += factor*load[1]
                        elif load[0] == 'FY':
                            P[ID*6 + 1, 0] += factor*load[1]
                        elif load[0] == 'FZ':
                            P[ID*6 + 2, 0] += factor*load[1]
                        elif load[0] == 'MX':
                            P[ID*6 + 3, 0] += factor*load[1]
                        elif load[0] == 'MY':
                            P[ID*6 + 4, 0] += factor*load[1]
                        elif load[0] == 'MZ':
                            P[ID*6 + 5, 0] += factor*load[1]
        
        # Return the global nodal force vector
        return P

    def D(self, combo_name='Combo 1'):
        """Returns the global displacement vector for the model.

        :param combo_name: The name of the load combination to get the results for. Defaults to
                           'Combo 1'.
        :type combo_name: str, optional
        :return: The global displacement vector for the model
        :rtype: array
        """
 
        # Return the global displacement vector
        return self._D[combo_name]

    def analyze(self, log=False, check_stability=True, check_statics=False, max_iter=30, sparse=True, combo_tags=None):
        """Performs first-order static analysis. Iterations are performed if tension-only members or compression-only members are present.

        :param log: Prints the analysis log to the console if set to True. Default is False.
        :type log: bool, optional
        :param check_stability: When set to `True`, checks for nodal instabilities. This slows down analysis a little. Default is `True`.
        :type check_stability: bool, optional
        :param check_statics: When set to `True`, causes a statics check to be performed
        :type check_statics: bool, optional
        :param max_iter: The maximum number of iterations to try to get convergence for tension/compression-only analysis. Defaults to 30.
        :type max_iter: int, optional
        :param sparse: Indicates whether the sparse matrix solver should be used. A matrix can be considered sparse or dense depening on how many zero terms there are. Structural stiffness matrices often contain many zero terms. The sparse solver can offer faster solutions for such matrices. Using the sparse solver on dense matrices may lead to slower solution times.
        :type sparse: bool, optional
        :raises Exception: _description_
        :raises Exception: _description_
        """

        if log:
            print('+-----------+')
            print('| Analyzing |')
            print('+-----------+')

        # Import `scipy` features if the sparse solver is being used
        if sparse == True:
            from PyNite.Solvers import spsolve

        # Prepare the model for analysis
        Analysis._prepare_model(self)

        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = Analysis._partition_D(self)

        # Identify which load combinations have the tags the user has given
        combo_list = Analysis._identify_combos(self, combo_tags)

        # Step through each load combination
        for combo in combo_list:

            if log:
                print('')
                print('- Analyzing load combination ' + combo.name)

            # Keep track of the number of iterations
            iter_count = 1
            convergence = False
            divergence = False

            # Iterate until convergence or divergence occurs
            while convergence == False and divergence == False:
                
                # Check for tension/compression-only divergence
                if iter_count > max_iter:
                    divergence = True
                    raise Exception('Model diverged during tension/compression-only analysis')
                
                # Get the partitioned global stiffness matrix K11, K12, K21, K22
                if sparse == True:
                    K11, K12, K21, K22 = Analysis._partition(self, self.K(combo.name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)
                else:
                    K11, K12, K21, K22 = Analysis._partition(self, self.K(combo.name, log, check_stability, sparse), D1_indices, D2_indices)

                # Get the partitioned global fixed end reaction vector
                FER1, FER2 = Analysis._partition(self, self.FER(combo.name), D1_indices, D2_indices)

                # Get the partitioned global nodal force vector       
                P1, P2 = Analysis._partition(self, self.P(combo.name), D1_indices, D2_indices)          

                # Calculate the global displacement vector
                if log: print('- Calculating global displacement vector')
                if K11.shape == (0, 0):
                    # All displacements are known, so D1 is an empty vector
                    D1 = []
                else:
                    try:
                        # Calculate the unknown displacements D1
                        if sparse == True:
                            # The partitioned stiffness matrix is in `lil` format, which is great
                            # for memory, but slow for mathematical operations. The stiffness
                            # matrix will be converted to `csr` format for mathematical operations.
                            # The `@` operator performs matrix multiplication on sparse matrices.
                            D1 = spsolve(K11.tocsr(), subtract(subtract(P1, FER1), K12.tocsr() @ D2))
                            D1 = D1.reshape(len(D1), 1)
                        else:
                            D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))
                    except:
                        # Return out of the method if 'K' is singular and provide an error message
                        raise Exception('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')

                # Store the calculated displacements to the model and the nodes in the model
                Analysis._store_displacements(self, D1, D2, D1_indices, D2_indices, combo)
                
                # Check for tension/compression-only convergence
                convergence = Analysis._check_TC_convergence(self, combo.name, log=log)

                if convergence == False:
                    if log: print('- Tension/compression-only analysis did not converge. Adjusting stiffness matrix and reanalyzing.')
                else:
                    if log: print('- Tension/compression-only analysis converged after ' + str(iter_count) + ' iteration(s)')

                # Keep track of the number of tension/compression only iterations
                iter_count += 1

        # Calculate reactions
        Analysis._calc_reactions(self, log, combo_tags)

        if log:
            print('')     
            print('- Analysis complete')
            print('')

        # Check statics if requested
        if check_statics == True:
            Analysis._check_statics(self, combo_tags)
        
        # Flag the model as solved
        self.solution = 'Linear TC'

    def analyze_linear(self, log=False, check_stability=True, check_statics=False, sparse=True, combo_tags=None):
        """Performs first-order static analysis. This analysis procedure is much faster since it only assembles the global stiffness matrix once, rather than once for each load combination. It is not appropriate when non-linear behavior such as tension/compression only analysis or P-Delta analysis are required.

        :param log: Prints the analysis log to the console if set to True. Default is False.
        :type log: bool, optional
        :param check_stability: When set to True, checks the stiffness matrix for any unstable degrees of freedom and reports them back to the console. This does add to the solution time. Defaults to True.
        :type check_stability: bool, optional
        :param check_statics: When set to True, causes a statics check to be performed. Defaults to False.
        :type check_statics: bool, optional
        :param sparse: Indicates whether the sparse matrix solver should be used. A matrix can be considered sparse or dense depening on how many zero terms there are. Structural stiffness matrices often contain many zero terms. The sparse solver can offer faster solutions for such matrices. Using the sparse solver on dense matrices may lead to slower solution times. Be sure ``scipy`` is installed to use the sparse solver. Default is True.
        :type sparse: bool, optional
        :raises Exception: Occurs when a singular stiffness matrix is found. This indicates an unstable structure has been modeled.
        """

        if log:
            print('+-------------------+')
            print('| Analyzing: Linear |')
            print('+-------------------+')
        
        # Import `scipy` features if the sparse solver is being used
        if sparse == True:
            from scipy.sparse.linalg import spsolve

        # Prepare the model for analysis
        Analysis._prepare_model(self)

        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = Analysis._partition_D(self)

        # Get the partitioned global stiffness matrix K11, K12, K21, K22
        # Note that for linear analysis the stiffness matrix can be obtained for any load combination, as it's the same for all of them
        combo_name = list(self.load_combos.keys())[0]
        if sparse == True:
            K11, K12, K21, K22 = Analysis._partition(self, self.K(combo_name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)
        else:
            K11, K12, K21, K22 = Analysis._partition(self, self.K(combo_name, log, check_stability, sparse), D1_indices, D2_indices)

        # Identify which load combinations have the tags the user has given
        combo_list = Analysis._identify_combos(self, combo_tags)

        # Step through each load combination
        for combo in combo_list:

            if log:
                print('')
                print('- Analyzing load combination ' + combo.name)

            # Get the partitioned global fixed end reaction vector
            FER1, FER2 = Analysis._partition(self, self.FER(combo.name), D1_indices, D2_indices)

            # Get the partitioned global nodal force vector       
            P1, P2 = Analysis._partition(self, self.P(combo.name), D1_indices, D2_indices)          

            # Calculate the global displacement vector
            if log: print('- Calculating global displacement vector')
            if K11.shape == (0, 0):
                # All displacements are known, so D1 is an empty vector
                D1 = []
            else:
                try:
                    # Calculate the unknown displacements D1
                    if sparse == True:
                        # The partitioned stiffness matrix is in `lil` format, which is great
                        # for memory, but slow for mathematical operations. The stiffness
                        # matrix will be converted to `csr` format for mathematical operations.
                        # The `@` operator performs matrix multiplication on sparse matrices.
                        D1 = spsolve(K11.tocsr(), subtract(subtract(P1, FER1), K12.tocsr() @ D2))
                        D1 = D1.reshape(len(D1), 1)
                    else:
                        D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))
                except:
                    # Return out of the method if 'K' is singular and provide an error message
                    raise Exception('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')

            # Store the calculated displacements to the model and the nodes in the model
            Analysis._store_displacements(self, D1, D2, D1_indices, D2_indices, combo)

        # Calculate reactions
        Analysis._calc_reactions(self, log, combo_tags)

        if log:
            print('')     
            print('- Analysis complete')
            print('')

        # Check statics if requested
        if check_statics == True:
            Analysis._check_statics(self, combo_tags)
        
        # Flag the model as solved
        self.solution = 'Linear'

    def analyze_PDelta(self, log=False, check_stability=True, max_iter=30, sparse=True, combo_tags=None):
        """Performs second order (P-Delta) analysis. This type of analysis is appropriate for most models using beams, columns and braces. Second order analysis is usually required by material specific codes. The analysis is iterative and takes longer to solve. Models with slender members and/or members with combined bending and axial loads will generally have more significant P-Delta effects. P-Delta effects in plates/quads are not considered.

        :param log: Prints updates to the console if set to True. Default is False.
        :type log: bool, optional
        :param check_stability: When set to True, checks the stiffness matrix for any unstable degrees of freedom and reports them back to the console. This does add to the solution time. Defaults to True.
        :type check_stability: bool, optional
        :param max_iter: The maximum number of iterations permitted. If this value is exceeded the program will report divergence. Defaults to 30.
        :type max_iter: int, optional
        :param sparse: Indicates whether the sparse matrix solver should be used. A matrix can be considered sparse or dense depening on how many zero terms there are. Structural stiffness matrices often contain many zero terms. The sparse solver can offer faster solutions for such matrices. Using the sparse solver on dense matrices may lead to slower solution times. Be sure ``scipy`` is installed to use the sparse solver. Default is True.
        :type sparse: bool, optional
        :raises ValueError: Occurs when there is a singularity in the stiffness matrix, which indicates an unstable structure.
        :raises Exception: Occurs when a model fails to converge.
        """
        
        if log:
            print('+--------------------+')
            print('| Analyzing: P-Delta |')
            print('+--------------------+')

        # Import `scipy` features if the sparse solver is being used
        if sparse == True:
            from scipy.sparse.linalg import spsolve

        # Prepare the model for analysis
        Analysis._prepare_model(self)
        
        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = Analysis._partition_D(self)

        # Identify which load combinations have the tags the user has given
        combo_list = Analysis._identify_combos(self, combo_tags)

        # Step through each load combination
        for combo in combo_list:

            # Get the partitioned global fixed end reaction vector
            FER1, FER2 = Analysis._partition(self, self.FER(combo.name), D1_indices, D2_indices)

            # Get the partitioned global nodal force vector       
            P1, P2 = Analysis._partition(self, self.P(combo.name), D1_indices, D2_indices)

            # Run the P-Delta analysis for this load combination
            Analysis._PDelta_step(self, combo.name, P1, FER1, D1_indices, D2_indices, D2, log, sparse, check_stability, max_iter, True)
        
        # Calculate reactions
        Analysis._calc_reactions(self, log, combo_tags)

        if log:
            print('')
            print('- Analysis complete')
            print('')
        
        # Flag the model as solved
        self.solution = 'P-Delta'
    
    def _not_ready_yet_analyze_pushover(self, log=False, check_stability=True, push_combo='Push', max_iter=30, tol=0.01, sparse=True, combo_tags=None):

        if log:
            print('+---------------------+')
            print('| Analyzing: Pushover |')
            print('+---------------------+')
        
        # Import `scipy` features if the sparse solver is being used
        if sparse == True:
            from scipy.sparse.linalg import spsolve

        # Prepare the model for analysis
        Analysis._prepare_model(self)
        
        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = Analysis._partition_D(self)

        # Identify and tag the primary load combinations the pushover load will be added to
        for combo in self.load_combos.values():

            # No need to tag the pushover combo
            if combo.name != push_combo:

                # Add 'primary' to the combo's tags if it's not already there
                if combo.combo_tags is None:
                    combo.combo_tags = ['primary']
                elif 'primary' not in combo.combo_tags:
                    combo.combo_tags.append('primary')

        # Identify which load combinations have the tags the user has given
        # TODO: Remove the pushover combo istelf from `combo_list`
        combo_list = Analysis._identify_combos(self, combo_tags)
        combo_list = [combo for combo in combo_list if combo.name != push_combo]

        # Step through each load combination
        for combo in combo_list:
            
            # Skip the pushover combo
            if combo.name == push_combo:
                continue

            if log:
                print('')
                print('- Analyzing load combination ' + combo.name)
            
            # Reset nonlinear material member end forces to zero
            for member in self.members.values():
                member._fxi, member._myi, member._mzi = 0, 0, 0
                member._fxj, member._myj, member._mzj = 0, 0, 0

            # Get the pushover load step and initialize the load factor
            load_step = list(self.load_combos[push_combo].factors.values())[0]
            load_factor = load_step
            step_num = 1
            
            # Get the partitioned global fixed end reaction vector
            FER1, FER2 = Analysis._partition(self, self.FER(combo.name), D1_indices, D2_indices)

            # Get the partitioned global nodal force vector       
            P1, P2 = Analysis._partition(self, self.P(combo.name), D1_indices, D2_indices)

            # Get the partitioned global fixed end reaction vector for a pushover load increment
            FER1_push, FER2_push = Analysis._partition(self, self.FER(push_combo), D1_indices, D2_indices)

            # Get the partitioned global nodal force vector for a pushover load increment
            P1_push, P2_push = Analysis._partition(self, self.P(push_combo), D1_indices, D2_indices)

            # Solve the current load combination without the pushover load applied
            Analysis._PDelta_step(self, combo.name, P1, FER1, D1_indices, D2_indices, D2, log, sparse, check_stability, max_iter, first_step=True)

            # Since a P-Delta analysis was just run, we'll need to correct the solution to flag it
            # as 'pushover' instead of 'PDelta'
            self.solution = 'Pushover'

            # Apply the pushover load in steps, summing deformations as we go, until the full
            # pushover load has been analyzed
            while load_factor <= 1:
                
                # Inform the user which pushover load step we're on
                if log:
                    print('- Beginning pushover load step #' + str(step_num))
                
                # Reset all member plastic load reversal flags to be `False`
                for member in self.members.values():
                    member.pl_reverse = False

                # Run/rerun this load step until no new unloaded member flags exist
                run_step = True
                while run_step == True:

                    # Assume this iteration will converge
                    run_step = False

                    # Store the model's current displacements in case we need to revert
                    D_temp = self._D

                    # Run or rerun the next pushover load step
                    d_Delta = Analysis._pushover_step(self, combo.name, push_combo, step_num, P1_push, FER1_push, D1_indices, D2_indices, D2, log, sparse, check_stability)

                # Update nonlinear material member end forces for each member
                for member in self.members.values():
                    member._fxi = member.f(combo.name, push_combo, step_num)[0, 0]
                    member._myi = member.f(combo.name, push_combo, step_num)[4, 0]
                    member._mzi = member.f(combo.name, push_combo, step_num)[5, 0]
                    member._fxj = member.f(combo.name, push_combo, step_num)[6, 0]
                    member._myj = member.f(combo.name, push_combo, step_num)[10, 0]
                    member._mzj = member.f(combo.name, push_combo, step_num)[11, 0]

                # Move on to the next load step
                load_factor += load_step
                step_num += 1

        # Calculate reactions for every primary load combination
        Analysis._calc_reactions(self, log, combo_tags=['primary'])

        if log:
            print('')
            print('- Analysis complete')
            print('')
        
        # Flag the model as solved
        self.solution = 'Pushover'
    
    def unique_name(self, dictionary, prefix):
        """Returns the next available unique name for a dictionary of objects.

        :param dictionary: The dictionary to get a unique name for.
        :type dictionary: dict
        :param prefix: The prefix to use for the unique name.
        :type prefix: str
        :return: A unique name for the dictionary.
        :rtype: str
        """

        # Select a trial value for the next available name
        name = prefix + str(len(dictionary) + 1)
        i = 2
        while name in dictionary.keys():
            name = prefix + str(len(dictionary) + i)
            i += 1
        
        # Return the next available name
        return name
    
    def rename(self):
        """
        Renames all the nodes and elements in the model.
        """

        # Rename each node in the model
        temp = self.nodes.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'N' + str(id)
            self.nodes[new_key] = self.nodes.pop(old_key)
            self.nodes[new_key].name = new_key
            id += 1
        
        # Rename each spring in the model
        temp = self.springs.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'S' + str(id)
            self.springs[new_key] = self.springs.pop(old_key)
            self.springs[new_key].name = new_key
            id += 1

        # Rename each member in the model
        temp = self.members.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'M' + str(id)
            self.members[new_key] = self.members.pop(old_key)
            self.members[new_key].name = new_key
            id += 1
        
        # Rename each plate in the model
        temp = self.plates.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'P' + str(id)
            self.plates[new_key] = self.plates.pop(old_key)
            self.plates[new_key].name = new_key
            id += 1
        
        # Rename each quad in the model
        temp = self.quads.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'Q' + str(id)
            self.quads[new_key] = self.quads.pop(old_key)
            self.quads[new_key].name = new_key
            id += 1

    def orphaned_nodes(self):
        """
        Returns a list of the names of nodes that are not attached to any elements.
        """

        # Initialize a list of orphaned nodes
        orphans = []

        # Step through each node in the model
        for node in self.nodes.values():

            orphaned = False

            # Check to see if the node is attached to any elements
            quads = [quad.name for quad in self.quads.values() if quad.i_node == node or quad.j_node == node or quad.m_node == node or quad.n_node == node]
            plates = [plate.name for plate in self.plates.values() if plate.i_node == node or plate.j_node == node or plate.m_node == node or plate.n_node == node]
            members = [member.name for member in self.members.values() if member.i_node == node or member.j_node == node]
            springs = [spring.name for spring in self.springs.values() if spring.i_node == node or spring.j_node == node]

            # Determine if the node is orphaned
            if quads == [] and plates == [] and members == [] and springs == []:
                orphaned = True
            
            # Add the orphaned nodes to the list of orphaned nodes
            if orphaned == True:
                orphans.append(node.name)
        
        return orphans
      
