#%%
# This is an example of a how plates can be used to create a wall panel for out-of-plane bending problems.
# A 'RectWall' class has been defined to automate the process of meshing, loading, and plotting results for a wall
# panel of any geometry and edge support conditions. Below the class definition is a brief script showing how the
# 'RectWall' class can be been implemented.

#%%
from PyNite import FEModel3D, Node3D, Plate3D

# Packages used for plotting contours
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

#%%
class RectWall():

    # Constructor
    def __init__(self, width, height, thickness, E=3824, nu=0.3, mesh_size=1, bot_support='Pinned', top_support='Pinned', left_support='Free', right_support='Free'):
        """
        Initialize a mesh.

        Args:
            self: (todo): write your description
            width: (int): write your description
            height: (int): write your description
            thickness: (todo): write your description
            E: (int): write your description
            nu: (int): write your description
            mesh_size: (int): write your description
            bot_support: (int): write your description
            top_support: (int): write your description
            left_support: (int): write your description
            right_support: (int): write your description
        """

        self.width = width
        self.height = height
        self.thickness = thickness
        self.mesh_size = mesh_size

        self.E = E
        self.nu = nu

        self.fem = FEModel3D()  # A finite element model for the wall

        self.loads = []   # A list of surface loads applied to the wall panel
        self.nodes = []   # A list of nodes that make up the wall panel
        self.plates = []  # A list of plates that make up the wall panel

        self.bot_support = bot_support
        self.top_support = top_support
        self.left_support = left_support
        self.right_support = right_support

        self.__analyzed = False
    
    # Adds a pressure load to the wall
    def add_load(self, y_bot, y_top, p_bot, p_top):
        """
        Add a load to the topology.

        Args:
            self: (todo): write your description
            y_bot: (str): write your description
            y_top: (str): write your description
            p_bot: (str): write your description
            p_top: (str): write your description
        """
        
        # Add the load to the wall
        self.loads.append((y_bot, y_top, p_bot, p_top))
        self.__analyzed = False
    
    # Removes all loads from the wall
    def clear_loads(self):
        """
        Clears all files from the database.

        Args:
            self: (todo): write your description
        """

        # Clear all the wall loads
        self.loads = []
        self.__analyzed = False
        self.fem.ClearLoads()
    
    # Descritizes the wall
    def __descritize(self):
        """
        Descriting description.

        Args:
            self: (todo): write your description
        """

        # Determine how many columns of plates are necessary
        num_cols = round(self.width/self.mesh_size)

        # Determine the plate width
        pl_width = self.width/num_cols

        # Create a list to store mesh control points along the height of the wall
        # and add the top and bottom of the wall to the list
        y = [0.0, self.height]

        # Add a control point for each height where the load changes
        for load_item in self.loads:

            # Add the load start and end locations
            y.append(load_item[0])
            y.append(load_item[1])

            # Remove duplicates from the list of control points
            y = list(set(y))
            
            # Sort the list
            y.sort()

        # Prepare to iterate through each control point after the first one
        iter_y = iter(y)
        next(iter_y)
        y_prev = 0.0

        # Initialize the plate and node ID's to 1
        pl_id = 1
        node_id = 1

        # Create the bottom row of nodes
        for i in range(num_cols + 1):
            self.fem.AddNode('N'+str(node_id), i*pl_width, 0.0, 0.0)
            node_id += 1

        # Add plates and nodes between each control point
        for control_point in iter_y:

            # Determine the height between control points ("lift" height)
            h = control_point - y_prev

            # Determine how many rows of plates are needed in the lift
            num_rows = max(1, round(h/self.mesh_size))

            # Determine the plate height
            pl_height = h/num_rows

            # Generate nodes
            for j in range(num_rows + 1):
                for i in range(num_cols + 1):
                    
                    # The first row of nodes in the lift has already been generated
                    if j != 0:

                        # Generate the node
                        self.fem.AddNode('N'+str(node_id), i*pl_width, j*pl_height + y_prev, 0.0)
                        node_id += 1
            
            # Generate plates
            for j in range(num_rows):
                for i in range(num_cols):

                    # Determine which nodes the plate will be attached to
                    ni = pl_id + max(0, int((pl_id-1)/num_cols))
                    nj = ni + 1
                    nm = nj + (num_cols + 1)
                    nn = nm - 1

                    # Add the plate to the list of plates
                    self.fem.AddPlate('P'+str(pl_id), 'N'+str(ni), 'N'+str(nj), 'N'+str(nm), 'N'+str(nn), self.thickness, self.E, self.nu) 

                    # Prepare to move to the next iteration
                    pl_id += 1
            
            # Prepare for the next iteration
            y_prev = control_point

    # Defines support conditions at each node
    def __define_supports(self):
        """
        Check whether or notclose is true if there.

        Args:
            self: (todo): write your description
        """

        # Step through each node in the model
        for node in self.fem.Nodes:

            # Determine if the node falls on any of the boundaries
            # Left edge
            if np.isclose(node.X, 0.0):
                if self.left_support == "Fixed":
                    node.SupportDX, node.SupportDY, node.SupportDZ, node.SupportRX, node.SupportRY = True, True, True, True, True
                elif self.left_support == "Pinned":
                    node.SupportDX, node.SupportDY, node.SupportDZ = True, True, True
            # Right edge
            elif np.isclose(node.X, self.width):
                if self.right_support == "Fixed":
                    node.SupportDX, node.SupportDY, node.SupportDZ, node.SupportRX, node.SupportRY = True, True, True, True, True
                elif self.right_support == "Pinned":
                    node.SupportDX, node.SupportDY, node.SupportDZ = True, True, True
            # Bottom edge
            elif np.isclose(node.Y, 0.0):
                if self.bot_support == "Fixed":
                    node.SupportDX, node.SupportDY, node.SupportDZ, node.SupportRX, node.SupportRY = True, True, True, True, True
                elif self.bot_support == "Pinned":
                    node.SupportDX, node.SupportDY, node.SupportDZ = True, True, True
            # Top edge
            elif np.isclose(node.Y, self.height):
                if self.top_support == "Fixed":
                    node.SupportDX, node.SupportDY, node.SupportDZ, node.SupportRX, node.SupportRY = True, True, True, True, True
                elif self.top_support == "Pinned":
                    node.SupportDX, node.SupportDY, node.SupportDZ = True, True, True

    # Analyzes the wall
    def analyze(self):
        """
        Analyze the network

        Args:
            self: (todo): write your description
        """

        # Descritize the wall
        self.__descritize()

        # Add supports to the wall
        self.__define_supports()

        # Add plate loads to the model
        for plate in self.fem.Plates:

            i_node = plate.iNode
            j_node = plate.jNode
            m_node = plate.mNode
            n_node = plate.nNode

            # Calculate the average Y coordinate of the four nodes
            avgY = (i_node.Y + j_node.Y + m_node.Y + n_node.Y)/4

            # Add plate surface loads to the model
            pressure = 0.0
            for load in self.loads:
                
                y1 = load[0]
                y2 = load[1]
                p1 = load[2]
                p2 = load[3]

                # Calculate the pressure on the plate apply it to the plate
                if avgY <= y2 and avgY >= y1:

                    # Calculate the pressure the plate
                    pressure += (p2 - p1)/(y2 - y1)*(avgY - y1) + p1

            # Add the surface pressure to the plate
            plate.pressures.append([pressure, 'Case 1'])

        # Analyze the model. The stiffness matrix will usually be sparse so some speed can usually
        # be gained by using the sparse solver.
        self.fem.Analyze(sparse=True)

        # Find the maximum displacement
        DZ = self.fem.Nodes[0].DZ['Combo 1']
        for node in self.fem.Nodes:
            if abs(node.DZ['Combo 1']) > abs(DZ):
                DZ = node.DZ['Combo 1']

        print('Max Displacement:', DZ)

    # Creates a contour plot of the wall forces
    def plot_forces(self, force_type):
        """
        Returns a plot of the wall's internal forces
        """
        
        # Determine the total number of nodes in the wall
        num_nodes = len(self.fem.Nodes)

        # Create a list of unique node X-coordinates
        x = []
        for i in range(num_nodes):
            x.append(self.fem.Nodes[i].X)
            if self.fem.Nodes[i+1].X < self.fem.Nodes[i].X:
                break
        
        # Get the number of columns of nodes
        num_cols = len(x)
        
        # Determine how many rows of nodes
        num_rows = int(round(num_nodes/num_cols, 0))

        # Initialize a list of unique Y-coordinates
        y = [self.fem.Nodes[0].Y]
        y_prev = self.fem.Nodes[0].Y

        # Initialize the list of node force results
        z = []

        # Determine which index in the 'M' vector we're interested in
        if force_type == 'Mx':
            index = 0
        elif force_type == 'My':
            index = 1
        elif force_type == 'Mxy':
            index = 2
        elif force_type == 'Qx':
            index = 0
        elif force_type == 'Qy':
            index = 1

        # Step through each node
        for node in self.fem.Nodes:

            # Add unique node Y-coordinates as we go
            if node.Y > y_prev:
                y.append(node.Y)
                y_prev = node.Y
            
            # Initialize the force at the node to zero
            force = 0
            count = 0

            # Find the plates that attach to the node
            for plate in self.fem.Plates:

                # Sum plate corner forces at the node
                if plate.iNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += plate.shear(0, 0, 'Combo 1')[index][0]
                    else:
                        force += plate.moment(0, 0, 'Combo 1')[index][0]
                    count += 1
                elif plate.nNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += plate.shear(0, plate.height(), 'Combo 1')[index][0]
                    else:
                        force += plate.moment(0, plate.height(), 'Combo 1')[index][0]
                    count += 1
                elif plate.mNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += plate.shear(plate.width(), plate.height(), 'Combo 1')[index][0]
                    else:
                        force += plate.moment(plate.width(), plate.height(), 'Combo 1')[index][0]
                    count += 1
                elif plate.jNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += plate.shear(plate.width(), 0, 'Combo 1')[index][0]
                    else:
                        force += plate.moment(plate.width(), 0, 'Combo 1')[index][0]
                    count += 1

            # Calculate the average force at the node
            force = force/count

            # Add the total force at the node to the list of forces
            z.append(force)

        # Convert the lists to numpy arrays
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        # Create a meshgrid for the X and Y-coordinates
        X, Y = np.meshgrid(x, y)

        # Reshape the node force results to fit the meshgrid
        z = z.reshape(num_rows, num_cols)

        fig, ax = plt.subplots()
        cs = ax.contourf(X, Y, z) #, 5)
        plt.colorbar(cs)
        ax.set(xlim=[0, max(x)], ylim=[0, max(y)], aspect=1)
        ax.set_title(force_type)
        plt.show()

    def plot_disp(self):
        """
        Plots the 2d plot

        Args:
            self: (todo): write your description
        """
        
        # Determine the total number of nodes in the wall
        num_nodes = len(self.fem.Nodes)

        # Create a list of unique node X-coordinates
        x = []
        for i in range(num_nodes):
            x.append(self.fem.Nodes[i].X)
            if self.fem.Nodes[i+1].X < self.fem.Nodes[i].X:
                break
        
        # Get the number of columns of nodes
        num_cols = len(x)
        
        # Determine how many rows of nodes
        num_rows = int(round(num_nodes/num_cols, 0))

        # Initialize a list of unique Y-coordinates
        y = [self.fem.Nodes[0].Y]
        y_prev = self.fem.Nodes[0].Y

        # Initialize the list of node displacement results
        z = []

        # Step through each node
        for node in self.fem.Nodes:

            # Add unique node Y-coordinates as we go
            if node.Y > y_prev:
                y.append(node.Y)
                y_prev = node.Y

            disp = node.DZ['Combo 1']

            # Add the total force at the node to the list of forces
            z.append(disp)

        # Convert the lists to numpy arrays
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        # Create a meshgrid for the X and Y-coordinates
        X, Y = np.meshgrid(x, y)

        # Reshape the node force results to fit the meshgrid
        z = z.reshape(num_rows, num_cols)

        fig, ax = plt.subplots()
        cs = ax.contourf(X, Y, z) #, 5)
        plt.colorbar(cs)
        ax.set(xlim=[0, max(x)], ylim=[0, max(y)], aspect=1)
        ax.set_title('Displacement')
        plt.show()

#%%
# Rectangular wall panel implementation example
# ---------------------------------------------
E = 57000*(4500)**0.5*144 # psf
t = 1 # ft
width = 10 # ft
height = 10 # ft
nu = 0.17
meshsize = 0.5 # ft
load = 250 # psf

myWall = RectWall(width, height, t, E, nu, meshsize, 'Fixed', 'Fixed', 'Fixed', 'Fixed')
myWall.add_load(0, height, load, load)

# Analyze the wall
import cProfile
cProfile.run('myWall.analyze()', sort='cumtime')

# Render the wall. The default load combination 'Combo 1' will be displayed since we're not specifying otherwise.
# The plates will be set to show the 'Mx' results
from PyNite import Visualization
Visualization.RenderModel(myWall.fem, text_height=meshsize/6, render_loads=True, color_map='Mx')

# Plot the displacement contour
myWall.plot_disp()

# Plot the shear force contours
myWall.plot_forces('Qx')
myWall.plot_forces('Qy')

# Plot the moment contours
myWall.plot_forces('Mx')
myWall.plot_forces('My')
myWall.plot_forces('Mxy')

# Results from Timoshenko's "Theory of Plates and Shells" textbook
D = E*t**3/(12*(1-nu**2))
print('Expected displacement from Timoshenko: ', 0.00126*load*width**4/D)

# Create a PDF report
# It will be output to the PyNite folder unless the 'output_path' variable below is modified
from PyNite import Reporting
Reporting.CreateReport(myWall.fem, output_filepath='.//PyNite Report.pdf', members=False, member_releases=False, \
                       member_end_forces=False, member_internal_forces=False)
