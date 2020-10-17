#%%
# This example works, but visualization and reporting for quads is not complete yet.

# This is an example of a how quads can be used to create a wall panel for out-of-plane bending problems.
# A 'RectWall' class has been defined to automate the process of meshing, loading, and plotting results for a wall
# panel of any geometry and edge support conditions. Below the class definition is a brief script showing how the
# 'RectWall' class can be been implemented.

#%%
from PyNite import FEModel3D, Node3D, Quad3D

# Packages used for plotting contours
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

#%%
class RectWall():

    # Constructor
    def __init__(self, width, height, thickness, E=3824, nu=0.3, mesh_size=1, bot_support='Pinned', top_support='Pinned', left_support='Free', right_support='Free'):

        self.width = width
        self.height = height
        self.thickness = thickness
        self.mesh_size = mesh_size

        self.E = E
        self.nu = nu

        self.fem = FEModel3D()  # A finite element model for the wall

        self.loads = []   # A list of surface loads applied to the wall panel
        self.nodes = []   # A list of nodes that make up the wall panel
        self.quads = []   # A list of quadrilaterals that make up the wall panel

        self.bot_support = bot_support
        self.top_support = top_support
        self.left_support = left_support
        self.right_support = right_support

        self.__analyzed = False
    
    # Adds a pressure load to the wall
    def add_load(self, y_bot, y_top, p_bot, p_top):
        
        # Add the load to the wall
        self.loads.append((y_bot, y_top, p_bot, p_top))
        self.__analyzed = False
    
    # Removes all loads from the wall
    def clear_loads(self):

        # Clear all the wall loads
        self.loads = []
        self.__analyzed = False
        self.fem.ClearLoads()
    
    # Descritizes the wall
    def __descritize(self):

        # Determine how many columns of quads are necessary
        num_cols = round(self.width/self.mesh_size)

        # Determine the quad width
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

        # Initialize the quad and node ID's to 1
        pl_id = 1
        node_id = 1

        # Create the bottom row of nodes
        for i in range(num_cols + 1):
            self.fem.AddNode('N'+str(node_id), i*pl_width, 0.0, 0.0)
            node_id += 1

        # Add quads and nodes between each control point
        for control_point in iter_y:

            # Determine the height between control points ("lift" height)
            h = control_point - y_prev

            # Determine how many rows of quads are needed in the lift
            num_rows = max(1, round(h/self.mesh_size))

            # Determine the quad height
            pl_height = h/num_rows

            # Generate nodes
            for j in range(num_rows + 1):
                for i in range(num_cols + 1):
                    
                    # The first row of nodes in the lift has already been generated
                    if j != 0:

                        # Generate the node
                        self.fem.AddNode('N'+str(node_id), i*pl_width, j*pl_height + y_prev, 0.0)
                        node_id += 1
            
            # Generate quadrilaterals
            for j in range(num_rows):
                for i in range(num_cols):

                    # Determine which nodes the quadrilateral will be attached to
                    ni = pl_id + max(0, int((pl_id-1)/num_cols))
                    nj = ni + 1
                    nm = nj + (num_cols + 1)
                    nn = nm - 1

                    # Add the quadrilateral to the list of quadrilaterals
                    self.fem.AddQuad('P'+str(pl_id), 'N'+str(ni), 'N'+str(nj), 'N'+str(nm), 'N'+str(nn), self.thickness, self.E, self.nu) 

                    # Prepare to move to the next iteration
                    pl_id += 1
            
            # Prepare for the next iteration
            y_prev = control_point

    # Defines support conditions at each node
    def __define_supports(self):

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

        # Descritize the wall
        self.__descritize()

        # Add supports to the wall
        self.__define_supports()

        # Add quad loads to the model
        for quad in self.fem.Quads:

            i_node = quad.iNode
            j_node = quad.jNode
            m_node = quad.mNode
            n_node = quad.nNode

            # Calculate the average Y coordinate of the four nodes
            avgY = (i_node.Y + j_node.Y + m_node.Y + n_node.Y)/4

            # Add quad surface loads to the model
            pressure = 0
            for load in self.loads:
                
                y1 = load[0]
                y2 = load[1]
                p1 = load[2]
                p2 = load[3]

                # Calculate the pressure on the quad and the load it applied to each of its nodes
                if avgY <= y2 and avgY >= y1:

                    # Calculate the pressure the quadrilateral
                    pressure += (p2 - p1)/(y2 - y1)*(avgY - y1) + p1
                    
                    # Add the surface pressure to the quadrilateral
                    quad.pressures.append([pressure, 'Case 1'])

        # Analyze the model
        self.fem.Analyze()

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

            # Find the quadrilaterals that attach to the node
            for quad in self.fem.Quads:

                # Sum quad corner forces at the node
                if quad.iNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += quad.shear(-1, -1, 'Combo 1')[index][0]
                    else:
                        force += quad.moment(-1, -1, 'Combo 1')[index][0]
                    count += 1
                elif quad.jNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += quad.shear(1, -1, 'Combo 1')[index][0]
                    else:
                        force += quad.moment(1, -1, 'Combo 1')[index][0]
                    count += 1
                elif quad.mNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += quad.shear(1, 1, 'Combo 1')[index][0]
                    else:
                        force += quad.moment(1, 1, 'Combo 1')[index][0]
                    count += 1
                elif quad.nNode.ID == node.ID:
                    if force_type == 'Qx' or force_type == 'Qy':
                        force += quad.shear(-1, 1, 'Combo 1')[index][0]
                    else:
                        force += quad.moment(-1, 1, 'Combo 1')[index][0]
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
E = 57*(4500)**0.5 # ksi
t = 12 # in
width = 10*12 # in
height = 20*12 #in
nu = 0.17
meshsize = 6 # in
load = 0.250/144 # ksi

myWall = RectWall(width, height, t, E, nu, meshsize, 'Fixed', 'Fixed', 'Fixed', 'Fixed')
myWall.add_load(0, height, load, load)

# Analyze the wall
myWall.analyze()

# Render the wall. The default load combination 'Combo 1' will be displayed since we're not specifying otherwise.
# from PyNite import Visualization
# Visualization.RenderModel(myWall.fem, text_height=meshsize/6, render_loads=True)

# Results from Timoshenko's "Theory of Plates and Shells" Table 35, p. 202
D = E*t**3/(12*(1-nu**2))
print('Solution from Timoshenko Table 35 for b/a = 2.0:')
print('Expected displacement: ', 0.00254*load*width**4/D)
print('Mx at Center:', 0.0412*load*width**2)
print('Mx at Edges:', -0.0829*load*width**2)
print('My at Center:', 0.0158*load*width**2)
print('My at Top & Bottom:', -0.0571*load*width**2)

print('')
print('Max Mx', max([max(quad.moment(-1, -1)[0], quad.moment(1, -1)[0], quad.moment(1, 1)[0], quad.moment(-1, 1)[0]) for quad in myWall.fem.Quads]))
print('Min Mx', min([min(quad.moment(-1, -1)[0], quad.moment(1, -1)[0], quad.moment(1, 1)[0], quad.moment(-1, 1)[0]) for quad in myWall.fem.Quads]))
print('Max My', max([max(quad.moment(-1, -1)[1], quad.moment(1, -1)[1], quad.moment(1, 1)[1], quad.moment(-1, 1)[1]) for quad in myWall.fem.Quads]))
print('Min My', min([min(quad.moment(-1, -1)[1], quad.moment(1, -1)[1], quad.moment(1, 1)[1], quad.moment(-1, 1)[1]) for quad in myWall.fem.Quads]))
print('Max Mxy', max([max(quad.moment(-1, -1)[2], quad.moment(1, -1)[2], quad.moment(1, 1)[2], quad.moment(-1, 1)[2]) for quad in myWall.fem.Quads]))
print('Min Mxy', min([min(quad.moment(-1, -1)[2], quad.moment(1, -1)[2], quad.moment(1, 1)[2], quad.moment(-1, 1)[2]) for quad in myWall.fem.Quads]))


print('Max Qx', max([max(quad.shear(-1, -1)[0], quad.shear(1, -1)[0], quad.shear(1, 1)[0], quad.shear(-1, 1)[0]) for quad in myWall.fem.Quads]))
print('Min Qx', min([min(quad.shear(-1, -1)[0], quad.shear(1, -1)[0], quad.shear(1, 1)[0], quad.shear(-1, 1)[0]) for quad in myWall.fem.Quads]))
print('Max Qy', max([max(quad.shear(-1, -1)[1], quad.shear(1, -1)[1], quad.shear(1, 1)[1], quad.shear(-1, 1)[1]) for quad in myWall.fem.Quads]))
print('Max Qy', min([min(quad.shear(-1, -1)[1], quad.shear(1, -1)[1], quad.shear(1, 1)[1], quad.shear(-1, 1)[1]) for quad in myWall.fem.Quads]))

# Plot the displacement contour
myWall.plot_disp()

# Plot the moment contours
myWall.plot_forces('Mx')
myWall.plot_forces('My')
myWall.plot_forces('Mxy')

# Plot the shear force contours
myWall.plot_forces('Qx')
myWall.plot_forces('Qy')

# Create a PDF report
# It will be output to the PyNite folder unless the 'output_path' variable below is modified
# from PyNite import Reporting
# Reporting.CreateReport(myWall.fem, output_filepath='.//PyNite Report.pdf', members=False, member_releases=False, \
#                        member_end_forces=False, member_internal_forces=False)