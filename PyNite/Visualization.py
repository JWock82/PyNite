
from json import load
import warnings

from IPython.display import Image
from numpy import array, empty, append, cross
from numpy.linalg import norm
import vtk

class Renderer():
    """Used to render finite element models.
    """

    scalar = None

    def __init__(self, model):

        self.model = model

        # Default settings for rendering
        self.annotation_size = 5
        self.deformed_shape = False
        self.deformed_scale = 30
        self.render_nodes = True
        self.render_loads = True
        self.show_plate_edges = True
        self.color_map = None
        self.num_contours = -1 #-1 will provide smooth contours, values greater than 0 will create banded contours
        self.contour_range = (None, None) #(minval, maxval)
        self.combo_name = 'Combo 1'
        self.case = None
        self.labels = True
        self.scalar_bar = False
        self.scalar_bar_text_size = 24
        self.theme = 'default'

        # Initialize VTK objects
        self.renderer = vtk.vtkRenderer()
        self.window = vtk.vtkRenderWindow()
        self.window.SetWindowName('PyNite - Simple Finite Element Analysis in Python')
        self.window.AddRenderer(self.renderer)

    @property
    def window_width(self):
        return self.window.GetSize()[0]
    
    @window_width.setter
    def window_width(self, width):
        height = self.window.GetSize()[1]
        self.window.SetSize(width, height)
    
    @property
    def window_height(self):
        return self.window.GetSize()[1]
    
    @window_height.setter
    def window_height(self, height):
        width = self.window.GetSize()[0]
        self.window.SetSize(width, height)

    def set_annotation_size(self, size=5):
        self.annotation_size = size
    
    def set_deformed_shape(self, deformed_shape=False):
        self.deformed_shape = deformed_shape
    
    def set_deformed_scale(self, scale=30):
        self.deformed_scale = scale
    
    def set_render_nodes(self, render_nodes=True):
        self.render_nodes = render_nodes

    def set_render_loads(self, render_loads=True):
        self.render_loads = render_loads
    
    def set_color_map(self, color_map=None):
        """Sets the color map for plate contours.

        :param color_map: The color map to use: Valid options are 'Qx', Qy', 'Mx', 'My', 'Mxy', 'Sx', 'Sy' 'Txy'. Qx and Qy are out-of-plane shear forces. Mx, My, and Mxy are local out-of-plane bending moments. 'Sx' and 'Sy' are membrane forces, and 'Txy' is an in-plane shear force. Defaults to None.
        :type color_map: str, optional
        """
        self.color_map = color_map
        
    def set_contour_range(self, minval=None, maxval=None):
        self.contour_range = (minval, maxval)
        
    def set_number_of_contours(self, num_contours = -1):
        self.num_contours = num_contours
    
    def set_combo_name(self, combo_name='Combo 1'):
        self.combo_name = combo_name
        self.case = None
    
    def set_case(self, case=None):
        self.case = case
        self.combo_name = None
    
    def set_show_labels(self, show_labels=True):
        self.labels = show_labels
    
    def set_scalar_bar(self, scalar_bar=False):
        self.scalar_bar = scalar_bar
    
    def set_scalar_bar_text_size(self, text_size=24):
        self.scalar_bar_text_size = text_size

    def render_model(self, interact=True, reset_camera=True):
        """
        Renders the model in a window

        Parameters
        ----------
        interact : bool
            Suppresses interacting with the window if set to `False`. This can be used to capture a
            screenshot without pausing the program for the user to interact. Default is `True`.
        reset_camera : bool
            Resets the camera if set to `True`. Default is `True`.
        """
        
        # Get the render window
        window = self.window
        
        # Update the renderer
        self.update(reset_camera)

        # Render the window
        window.Render()

        # Handle user interaction if requested by the user
        if interact == True:
        
            # Set up an interactor. The interactor style determines how user interactions affect the
            # view. The trackball camera style behaves much like popular commercial CAD programs.
            interactor = vtk.vtkRenderWindowInteractor()
            style = vtk.vtkInteractorStyleTrackballCamera()
            interactor.SetInteractorStyle(style)
            interactor.SetRenderWindow(self.window)
            
            # Start the interactor. Code execution will pause here until the user closes the window.
            interactor.Start()
            
            # Finalize the render window once the user closes out of it. I don't understand everything
            # this does, but I've found screenshots will cause the program to crash if this line is
            # omitted. I have noticed it will shut down the interactor.
            window.Finalize()

        # Return the window
        return window

    def screenshot(self, filepath='console', interact=True, reset_camera=True):
        """
        Renders the model in a window. When the window is closed a screenshot is captured.

        Parameters
        ----------
        filepath : string
            Sends a screenshot to the specified filepath. The screenshot will be taken when the
            user closes out of the render window. If `filepath` is set to 'console' the screenshot
            will be returned as an IPython image. If set to 'BytesIO' it will return the image as a
            `BytesIO` object. Default is 'console'.
        interact : bool
            Suppresses interacting with the window if set to `False`. This can be used to capture a
            screenshot without pausing the program for the user to interact. Default is `True`.
        reset_camera : boolju
            Resets the camera if set to `True`. Default is `True`.
        """

        # Render the model in a window and save the window
        window = self.render_model(interact, reset_camera)

        # Screenshot code
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(window)
        w2if.SetInputBufferTypeToRGB()
        w2if.ReadFrontBufferOff()

        # These next two lines are in the examples and documentation for VTK, but don't seem to do
        # anything. I've left them here in case I find a bug somewhere down the line that needs
        # fixing.
        # w2if.Update()
        # w2if.Modified()

        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(w2if.GetOutputPort())

        if filepath == 'console' or filepath == 'BytesIO':

            writer.SetWriteToMemory(1)
            writer.Write()
            fig_file = memoryview(writer.GetResult()).tobytes()

            # Now that we're done with the render window, finalize it
            window.Finalize()

            if filepath == 'console':
                return Image(fig_file)
            elif filepath == 'BytesIO':
                from io import BytesIO
                return BytesIO(fig_file)
        else:

            writer.SetFileName(filepath)
            writer.Write()

            # Now that we're done with the render window, finalize it
            window.Finalize()

            return

    def update(self, reset_camera=True):
        """
        Builds or rebuilds the VTK renderer

        Parameters
        ----------
        reset_camera : bool
            Resets the camera if set to `True`. Default is `True`.
        """

        # Input validation
        if self.deformed_shape and self.case != None:
            raise Exception('Deformed shape is only available for load combinations,'
                            ' not load cases.')
        if self.model.LoadCombos == {} and self.render_loads == True and self.case == None:
            self.render_loads = False
            warnings.warn('Unable to render load combination. No load combinations defined.', UserWarning)
        
        # Check if nodes are to be rendered
        if self.render_nodes == True:

            if self.theme == 'print':
                color = 'black'
            else:
                color = None

            # Create a visual node for each node in the model
            vis_nodes = []
            for node in self.model.Nodes.values():
                vis_nodes.append(VisNode(node, self.annotation_size, color))
            
            # Create a visual auxiliary node for each auxiliary node in the model
            vis_aux_nodes = []
            for aux_node in self.model.AuxNodes.values():
                vis_aux_nodes.append(VisNode(aux_node, self.annotation_size, color))
        
        # Create a visual spring for each spring in the model
        vis_springs = []
        for spring in self.model.Springs.values():
            vis_springs.append(VisSpring(spring, self.model.Nodes, self.annotation_size))    
    
        # Create a visual member for each member in the model
        vis_members = []
        for member in self.model.Members.values():
            vis_members.append(VisMember(member, self.model.Nodes, self.annotation_size, self.theme))
        
        # Get the renderer
        renderer = self.renderer

        # Clear out all the old actors from any previous renderings
        for actor in renderer.GetActors():
            renderer.RemoveActor(actor)

        # Add actors for each spring
        for vis_spring in vis_springs:
        
            # Add the actor for the spring
            renderer.AddActor(vis_spring.actor)

            if self.labels == True:
                # Add the actor for the spring label
                renderer.AddActor(vis_spring.lblActor)
            
                # Set the text to follow the camera as the user interacts. This will
                # require a reset of the camera (see below)
                vis_spring.lblActor.SetCamera(renderer.GetActiveCamera())    

        # Add actors for each member
        for vis_member in vis_members:
            
            # Add the actor for the member
            renderer.AddActor(vis_member.actor)

            if self.labels == True:
                
                # Add the actor for the member label
                renderer.AddActor(vis_member.lblActor)
            
                # Set the text to follow the camera as the user interacts. This will
                # require a reset of the camera (see below)
                vis_member.lblActor.SetCamera(renderer.GetActiveCamera())
        
        # Check if nodes are to be rendered
        if self.render_nodes == True:
            
            # Combine the polydata from each node

            # Create an append filter for combining node polydata
            node_polydata = vtk.vtkAppendPolyData()

            for vis_node in vis_nodes:
                
                # Add the node's polydata
                node_polydata.AddInputData(vis_node.polydata.GetOutput())

                if self.labels == True:

                    if self.theme == 'print':
                        vis_node.lblActor.GetProperty().SetColor(0, 0, 0)  # black
                    
                    # Add the actor for the node label
                    renderer.AddActor(vis_node.lblActor)
                
                    # Set the text to follow the camera as the user interacts. This will
                    # require a reset of the camera (see below)
                    vis_node.lblActor.SetCamera(renderer.GetActiveCamera())
            
            # Update the node polydata in the append filter
            node_polydata.Update()
            
            # Create a mapper and actor for the nodes
            node_mapper = vtk.vtkPolyDataMapper()
            node_mapper.SetInputConnection(node_polydata.GetOutputPort())
            node_actor = vtk.vtkActor()
            node_actor.SetMapper(node_mapper)

            # Adjust the color of the nodes here.
            if self.theme == 'print':
                node_actor.GetProperty().SetColor(0, 0, 0)  # Black
            
            # Add the node actor to the renderer
            renderer.AddActor(node_actor)

            # Add actors for each auxiliary node. There are usually only a few of these, so they shouldn't bog down interaction if we do them individually rather than with a `vtkAppendFilter`
            for vis_aux_node in vis_aux_nodes:

                # Add the actor for the auxiliary node
                renderer.AddActor(vis_aux_node.actor)

                if self.labels == True:
                    
                    # Add the actor for the auxiliary node label
                    renderer.AddActor(vis_aux_node.lblActor)
            
                    # Set the text to follow the camera as the user interacts. This will
                    # require a reset of the camera (see below)
                    vis_aux_node.lblActor.SetCamera(renderer.GetActiveCamera())

        # Render the deformed shape if requested
        if self.deformed_shape == True:
            _DeformedShape(self.model, renderer, self.deformed_scale, self.annotation_size, self.combo_name, self.render_nodes, self.theme)

        # Render the loads if requested
        if (self.combo_name != None or self.case != None) and self.render_loads != False:
            _RenderLoads(self.model, renderer, self.annotation_size, self.combo_name, self.case, self.theme)
        
        # Render the plates and quads, if present
        if self.model.Quads or self.model.Plates:
            _RenderContours(self.model, renderer, self.deformed_shape, self.deformed_scale,
                            self.color_map, self.scalar_bar, self.scalar_bar_text_size,
                            self.combo_name, self.show_plate_edges, self.num_contours, self.contour_range,
                            self.theme)

        # Set the window's background color
        if self.theme == 'default':
            renderer.SetBackground(0, 0, 128)  # Blue
        elif self.theme == 'print':
            renderer.SetBackground(255, 255, 255)  # White
        
        # Reset the camera
        if reset_camera: renderer.ResetCamera()

#%%
# The code in this section will be deprecated at some point

def RenderModel(model, annotation_size=5, deformed_shape=False, deformed_scale=30,
                render_loads=True, color_map=None, combo_name='Combo 1', case=None, labels=True,
                screenshot=None):
    warnings.warn('`RenderModel` will be replaced with `render_model` in a future version of PyNite.', FutureWarning)
    render_model(model, annotation_size, deformed_shape, deformed_scale, render_loads, color_map,
                 True, combo_name, case, labels, screenshot)

def render_model(model, annotation_size=5, deformed_shape=False, deformed_scale=30, render_loads=True, color_map=None, scalar_bar=True, combo_name='Combo 1', case=None, labels=True, screenshot=None, theme='default'):
    '''
    Renders a finite element model using VTK.
    
    Parameters
    ----------
    model : FEModel3D
        Finite element model to be rendered.
    annotation_size : number, optional
        Controls the height of text displayed with the model. The units used for
        `annotation_size` are the same as those used for lengths in the model. Sizes
        of other objects (such as nodes) are related to this value. The default is
        5.
    deformed_shape : boolean, optional
        Determines whether the deformed shape will be rendered or not. The model
        must be solved to use this feature. Deformed shapes are not available for
        load cases, only load combinations. The default is False.
    deformed_scale : boolean, optional
        Determines what magnification factor will be applied to the deformed
        shape. The default is 30.
    render_loads : boolean, optional
        Determines if loads will be rendered with the model. The default is True.
    color_map : string, optional
        The type of contour to plot. Acceptable values are: 'dz', 'Mx', 'My', 'Mxy', 'Qx', 'Qy',
        'Sx', 'Sy', 'Txy'
        If no value is specified the default is None.
    combo_name : string, optional
        The load combination used for rendering the deformed shape and the loads.
        The default is 'Combo 1'.
    case : string, optional
        The load case used for rendering loads. The default is None.
    labels : boolean, optional
        Determines if labels will be rendered. Each label is a single actor in VTK, which slows down
        rendering on models with thousands of labels. Set this option to `False` if you want more
        speed when rendering and interacting with a large model. This can be very useful on plate
        models with large meshes.
    screenshot : string, optional
        Sends a screenshot to the specified filepath unless set to None. The screenshot will be taken
        when the user closes out of the render window. If screenshot is set to 'console' the
        screenshot will be returned as an IPython image. Default is None.
      
    Raises
    ------
    Exception
        A deformed shape is requested and a load case has been specified.
    
    Returns
    -------
    Ipython Image (if screenshot is set to 'console')
    
    '''
    
    # Input validation
    if deformed_shape and case != None:
        raise Exception('Deformed shape is only available for load combinations,'
                        ' not load cases.')
    if model.LoadCombos == {} and render_loads == True and case == None:
        raise Exception('Unable to render load combination. No load combinations defined.')
    
    # Create a visual node for each node in the model
    if theme == 'print': color = 'black'
    else: color = None
    vis_nodes = []
    for node in model.Nodes.values():
        vis_nodes.append(VisNode(node, annotation_size, color))
    
    # Create a visual auxiliary node for each auxiliary node in the model
    vis_aux_nodes = []
    for aux_node in model.AuxNodes.values():
        vis_aux_nodes.append(VisNode(aux_node, annotation_size, color='red'))
    
    # Create a visual spring for each spring in the model
    vis_springs = []
    for spring in model.Springs.values():
        vis_springs.append(VisSpring(spring, model.Nodes, annotation_size))    
  
    # Create a visual member for each member in the model
    vis_members = []
    for member in model.Members.values():
        vis_members.append(VisMember(member, model.Nodes, annotation_size, theme))
    
    # Create a window
    window = vtk.vtkRenderWindow()
    
    # Set the pixel width and length of the window
    window.SetSize(750, 750)
    
    # Set up the interactor. The interactor style determines how user
    # interactions affect the view. The trackball camera style behaves much
    # like popular commercial CAD programs.
    interactor = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)
    interactor.SetRenderWindow(window)
    
    # Create a renderer object and add it to the window
    renderer = vtk.vtkRenderer()
    window.AddRenderer(renderer)

    # Add actors for each spring
    for vis_spring in vis_springs:
      
        # Add the actor for the spring
        renderer.AddActor(vis_spring.actor)

        if labels == True:
            # Add the actor for the spring label
            renderer.AddActor(vis_spring.lblActor)
          
            # Set the text to follow the camera as the user interacts. This will
            # require a reset of the camera (see below)
            vis_spring.lblActor.SetCamera(renderer.GetActiveCamera())    

    # Add actors for each member
    for vis_member in vis_members:
        
        # Add the actor for the member
        renderer.AddActor(vis_member.actor)

        if labels == True:
            # Add the actor for the member label
            renderer.AddActor(vis_member.lblActor)
          
            # Set the text to follow the camera as the user interacts. This will
            # require a reset of the camera (see below)
            vis_member.lblActor.SetCamera(renderer.GetActiveCamera())

    # Combine the polydata from each node

    # Create an append filter for combining node polydata
    node_polydata = vtk.vtkAppendPolyData()

    for vis_node in vis_nodes:
        
        # Add the node's polydata
        node_polydata.AddInputData(vis_node.polydata.GetOutput())

        if labels == True:
            
            # Add the actor for the node label
            renderer.AddActor(vis_node.lblActor)
        
            # Set the text to follow the camera as the user interacts. This will
            # require a reset of the camera (see below)
            vis_node.lblActor.SetCamera(renderer.GetActiveCamera())
      
    # Update the node polydata in the append filter
    node_polydata.Update()
    
    # Create a mapper and actor for the nodes
    node_mapper = vtk.vtkPolyDataMapper()
    node_mapper.SetInputConnection(node_polydata.GetOutputPort())
    node_actor = vtk.vtkActor()
    node_actor.SetMapper(node_mapper)
    
    # Add the node actor to the renderer
    renderer.AddActor(node_actor)

    # Add actors for each auxiliary node
    for vis_aux_node in vis_aux_nodes:
      
        # Add the actor for the auxiliary node
        renderer.AddActor(vis_aux_node.actor)

        if labels == True:
            
            # Add the actor for the auxiliary node label
            renderer.AddActor(vis_aux_node.lblActor)
      
            # Set the text to follow the camera as the user interacts. This will
            # require a reset of the camera (see below)
            vis_aux_node.lblActor.SetCamera(renderer.GetActiveCamera())

    # Render the deformed shape if requested
    if deformed_shape == True:
        _DeformedShape(model, renderer, deformed_scale, annotation_size, combo_name, render_nodes=True, theme=theme)

    # Render the loads if requested
    if (combo_name != None or case != None) and render_loads != False:
        _RenderLoads(model, renderer, annotation_size, combo_name, case)
    
    # Render the plates and quads, if present
    if model.Quads or model.Plates:
        _RenderContours(model, renderer, deformed_shape, deformed_scale, color_map, scalar_bar,
                        12, combo_name)

    # Set the window's background to gray
    renderer.SetBackground(0/255, 0/255, 128/255)
    
    # Reset the camera
    renderer.ResetCamera()
    
    # Render the window
    window.Render()

    # Start the interactor. Code execution will pause here until the user closes the window.
    interactor.Start()

    # Finalize the render window once the user closes out of it. I don't understand everything
    # this does, but I've found screenshots will cause the program to crash if this line is
    # omitted.
    window.Finalize()

    # Create a screenshot of the last view in the window (if requested)
    if screenshot != None:
        
        # Screenshot code
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(window)
        w2if.SetInputBufferTypeToRGB()
        w2if.ReadFrontBufferOff()

        # These next two lines are in the examples and documentation for VTK, but don't seem to do
        # anything. I've left them here in case I find a bug somewhere down the line that needs
        # fixing.
        # w2if.Update()
        # w2if.Modified()

        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(w2if.GetOutputPort())

        if screenshot == 'console':
            writer.SetWriteToMemory(1)
            writer.Write()
            fig_file = memoryview(writer.GetResult()).tobytes()
            return Image(fig_file)
        else:
            writer.SetFileName(screenshot)
            writer.Write()

#%%
# Converts a node object into a node for the viewer
class VisNode():

    # Constructor
    def __init__(self, node, annotation_size=5, color=None):
      
        # Create an append filter to append all the sources related to the node into a single 'PolyData' object
        self.polydata = vtk.vtkAppendPolyData()
      
        # Get the node's position
        X = node.X # Global X coordinate
        Y = node.Y # Global Y coordinate
        Z = node.Z # Global Z coordinate
      
        # Generate a sphere source for the node
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(X, Y, Z)
        sphere.SetRadius(0.6*annotation_size)
        sphere.Update()
        self.polydata.AddInputData(sphere.GetOutput())
      
        # Create the text for the node label
        label = vtk.vtkVectorText()
        label.SetText(node.name)
        
        # Set up a mapper for the node label
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())
      
        # Set up an actor for the node label
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition(X + 0.6*annotation_size, Y + 0.6*annotation_size, Z)
      
        # Generate any supports that occur at the node
        # Check for a fixed suppport
        if node.support_DX == True and node.support_DY == True and node.support_DZ == True \
        and node.support_RX == True and node.support_RY == True and node.support_RZ == True:
      
            # Create the fixed support
            support = vtk.vtkCubeSource()
            support.SetCenter(node.X, node.Y, node.Z)
            support.SetXLength(annotation_size*1.2)
            support.SetYLength(annotation_size*1.2)
            support.SetZLength(annotation_size*1.2)
        
            # Copy and append the support data to the append filter
            support.Update()
            self.polydata.AddInputData(support.GetOutput())
        
        # Check for a pinned support
        elif node.support_DX == True and node.support_DY == True and node.support_DZ == True \
        and node.support_RX == False and node.support_RY == False and node.support_RZ == False:
          
            # Create the pinned support
            support = vtk.vtkConeSource()
            support.SetCenter(node.X, node.Y-0.6*annotation_size, node.Z)
            support.SetDirection((0, 1, 0))
            support.SetHeight(annotation_size*1.2)
            support.SetRadius(annotation_size*1.2)
        
            # Copy and append the support data to the append filter
            support.Update()
            self.polydata.AddInputData(support.GetOutput())
        
        # Other support conditions
        else:
      
            # Restrained against X translation
            if node.support_DX == True:
                  
                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X-annotation_size, node.Y, node.Z)
                support1.SetPoint2(node.X+annotation_size, node.Y, node.Z)
          
                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())
          
                support2 = vtk.vtkConeSource()
                support2.SetCenter(node.X-annotation_size, node.Y, node.Z)
                support2.SetDirection((1, 0, 0))
                support2.SetHeight(annotation_size*0.6)
                support2.SetRadius(annotation_size*0.3)
          
                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())
          
                support3 = vtk.vtkConeSource()
                support3.SetCenter(node.X+annotation_size, node.Y, node.Z)
                support3.SetDirection((-1, 0, 0))
                support3.SetHeight(annotation_size*0.6)
                support3.SetRadius(annotation_size*0.3)
          
                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())
          
            # Restrained against Y translation
            if node.support_DY == True:
              
                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y-annotation_size, node.Z)
                support1.SetPoint2(node.X, node.Y+annotation_size, node.Z)
          
                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())
          
                support2 = vtk.vtkConeSource()
                support2.SetCenter(node.X, node.Y-annotation_size, node.Z)
                support2.SetDirection((0, 1, 0))
                support2.SetHeight(annotation_size*0.6)
                support2.SetRadius(annotation_size*0.3)
          
                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())
          
                support3 = vtk.vtkConeSource()
                support3.SetCenter(node.X, node.Y+annotation_size, node.Z)
                support3.SetDirection((0, -1, 0))
                support3.SetHeight(annotation_size*0.6)
                support3.SetRadius(annotation_size*0.3)
          
                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())
          
            # Restrained against Z translation
            if node.support_DZ == True:
            
                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y, node.Z-annotation_size)
                support1.SetPoint2(node.X, node.Y, node.Z+annotation_size)
          
                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())
          
                support2 = vtk.vtkConeSource()
                support2.SetCenter(node.X, node.Y, node.Z-annotation_size)
                support2.SetDirection((0, 0, 1))
                support2.SetHeight(annotation_size*0.6)
                support2.SetRadius(annotation_size*0.3)
          
                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())
          
                support3 = vtk.vtkConeSource()
                support3.SetCenter(node.X, node.Y, node.Z+annotation_size)
                support3.SetDirection((0, 0, -1))
                support3.SetHeight(annotation_size*0.6)
                support3.SetRadius(annotation_size*0.3)
          
                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())
      
            # Restrained against rotation about the X-axis
            if node.support_RX == True:
            
                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X-1.6*annotation_size, node.Y, node.Z)
                support1.SetPoint2(node.X+1.6*annotation_size, node.Y, node.Z)
          
                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())
          
                support2 = vtk.vtkCubeSource()
                support2.SetCenter(node.X-1.9*annotation_size, node.Y, node.Z)
                support2.SetXLength(annotation_size*0.6)
                support2.SetYLength(annotation_size*0.6)
                support2.SetZLength(annotation_size*0.6)
          
                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())
          
                support3 = vtk.vtkCubeSource()
                support3.SetCenter(node.X+1.9*annotation_size, node.Y, node.Z)
                support3.SetXLength(annotation_size*0.6)
                support3.SetYLength(annotation_size*0.6)
                support3.SetZLength(annotation_size*0.6)
          
                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())
          
            # Restrained against rotation about the Y-axis
            if node.support_RY == True:
            
                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y-1.6*annotation_size, node.Z)
                support1.SetPoint2(node.X, node.Y+1.6*annotation_size, node.Z)
          
                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())
          
                support2 = vtk.vtkCubeSource()
                support2.SetCenter(node.X, node.Y-1.9*annotation_size, node.Z)
                support2.SetXLength(annotation_size*0.6)
                support2.SetYLength(annotation_size*0.6)
                support2.SetZLength(annotation_size*0.6)
          
                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())
          
                support3 = vtk.vtkCubeSource()
                support3.SetCenter(node.X, node.Y+1.9*annotation_size, node.Z)
                support3.SetXLength(annotation_size*0.6)
                support3.SetYLength(annotation_size*0.6)
                support3.SetZLength(annotation_size*0.6)
          
                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())
          
            # Restrained against rotation about the Z-axis
            if node.support_RZ == True:
            
                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y, node.Z-1.6*annotation_size)
                support1.SetPoint2(node.X, node.Y, node.Z+1.6*annotation_size)
          
                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())
          
                support2 = vtk.vtkCubeSource()
                support2.SetCenter(node.X, node.Y, node.Z-1.9*annotation_size)
                support2.SetXLength(annotation_size*0.6)
                support2.SetYLength(annotation_size*0.6)
                support2.SetZLength(annotation_size*0.6)
          
                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())
          
                support3 = vtk.vtkCubeSource()
                support3.SetCenter(node.X, node.Y, node.Z+1.9*annotation_size)
                support3.SetXLength(annotation_size*0.6)
                support3.SetYLength(annotation_size*0.6)
                support3.SetZLength(annotation_size*0.6)
          
                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())
        
        # Update the append filter
        self.polydata.Update()
      
        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.polydata.GetOutputPort())
        self.actor = vtk.vtkActor()
        
        # Set the mapper for the node's actor
        self.actor.SetMapper(mapper)

        # Color will be added to the node actors outside of this class once they are all combined into one `vtkAppendFilter` object.

        # TODO: Delete legacy code below once it's been proven to be unnecessary over time.

        # Add color to the node and label actors if specified
        # if color == 'red':
        #     self.actor.GetProperty().SetColor(255, 0, 0)  # Red
        #     self.lblActor.GetProperty().SetColor(255, 0, 0)  # Red
        # elif color == 'yellow':
        #     self.actor.GetProperty().SetColor(255, 255, 0)  # Yellow
        #     self.lblActor.GetProperty().SetColor(255, 255, 0)  # Yellow
        # elif color == 'black':
        #     self.actor.GetProperty().SetColor(0, 0, 0)  # Black
        #     self.lblActor.GetProperty().SetColor(0, 0, 0)  # Black

class VisSpring():
    
    def __init__(self, spring, nodes, annotation_size=5, color=None):

        # Generate a line source for the spring
        line = vtk.vtkLineSource()

        # Step through each node in the model and find the position of the
        # i-node and j-node
        for node in nodes.values():

            # Check to see if the current node is the i-node
            if node.name == spring.i_node.name:
                Xi = node.X
                Yi = node.Y
                Zi = node.Z
                line.SetPoint1(Xi, Yi, Zi)

            # Check to see if the current node is the j-node
            elif node.name == spring.j_node.name:
                Xj = node.X
                Yj = node.Y
                Zj = node.Z
                line.SetPoint2(Xj, Yj, Zj)
    
        # Set up a mapper for the spring
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(line.GetOutputPort())

        # Set up an actor and a mapper for the spring
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)

        # Create the text for the spring label
        label = vtk.vtkVectorText()
        label.SetText(spring.name)

        # Set up a mapper for the spring label
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())

        # Set up an actor for the spring label
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition((Xi+Xj)/2, (Yi+Yj)/2, (Zi+Zj)/2)

        # Add some color
        if color is None:
            self.actor.GetProperty().SetColor(255, 0, 255)     # Magenta
            self.lblActor.GetProperty().SetColor(255, 0, 255)
        elif color == 'black':
            self.actor.GetProperty().SetColor(0, 0, 0)         # Black
            self.lblActor.GetProperty().SetColor(0, 0, 0)
      
# Converts a member object into a member for the viewer
class VisMember():

    # Constructor
    def __init__(self, member, nodes, annotation_size=5, theme='default'):
    
        # Generate a line for the member
        line = vtk.vtkLineSource()
      
        # Step through each node in the model and find the position of the i-node and j-node
        for node in nodes.values():
        
            # Check to see if the current node is the i-node
            if node.name == member.i_node.name:
                Xi = node.X
                Yi = node.Y
                Zi = node.Z
                line.SetPoint1(Xi, Yi, Zi)
        
            # Check to see if the current node is the j-node
            elif node.name == member.j_node.name:
                Xj = node.X
                Yj = node.Y
                Zj = node.Z
                line.SetPoint2(Xj, Yj, Zj)
        
        # Set up a mapper for the member
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(line.GetOutputPort())
      
        # Set up an actor for the member
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)
      
        # Create the text for the member label
        label = vtk.vtkVectorText()
        label.SetText(member.name)
      
        # Set up a mapper for the member label
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())
      
        # Set up an actor for the member label
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition((Xi+Xj)/2, (Yi+Yj)/2, (Zi+Zj)/2)

        # Adjust the color of the member if the theme is 'print'
        if theme == 'print':
            self.actor.GetProperty().SetColor(0/255, 0/255, 0/255)  # Black
            self.lblActor.GetProperty().SetColor(0/255, 0/255, 0/255) # Black

# Converts a node object into a node in its deformed position for the viewer
class VisDeformedNode():
    
    def __init__(self, node, scale_factor, annotation_size=5, combo_name='Combo 1'):
          
        # Calculate the node's deformed position
        newX = node.X + scale_factor*(node.DX[combo_name])
        newY = node.Y + scale_factor*(node.DY[combo_name])
        newZ = node.Z + scale_factor*(node.DZ[combo_name])
            
        # Generate a sphere source for the node in its deformed position
        self.source = vtk.vtkSphereSource()
        self.source.SetCenter(newX, newY, newZ)
        self.source.SetRadius(0.6*annotation_size)
        self.source.Update()

class VisDeformedMember():
    
    def __init__(self, member, nodes, scale_factor, combo_name='Combo 1'):
        
        # Determine if this member is active for each load combination
        self.active = member.active
      
        L = member.L() # Member length
        T = member.T() # Member local transformation matrix
      
        cos_x = array([T[0,0:3]]) # Direction cosines of local x-axis
        cos_y = array([T[1,0:3]]) # Direction cosines of local y-axis
        cos_z = array([T[2,0:3]]) # Direction cosines of local z-axis
      
        # Find the initial position of the local i-node
        # Step through each node
        for node in nodes.values():
          
            # Check to see if the current node is the i-node
            if node.name == member.i_node.name:
                Xi = node.X
                Yi = node.Y
                Zi = node.Z
      
        # Calculate the local y-axis displacements at 20 points along the member's
        # length
        DY_plot = empty((0, 3))
        for i in range(20):
                
            # Calculate the local y-direction displacement
            dy_tot = member.deflection('dy', L/19*i, combo_name)
        
            # Calculate the scaled displacement in global coordinates
            DY_plot = append(DY_plot, dy_tot*cos_y*scale_factor, axis=0)
      
        # Calculate the local z-axis displacements at 20 points along the member's
        # length
        DZ_plot = empty((0, 3)) 
        for i in range(20):
                
            # Calculate the local z-direction displacement
            dz_tot = member.deflection('dz', L/19*i, combo_name)
        
            # Calculate the scaled displacement in global coordinates
            DZ_plot = append(DZ_plot, dz_tot*cos_z*scale_factor, axis=0)
      
        # Calculate the local x-axis displacements at 20 points along the member's
        # length
        DX_plot = empty((0, 3)) 
        for i in range(20):
                
            # Displacements in local coordinates
            dx_tot = [[Xi, Yi, Zi]] + (L/19*i + member.deflection('dx', L/19*i, combo_name)*scale_factor)*cos_x
                  
            # Magnified displacements in global coordinates
            DX_plot = append(DX_plot, dx_tot, axis=0)
        
        # Sum the component displacements to obtain overall displacement
        D_plot = DY_plot + DZ_plot + DX_plot
      
        # Generate vtk points
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(len(D_plot))
      
        for i in range(len(D_plot)):
            points.SetPoint(i, D_plot[i, 0], D_plot[i, 1], D_plot[i, 2])
      
        # Generate vtk lines
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(len(D_plot))
      
        for i in range(len(D_plot)):
            lines.InsertCellPoint(i)
      
        # Create a polyline source from the defined points and lines
        self.source = vtk.vtkPolyData()
        self.source.SetPoints(points)
        self.source.SetLines(lines)   
                                  
class VisDeformedSpring():
    
    def __init__(self, spring, nodes, scale_factor, combo_name='Combo 1'):

        # Determine if this spring is active for each load combination
        self.active = spring.active
        
        # Generate a line source for the spring
        self.source = vtk.vtkLineSource()
        
        # Find the deformed position of the local i-node
        # Step through each node
        for node in nodes.values():
      
            # Check to see if the current node is the i-node
            if node.name == spring.i_node.name:
                Xi = node.X + node.DX[combo_name]*scale_factor
                Yi = node.Y + node.DY[combo_name]*scale_factor
                Zi = node.Z + node.DZ[combo_name]*scale_factor
                self.source.SetPoint1(Xi, Yi, Zi)
        
            # Check to see if the current node is the i-node
            if node.name == spring.j_node.name:
                Xj = node.X + node.DX[combo_name]*scale_factor
                Yj = node.Y + node.DY[combo_name]*scale_factor
                Zj = node.Z + node.DZ[combo_name]*scale_factor
                self.source.SetPoint2(Xj, Yj, Zj)
        
        self.source.Update()
    
class VisPtLoad():
    '''
    Creates a point load for the viewer
    '''
    
    def __init__(self, position, direction, length, label_text=None, annotation_size=5, color=None):
        '''
        Constructor.
      
        Parameters
        ----------
        position : tuple
          A tuple of X, Y and Z coordinates for the point of the load arrow: (X, Y, Z).
        direction : tuple
          A tuple indicating the direction vector for the load arrow: (i, j, k).
        length : number
          The length of the load arrow.
        tip_length : number
          The height of the arrow head.
        label_text : string
          Text that will show up at the tail of the arrow. If set to 'None' no text will be displayed.
        '''
      
        # Create a unit vector in the direction of the 'direction' vector
        unitVector = direction/norm(direction)
      
        # Create a 'vtkAppendPolyData' filter to append the tip and shaft together into a single dataset
        self.polydata = vtk.vtkAppendPolyData()
      
        # Determine if the load is positive or negative
        if length == 0:
            sign = 1
        else:
            sign = abs(length)/length
      
        # Generate the tip of the load arrow
        tip_length = abs(length)/4
        radius = abs(length)/16
        tip = vtk.vtkConeSource()
        tip.SetCenter(position[0] - tip_length*sign*0.5*unitVector[0], \
                      position[1] - tip_length*sign*0.5*unitVector[1], \
                      position[2] - tip_length*sign*0.5*unitVector[2])
        tip.SetDirection([direction[0]*sign, direction[1]*sign, direction[2]*sign])
        tip.SetHeight(tip_length)
        tip.SetRadius(radius)
        tip.Update()
      
        # Add the arrow tip to the append filter
        self.polydata.AddInputData(tip.GetOutput())
        
        # Create the shaft
        shaft = vtk.vtkLineSource()
        shaft.SetPoint1(position)
        shaft.SetPoint2((position[0]-length*unitVector[0], position[1]-length*unitVector[1], position[2]-length*unitVector[2]))
        shaft.Update()
      
        # Copy and append the shaft data to the append filter
        self.polydata.AddInputData(shaft.GetOutput())
        self.polydata.Update()
      
        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.polydata.GetOutputPort())
        self.actor = vtk.vtkActor()
        if color is None: self.actor.GetProperty().SetColor(0, 255, 0) # Green
        elif color == 'black': self.actor.GetProperty().SetColor(0, 0, 0)  # Black
        self.actor.SetMapper(mapper)
      
        # Create the label if needed
        if label_text != None:
      
            # Create the label and set its text
            self.label = vtk.vtkVectorText()
            self.label.SetText(label_text)
        
            # Set up a mapper for the label
            lblMapper = vtk.vtkPolyDataMapper()
            lblMapper.SetInputConnection(self.label.GetOutputPort())
        
            # Set up an actor for the label
            self.lblActor = vtk.vtkFollower()
            self.lblActor.SetMapper(lblMapper)
            self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
            self.lblActor.SetPosition(position[0] - (length - 0.6*annotation_size)*unitVector[0], \
                                      position[1] - (length - 0.6*annotation_size)*unitVector[1], \
                                      position[2] - (length - 0.6*annotation_size)*unitVector[2])
            if color is None: self.lblActor.GetProperty().SetColor(0, 255, 0)  # Green
            elif color == 'black': self.lblActor.GetProperty().SetColor(0, 0, 0)  # Black
      
class VisDistLoad():
    '''
    Creates a distributed load for the viewer
    '''
    
    def __init__(self, position1, position2, direction, length1, length2, label_text1, label_text2, annotation_size=5, color=None):
        '''
        Constructor.
        '''
      
        # Calculate the length of the distributed load
        loadLength = ((position2[0]-position1[0])**2 + (position2[1]-position1[1])**2 + (position2[2]-position1[2])**2)**0.5
      
        # Find the direction cosines for the line the load acts on
        lineDirCos = [(position2[0]-position1[0])/loadLength, (position2[1]-position1[1])/loadLength, (position2[2]-position1[2])/loadLength]
      
        # Find the direction cosines for the direction the load acts in
        dirDirCos = direction/norm(direction)
      
        # Create point loads at intervals roughly equal to 75% of the load's largest length (magnitude)
        # Add text labels to the first and last load arrow
        if loadLength > 0:
            num_steps = int(round(0.75*loadLength/max(abs(length1), abs(length2)), 0))
        else:
            num_steps = 0

        num_steps = max(num_steps, 1)
        step = loadLength/num_steps
        ptLoads = []
      
        for i in range(num_steps + 1):
      
            # Calculate the position (X, Y, Z) of this load arrow's point
            position = (position1[0] + i*step*lineDirCos[0], position1[1] + i*step*lineDirCos[1], position1[2] + i*step*lineDirCos[2])
        
            # Determine the length of this load arrow
            length = length1 + (length2 - length1)/loadLength*i*step
        
            # Determine the label's text
            if i == 0:
                label_text = label_text1
            elif i == num_steps:
                label_text = label_text2
        
            # Create the load arrow
            ptLoads.append(VisPtLoad(position, direction, length, label_text, annotation_size=annotation_size))
          
        # Draw a line between the first and last load arrow's tails
        tail_line = vtk.vtkLineSource()
        tail_line.SetPoint1((position1[0] - length1*dirDirCos[0], position1[1] - length1*dirDirCos[1], position1[2] - length1*dirDirCos[2]))
        tail_line.SetPoint2((position2[0] - length2*dirDirCos[0], position2[1] - length2*dirDirCos[1], position2[2] - length2*dirDirCos[2]))
      
        # Combine all the geometry into one 'vtkPolyData' object
        self.polydata = vtk.vtkAppendPolyData()
        for arrow in ptLoads:
            arrow.polydata.Update()
            self.polydata.AddInputData(arrow.polydata.GetOutput())
        
        tail_line.Update()
        self.polydata.AddInputData(tail_line.GetOutput())
        self.polydata.Update()
      
        # Create a mapper and actor for the geometry
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.polydata.GetOutputPort())
        self.actor = vtk.vtkActor()
        if color is None: self.actor.GetProperty().SetColor(0, 255, 0)  # Green
        elif color == 'black': self.actor.GetProperty().SetColor(0, 0, 0)  # Black
        self.actor.SetMapper(mapper)
      
        # Get the actors for the labels
        self.lblActors = [ptLoads[0].lblActor, ptLoads[len(ptLoads) - 1].lblActor]

class VisMoment():
    '''
    Creates a concentrated moment for the viewer
    '''
    
    def __init__(self, center, direction, radius, label_text=None, annotation_size=5, color=None):
        '''
        Constructor.
      
        Parameters
        ----------
        center : tuple
          A tuple of X, Y and Z coordinates for center of the moment: (X, Y, Z).
        direction : tuple
          A tuple indicating the direction vector for the moment: (i, j, k).
        radius : number
          The radius of the moment.
        tip_length : number
          The height of the arrow head.
        label_text : string
          Text that will show up at the tail of the moment. If set to 'None' no text will be displayed.
        '''
      
        # Create an append filter to store load polydata in
        self.polydata = vtk.vtkAppendPolyData()
        
        # Find a vector perpendicular to the directional unit vector
        v1 = direction/norm(direction)  # v1 = The directional unit vector for the moment
        v2 = _PerpVector(v1)             # v2 = A unit vector perpendicular to v1
        v3 = cross(v1, v2)
        v3 = v3/norm(v3)                # v3 = A unit vector perpendicular to v1 and v2
      
        # Generate an arc for the moment
        Xc, Yc, Zc = center
        arc = vtk.vtkArcSource()
        arc.SetCenter(Xc, Yc, Zc)
        arc.SetPoint1(Xc + v2[0]*radius, Yc + v2[1]*radius, Zc + v2[2]*radius)
        arc.SetPoint2(Xc + v3[0]*radius, Yc + v3[1]*radius, Zc + v3[2]*radius)
        arc.SetNegative(True)
        arc.SetResolution(20)
        arc.Update()
        self.polydata.AddInputData(arc.GetOutput())
      
        # Generate the arrow tip at the end of the arc
        tip_length = radius/2
        cone_radius = radius/8
        tip = vtk.vtkConeSource()
        tip.SetCenter(arc.GetPoint1()[0], arc.GetPoint1()[1], arc.GetPoint1()[2])
        tip.SetDirection(cross(v1, v2))
        tip.SetHeight(tip_length)
        tip.SetRadius(cone_radius)
        tip.Update()
        self.polydata.AddInputData(tip.GetOutput())
      
        # Update the polydata one last time now that we're done appending items to it
        self.polydata.Update()
      
        # Create the text label
        label = vtk.vtkVectorText()
        label.SetText(label_text)
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition(Xc + v3[0]*(radius + 0.25*annotation_size), \
                                  Yc + v3[1]*(radius + 0.25*annotation_size), \
                                  Zc + v3[2]*(radius + 0.25*annotation_size))
        if color is None: self.lblActor.GetProperty().SetColor(0, 255, 0)  # Green
        elif color == 'black': self.lblActor.GetProperty().SetColor(0, 0, 0)  # Black

class VisAreaLoad():
    '''
    Creates an area load for the viewer
    '''
    
    def __init__(self, position0, position1, position2, position3, direction, length, label_text, annotation_size=5, theme='default'):
        '''
        Constructor
        '''
      
        # Create a point load for each corner of the area load
        ptLoads = []
        ptLoads.append(VisPtLoad(position0, direction, length, label_text, annotation_size=annotation_size))
        ptLoads.append(VisPtLoad(position1, direction, length, label_text, annotation_size=annotation_size))
        ptLoads.append(VisPtLoad(position2, direction, length, label_text, annotation_size=annotation_size))
        ptLoads.append(VisPtLoad(position3, direction, length, label_text, annotation_size=annotation_size))
      
        # Find the direction cosines for the direction the load acts in
        dirDirCos = direction/norm(direction)
      
        # Find the positions of the tails of all the arrows at the corners of the area load. This is
        # where we will place the polygon.
        self.p0 = position0 - dirDirCos*length
        self.p1 = position1 - dirDirCos*length
        self.p2 = position2 - dirDirCos*length
        self.p3 = position3 - dirDirCos*length
      
        # Combine all geometry into one 'vtkPolyData' object
        self.polydata = vtk.vtkAppendPolyData()
        for arrow in ptLoads:
            self.polydata.AddInputData(arrow.polydata.GetOutput())
        self.polydata.Update()
      
        # Add a label
        self.label_actor = ptLoads[0].lblActor

        # Add color to the area load label
        if theme == 'print':
            self.label_actor.GetProperty().SetColor(255/255, 0/255, 0/255)  # red
        elif theme == 'default':
            self.label_actor.GetProperty().SetColor(0/255, 255/255, 0/255)  # green

def _PerpVector(v):
    '''
    Returns a unit vector perpendicular to v=[i, j, k]
    '''
    
    i = v[0]
    j = v[1]
    k = v[2]
    
    # Find a vector in a direction perpendicular to <i, j, k>
    if i == 0:
        i2 = 1
        j2 = 0
        k2 = 0
    elif j == 0:
        i2 = 0
        j2 = 1
        k2 = 0
    elif k == 0:
        i2 = 0
        j2 = 0
        k2 = 1
    else:
        i2 = 1
        j2 = 1
        k2 = -(i*i2+j*j2)/k
    
    # Return the unit vector
    return [i2, j2, k2]/norm([i2, j2, k2])

def _PrepContour(model, stress_type='Mx', combo_name='Combo 1'):

    if stress_type != None:
    
        # Erase any previous contours
        for node in model.Nodes.values():
            node.contour = []
    
        # Step through each element in the model
        for element in list(model.Quads.values()) + list(model.Plates.values()):
            
            # Rectangular elements and quadrilateral elements have different local coordinate systems.
            # Rectangles are based on a traditional (x, y) system, while quadrilaterals are based on a
            # 'natural' (r, s) coordinate system. To reduce duplication of code for both these elements
            # we'll define the edges of the plate here for either element using the (r, s) terminology.
            if element.type == 'Rect':
                r_left = 0
                r_right = element.width()
                s_bot = 0
                s_top = element.height()
            else:
                r_left = -1
                r_right = 1
                s_bot = -1
                s_top = 1
      
            # Determine which stress result has been requested by the user
            if stress_type == 'dz':
                # Internally PyNite defines the nodes for a rectangular element in the order (i, j, m, n),
                # while it defines the nodes for a quadrilateral element in the order (m, n, i, j)
                if element.type == 'Rect':
                    i, j, m, n = element.d(combo_name)[[2, 8, 14, 20], :]
                else:
                    i, j, m, n = element.d(combo_name)[[14, 20, 2, 8], :]
                element.i_node.contour.append(i)
                element.j_node.contour.append(j)
                element.m_node.contour.append(m)
                element.n_node.contour.append(n)
            elif stress_type == 'Mx':
                element.i_node.contour.append(element.moment(r_left, s_bot, combo_name)[0])
                element.j_node.contour.append(element.moment(r_right, s_bot, combo_name)[0])
                element.m_node.contour.append(element.moment(r_right, s_top, combo_name)[0])
                element.n_node.contour.append(element.moment(r_left, s_top, combo_name)[0])
            elif stress_type == 'My':
                element.i_node.contour.append(element.moment(r_left, s_bot, combo_name)[1])
                element.j_node.contour.append(element.moment(r_right, s_bot, combo_name)[1])
                element.m_node.contour.append(element.moment(r_right, s_top, combo_name)[1])
                element.n_node.contour.append(element.moment(r_left, s_top, combo_name)[1])
            elif stress_type == 'Mxy':
                element.i_node.contour.append(element.moment(r_left, s_bot, combo_name)[2])
                element.j_node.contour.append(element.moment(r_right, s_bot, combo_name)[2])
                element.m_node.contour.append(element.moment(r_right, s_top, combo_name)[2])
                element.n_node.contour.append(element.moment(r_left, s_top, combo_name)[2])
            elif stress_type == 'Qx':
                element.i_node.contour.append(element.shear(r_left, s_bot, combo_name)[0])
                element.j_node.contour.append(element.shear(r_right, s_bot, combo_name)[0])
                element.m_node.contour.append(element.shear(r_right, s_top, combo_name)[0])
                element.n_node.contour.append(element.shear(r_left, s_top, combo_name)[0])
            elif stress_type == 'Qy':
                element.i_node.contour.append(element.shear(r_left, s_bot, combo_name)[1])
                element.j_node.contour.append(element.shear(r_right, s_bot, combo_name)[1])
                element.m_node.contour.append(element.shear(r_right, s_top, combo_name)[1])
                element.n_node.contour.append(element.shear(r_left, s_top, combo_name)[1])
            elif stress_type == 'Sx':
                element.i_node.contour.append(element.membrane(r_left, s_bot, combo_name)[0])
                element.j_node.contour.append(element.membrane(r_right, s_bot, combo_name)[0])
                element.m_node.contour.append(element.membrane(r_right, s_top, combo_name)[0])
                element.n_node.contour.append(element.membrane(r_left, s_top, combo_name)[0])
            elif stress_type == 'Sy':
                element.i_node.contour.append(element.membrane(r_left, s_bot, combo_name)[1])
                element.j_node.contour.append(element.membrane(r_right, s_bot, combo_name)[1])
                element.m_node.contour.append(element.membrane(r_right, s_top, combo_name)[1])
                element.n_node.contour.append(element.membrane(r_left, s_top, combo_name)[1])
            elif stress_type == 'Txy':
                element.i_node.contour.append(element.membrane(r_left, s_bot, combo_name)[2])
                element.j_node.contour.append(element.membrane(r_right, s_bot, combo_name)[2])
                element.m_node.contour.append(element.membrane(r_right, s_top, combo_name)[2])
                element.n_node.contour.append(element.membrane(r_left, s_top, combo_name)[2])                 
      
        # Average the values at each node to obtain a smoothed contour
        for node in model.Nodes.values():
            # Prevent divide by zero errors for nodes with no contour values
            if node.contour != []:
                node.contour = sum(node.contour)/len(node.contour)

def _DeformedShape(model, vtk_renderer, scale_factor, annotation_size, combo_name, render_nodes=True, theme='default'):
    '''
    Renders the deformed shape of a model.
    
    Parameters
    ----------
    model : FEModel3D
        Finite element model to be rendered.
    renderer : vtk.vtkRenderer
        The VTK renderer object that will render the model.
    scale_factor : number
        The scale factor to apply to the model deformations.
    annotation_size : number
        Controls the height of text displayed with the model. The units used for `annotation_size` are
        the same as those used for lengths in the model. Sizes of other objects (such as nodes) are
        related to this value.
    combo_name : string
        The load case used for rendering the deflected shape.
    
    Returns
    -------
    None.
    '''
    
    # Create an append filter to add all the shape polydata to
    append_filter = vtk.vtkAppendPolyData()
    
    # Check if nodes are to be rendered
    if render_nodes == True:
        
        # Add the deformed nodes to the append filter
        for node in model.Nodes.values():
            
            vis_node = VisDeformedNode(node, scale_factor, annotation_size, combo_name)
            append_filter.AddInputData(vis_node.source.GetOutput())
        
    # Add the springs to the append filter
    for spring in model.Springs.values():
        
        # Only add the spring if it is active for the given load combination
        if spring.active[combo_name] == True:
            
            vis_spring = VisDeformedSpring(spring, model.Nodes, scale_factor, combo_name)
            append_filter.AddInputData(vis_spring.source.GetOutput())
            
    # Add the members to the append filter
    for member in model.Members.values():
        
        # Only add the member if it is active for the given load combination.
        if member.active[combo_name] == True:
            
            vis_member = VisDeformedMember(member, model.Nodes, scale_factor, combo_name)
            append_filter.AddInputData(vis_member.source)
            
    # Create a mapper and actor for the append filter
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(append_filter.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Adjust the color
    if theme == 'default':
        actor.GetProperty().SetColor(255/255, 255/255, 0/255)  # Yellow
    elif theme == 'print':
        actor.GetProperty().SetColor(26/255, 26/255, 26/255)  # Dark Grey
    
    # Add the actor to the renderer
    vtk_renderer.AddActor(actor)

def _RenderLoads(model, renderer, annotation_size, combo_name, case, theme='default'):

    # Create an append filter to store all the polydata in. This will allow us to use fewer actors to
    # display all the loads, which will greatly improve rendering speed as the user interacts. VTK
    # becomes very slow when a large number of actors are used.
    polydata = vtk.vtkAppendPolyData()
    
    # Polygons are treated as cells in VTK. Create a cell array to store all the area load polygons
    # in. We'll also create a list of points to store the polygon points in. The polydata for these
    # polygons will be stored separately from the other load data.
    polygons = vtk.vtkCellArray()
    polygon_points = vtk.vtkPoints()
    polygon_polydata = vtk.vtkPolyData()
    
    # Get the maximum load magnitudes that will be used to normalize the display scale
    max_pt_load, max_moment, max_dist_load, max_area_load = _MaxLoads(model, combo_name, case)
    
    # Display the requested load combination, or 'Combo 1' if no load combo or case has been
    # specified
    if case == None:
        # Store model.LoadCombos[combo].factors under a simpler name for use below
        load_factors = model.LoadCombos[combo_name].factors
    else:
        # Set up a load combination dictionary that represents the load case
        load_factors = {case: 1}
    
    # Step through each node
    for node in model.Nodes.values():
        
        # Step through and display each nodal load
        for load in node.NodeLoads:
            
            # Determine if this load is part of the requested LoadCombo or case
            if load[2] in load_factors:
              
                # Calculate the factored value for this load and it's sign (positive or negative)
                load_value = load[1]*load_factors[load[2]]
                if load_value != 0:
                    sign = load_value/abs(load_value)
                else:
                    sign = 1
                
                # Display the load
                if load[0] == 'FX':
                    ptLoad = VisPtLoad((node.X - 0.6*annotation_size*sign, node.Y, node.Z), [1, 0, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'FY':
                    ptLoad = VisPtLoad((node.X, node.Y - 0.6*annotation_size*sign, node.Z), [0, 1, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'FZ':
                    ptLoad = VisPtLoad((node.X, node.Y, node.Z - 0.6*annotation_size*sign), [0, 0, 1], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'MX':
                    ptLoad = VisMoment((node.X, node.Y, node.Z), (1*sign, 0, 0), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'MY':
                    ptLoad = VisMoment((node.X, node.Y, node.Z), (0, 1*sign, 0), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'MZ':
                    ptLoad = VisMoment((node.X, node.Y, node.Z), (0, 0, 1*sign), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                
                polydata.AddInputData(ptLoad.polydata.GetOutput())
                renderer.AddActor(ptLoad.lblActor)
                ptLoad.lblActor.SetCamera(renderer.GetActiveCamera())
                
    # Step through each member
    for member in model.Members.values():
    
        # Get the direction cosines for the member's local axes
        dir_cos = member.T()[0:3, 0:3]
      
        # Get the starting point for the member
        x_start, y_start, z_start = member.i_node.X, member.i_node.Y, member.i_node.Z
      
        # Step through each member point load
        for load in member.PtLoads:
      
            # Determine if this load is part of the requested load combination
            if load[3] in load_factors:
        
                # Calculate the factored value for this load and it's sign (positive or negative)
                load_value = load[1]*load_factors[load[3]]
                sign = load_value/abs(load_value)
          
                # Calculate the load's location in 3D space
                x = load[2]
                position = [x_start + dir_cos[0, 0]*x, y_start + dir_cos[0, 1]*x, z_start + dir_cos[0, 2]*x]
          
                # Display the load
                if load[0] == 'Fx':
                    ptLoad = VisPtLoad(position, dir_cos[0, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'Fy':
                    ptLoad = VisPtLoad(position, dir_cos[1, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'Fz':
                    ptLoad = VisPtLoad(position, dir_cos[2, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'Mx':
                    ptLoad = VisMoment(position, dir_cos[0, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'My':
                    ptLoad = VisMoment(position, dir_cos[1, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'Mz':
                    ptLoad = VisMoment(position, dir_cos[2, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'FX':
                    ptLoad = VisPtLoad(position, [1, 0, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'FY':
                    ptLoad = VisPtLoad(position, [0, 1, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'FZ':
                    ptLoad = VisPtLoad(position, [0, 0, 1], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'MX':
                    ptLoad = VisMoment(position, [1*sign, 0, 0], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'MY':
                    ptLoad = VisMoment(position, [0, 1*sign, 0], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                elif load[0] == 'MZ':
                    ptLoad = VisMoment(position, [0, 0, 1*sign], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
            
                polydata.AddInputData(ptLoad.polydata.GetOutput())
                renderer.AddActor(ptLoad.lblActor)
                ptLoad.lblActor.SetCamera(renderer.GetActiveCamera())
    
        # Step through each member distributed load
        for load in member.DistLoads:
    
            # Determine if this load is part of the requested load combination
            if load[5] in load_factors:
        
                # Calculate the factored value for this load and it's sign (positive or negative)
                w1 = load[1]*load_factors[load[5]]
                w2 = load[2]*load_factors[load[5]]
          
                # Calculate the loads location in 3D space
                x1 = load[3]
                x2 = load[4]
                position1 = [x_start + dir_cos[0, 0]*x1, y_start + dir_cos[0, 1]*x1, z_start + dir_cos[0, 2]*x1]
                position2 = [x_start + dir_cos[0, 0]*x2, y_start + dir_cos[0, 1]*x2, z_start + dir_cos[0, 2]*x2]
                
                # Display the load
                if load[0] == 'Fx':
                    distLoad = VisDistLoad(position1, position2, dir_cos[0, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'Fy':
                    distLoad = VisDistLoad(position1, position2, dir_cos[1, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'Fz':
                    distLoad = VisDistLoad(position1, position2, dir_cos[2, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'FX':
                    distLoad = VisDistLoad(position1, position2, [1, 0, 0], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'FY':
                    distLoad = VisDistLoad(position1, position2, [0, 1, 0], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'FZ':
                    distLoad = VisDistLoad(position1, position2, [0, 0, 1], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
               
                polydata.AddInputData(distLoad.polydata.GetOutput())
                renderer.AddActor(distLoad.lblActors[0])
                renderer.AddActor(distLoad.lblActors[1])
                distLoad.lblActors[0].SetCamera(renderer.GetActiveCamera())
                distLoad.lblActors[1].SetCamera(renderer.GetActiveCamera())   
    
    # Step through each plate
    i = 0
    for plate in list(model.Plates.values()) + list(model.Quads.values()):
    
        # Get the direction cosines for the plate's local z-axis
        dir_cos = plate.T()[0:3, 0:3]
        dir_cos = dir_cos[2]
      
        # Step through each plate load
        for load in plate.pressures:
      
            # Determine if this load is part of the requested load combination
            if load[1] in load_factors:
        
                # Calculate the factored value for this load
                load_value = load[0]*load_factors[load[1]]
                
                # Find the sign for this load. Intercept any divide by zero errors
                if load[0] == 0:
                    sign = 1
                else:
                    sign = abs(load[0])/load[0]
          
                # Find the position of the load's 4 corners
                position0 = [plate.i_node.X, plate.i_node.Y, plate.i_node.Z]
                position1 = [plate.j_node.X, plate.j_node.Y, plate.j_node.Z]
                position2 = [plate.m_node.X, plate.m_node.Y, plate.m_node.Z]
                position3 = [plate.n_node.X, plate.n_node.Y, plate.n_node.Z]
          
                # Create an area load and get its data
                area_load = VisAreaLoad(position0, position1, position2, position3, dir_cos*sign, abs(load_value)/max_area_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
          
                # Add the area load's arrows to the overall load polydata
                polydata.AddInputData(area_load.polydata.GetOutput())
          
                # Add the 4 points at the corners of this area load to the list of points
                polygon_points.InsertNextPoint(area_load.p0[0], area_load.p0[1], area_load.p0[2])
                polygon_points.InsertNextPoint(area_load.p1[0], area_load.p1[1], area_load.p1[2])
                polygon_points.InsertNextPoint(area_load.p2[0], area_load.p2[1], area_load.p2[2])
                polygon_points.InsertNextPoint(area_load.p3[0], area_load.p3[1], area_load.p3[2])
          
                # Create a polygon based on the four points we just defined.
                # The 1st number in `SetId()` is the local point id
                # The 2nd number in `SetId()` is the global point id
                polygon = vtk.vtkPolygon()
                polygon.GetPointIds().SetNumberOfIds(4)
                polygon.GetPointIds().SetId(0, i*4)
                polygon.GetPointIds().SetId(1, i*4 + 1)
                polygon.GetPointIds().SetId(2, i*4 + 2)
                polygon.GetPointIds().SetId(3, i*4 + 3)
          
                # Add the polygon to the list of polygons
                polygons.InsertNextCell(polygon)
                
                # Add the load label
                renderer.AddActor(area_load.label_actor)
          
                # Set the text to follow the camera as the user interacts
                area_load.label_actor.SetCamera(renderer.GetActiveCamera())
          
                # `i` keeps track of the next polygon's ID. We've just added a polygon, so `i` needs to
                # go up 1.
                i += 1
        
                # Create polygon polydata from all the points and polygons we just defined
                polygon_polydata.SetPoints(polygon_points)
                polygon_polydata.SetPolys(polygons)
    
    # Set up an actor and mapper for the loads
    load_mapper = vtk.vtkPolyDataMapper()
    load_mapper.SetInputConnection(polydata.GetOutputPort())
    load_actor = vtk.vtkActor()
    load_actor.SetMapper(load_mapper)

    # Colorize the loads
    if theme == 'default':
        load_actor.GetProperty().SetColor(0/255, 255/255, 0/255)  # Green
    elif theme == 'print':
        load_actor.GetProperty().SetColor(255/255, 0/255, 0/255)  # Red

    # Add the load actor to the renderer
    renderer.AddActor(load_actor)
    
    # Set up an actor and a mapper for the area load polygons
    polygon_mapper = vtk.vtkPolyDataMapper()
    polygon_mapper.SetInputData(polygon_polydata)
    polygon_actor = vtk.vtkActor()

    # polygon_actor.GetProperty().SetOpacity(0.5)      # 50% opacity
    polygon_actor.SetMapper(polygon_mapper)
    renderer.AddActor(polygon_actor)
    
    # Set the color of the area load polygons
    if theme == 'default':
        polygon_actor.GetProperty().SetColor(0/255, 255/255, 0/255)  # Green
    elif theme == 'print':
        polygon_actor.GetProperty().SetColor(255/255, 0/255, 0/255)  # Red

def _RenderContours(model, renderer, deformed_shape, deformed_scale, color_map,
                    scalar_bar, scalar_bar_text_size, combo_name, show_plate_edges = True,
                    num_contours = -1, contour_range=(None,None),
                    theme='default'):
  
    # Create a new `vtkCellArray` object to store the elements
    plates = vtk.vtkCellArray()
    
    # Create a `vtkPoints` object to store the coordinates of the corners of the elements
    plate_points = vtk.vtkPoints()
    
    # Create 2 lists to store plate result
    # `results` will store the results in a Python iterable list
    # `plate_results` will store the results in a `vtkDoubleArray` for VTK
    results = []
    plate_results = vtk.vtkDoubleArray()
    plate_results.SetNumberOfComponents(1)
    
    # Each element will be assigned a unique element number `i` beginning at 0
    i = 0
    
    # Calculate the smoothed contour results at each node
    _PrepContour(model, color_map, combo_name)
    
    # Add each plate and quad in the model to the cell array we just created
    for item in list(model.Plates.values()) + list(model.Quads.values()):
        
        # Create a point for each corner (must be in counter clockwise order)
        if deformed_shape == True:
            p0 = [item.i_node.X + item.i_node.DX[combo_name]*deformed_scale,
                  item.i_node.Y + item.i_node.DY[combo_name]*deformed_scale,
                  item.i_node.Z + item.i_node.DZ[combo_name]*deformed_scale]
            p1 = [item.j_node.X + item.j_node.DX[combo_name]*deformed_scale,
                  item.j_node.Y + item.j_node.DY[combo_name]*deformed_scale,
                  item.j_node.Z + item.j_node.DZ[combo_name]*deformed_scale]
            p2 = [item.m_node.X + item.m_node.DX[combo_name]*deformed_scale,
                  item.m_node.Y + item.m_node.DY[combo_name]*deformed_scale,
                  item.m_node.Z + item.m_node.DZ[combo_name]*deformed_scale]
            p3 = [item.n_node.X + item.n_node.DX[combo_name]*deformed_scale,
                  item.n_node.Y + item.n_node.DY[combo_name]*deformed_scale,
                  item.n_node.Z + item.n_node.DZ[combo_name]*deformed_scale]
        else:
            p0 = [item.i_node.X, item.i_node.Y, item.i_node.Z]
            p1 = [item.j_node.X, item.j_node.Y, item.j_node.Z]
            p2 = [item.m_node.X, item.m_node.Y, item.m_node.Z]
            p3 = [item.n_node.X, item.n_node.Y, item.n_node.Z]
      
        # Add the points to the `vtkPoints` object we created earlier
        plate_points.InsertNextPoint(p0)
        plate_points.InsertNextPoint(p1)
        plate_points.InsertNextPoint(p2)
        plate_points.InsertNextPoint(p3)
      
        # Create a `vtkQuad` based on the four points we just defined
        # The 1st number in `SetId()` is the local point id
        # The 2nd number in `SetId()` is the global point id
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0, i*4)
        quad.GetPointIds().SetId(1, i*4 + 1)
        quad.GetPointIds().SetId(2, i*4 + 2)
        quad.GetPointIds().SetId(3, i*4 + 3)
      
        # Get the contour value for each node
        r0 = item.i_node.contour
        r1 = item.j_node.contour
        r2 = item.m_node.contour
        r3 = item.n_node.contour
            
        if color_map != None:
              
            # Save the results to the Python list of results we created earlier
            results.append(r0)
            results.append(r1)
            results.append(r2)
            results.append(r3)
                  
            # Save the results to the `vtkDoubleArray` list of results for VTK
            plate_results.InsertNextTuple([r0])
            plate_results.InsertNextTuple([r1])
            plate_results.InsertNextTuple([r2])
            plate_results.InsertNextTuple([r3])
      
        # Insert the quad into the cell array
        plates.InsertNextCell(quad)
      
        # Increment `i` for the next plate
        i += 1

    # Create a `vtkPolyData` object to store plate data in
    plate_polydata = vtk.vtkPolyData()
      
    # Add the points and plates to the dataset
    plate_polydata.SetPoints(plate_points)
    plate_polydata.SetPolys(plates)

    #Extract plate edges
    if show_plate_edges:
        polyedges = vtk.vtkExtractEdges()
        polyedges.SetInputData(plate_polydata)
        polyedges.Update()
        
        polyedgeMapper = vtk.vtkPolyDataMapper()
        polyedgeMapper.SetInputData(polyedges.GetOutput())
        
        polyedgeActor = vtk.vtkActor()
        polyedgeActor.SetMapper(polyedgeMapper)
        polyedgeActor.GetProperty().SetColor(0, 0, 0)
                
    #Add scalar data
    if color_map: plate_polydata.GetPointData().SetScalars(plate_results)
    
    # Map the results to the plates
    if color_map is None:
        # Setup actor and mapper for the plates
        plate_mapper = vtk.vtkPolyDataMapper()
        plate_mapper.SetInputData(plate_polydata)
                    
    else:
        #Convert mesh to triangles and subdivide for smoother results
        tri_converter = vtk.vtkTriangleFilter()
        tri_converter.SetInputDataObject(plate_polydata)
        tri_converter.Update()
        
        subdivisionFilter = vtk.vtkLinearSubdivisionFilter()
        subdivisionFilter.SetNumberOfSubdivisions(3)
        subdivisionFilter.SetInputData(tri_converter.GetOutput())
        subdivisionFilter.Update()
    
        min_contour = contour_range[0] if contour_range[0] is not None else min(results)
        max_contour = contour_range[1] if contour_range[1] is not None else max(results)
            
        if num_contours<0: 
            #Create smooth contours (256 values)
            
            # Create a `vtkLookupTable` for the colors used to map results
            lut = vtk.vtkLookupTable()
            lut.SetTableRange(min_contour, max_contour)
            lut.SetNumberOfColors(256) 
            # The commented code below can be uncommented and modified to change the color scheme
            # ctf = vtk.vtkColorTransferFunction()
            # ctf.SetColorSpaceToDiverging()
            # ctf.AddRGBPoint(min_contour, 255, 0, 255)  # Purple
            # ctf.AddRGBPoint(max_contour, 255, 0, 0)    # Red
            # for i in range(256):
            #   rgb = list(ctf.GetColor(float(i)/256))
            #   rgb.append(1.0)
            #   lut.SetTableValue(i, *rgb)
            lut.Build()
            
            # Setup actor and mapper for the plates
            plate_mapper = vtk.vtkPolyDataMapper()
            plate_mapper.SetInputData(subdivisionFilter.GetOutput())
            plate_mapper.SetLookupTable(lut)
            plate_mapper.SetUseLookupTableScalarRange(True)
            plate_mapper.SetScalarModeToUsePointData()
            
        else:
            #Create banded contours
            if num_contours < 2:
                raise ValueError("At least two contours are required to generate banded contours.")
            
            bf = vtk.vtkBandedPolyDataContourFilter()
            bf.SetInputData(subdivisionFilter.GetOutput())
            bf.SetGenerateContourEdges(True)
            bf.SetScalarModeToValue()
            bf.GenerateValues(num_contours, (min_contour, max_contour))
            bf.Update()
            
            lut = vtk.vtkLookupTable()
            lut.SetNumberOfTableValues(num_contours)
            lut.SetTableRange((min_contour,max_contour))
            lut.SetHueRange(.667, 0.0)
            lut.SetAlphaRange(1., 1.)
            lut.IndexedLookupOff()
            lut.Build()
            
            edgeMapper = vtk.vtkPolyDataMapper()
            edgeMapper.SetInputData(bf.GetContourEdgesOutput())
            edgeMapper.SetResolveCoincidentTopologyToPolygonOffset()
            
            edgeActor = vtk.vtkActor()
            edgeActor.SetMapper(edgeMapper)
            edgeActor.GetProperty().SetColor(0, 0, 0)
            renderer.AddActor(edgeActor)
            
            plate_mapper = vtk.vtkPolyDataMapper()
            plate_mapper.SetInputConnection(bf.GetOutputPort())
            plate_mapper.SetLookupTable(lut)
            plate_mapper.SetScalarRange((min_contour,max_contour))
            plate_mapper.SetScalarModeToUseCellData()
            
            plate_mapper.Update() 
      
        # Add the scalar bar for the contours.
        if scalar_bar:
            
            if Renderer.scalar == None:
                Renderer.scalar = vtk.vtkScalarBarActor()
            
            scalar = Renderer.scalar

            # This next group of lines controls the font on the scalar bar
            scalar.SetUnconstrainedFontSize(1)
            scalar_text = vtk.vtkTextProperty()
            scalar_text.SetFontSize(max(int(scalar_bar_text_size), 1))
            scalar_text.SetBold(1)

            # The `vtkTextProperty` object is white by default
            if theme == 'print':
                scalar_text.SetColor(255/255, 255/255, 255/255)  # Black
            
            scalar.SetLabelTextProperty(scalar_text)

            scalar.SetMaximumWidthInPixels(100)

            scalar.SetTextPositionToPrecedeScalarBar()

            scalar.SetLookupTable(lut)

            renderer.AddActor(scalar)

    # Add the actor for the plates to the renderer    
    plate_actor = vtk.vtkActor()
    plate_actor.SetMapper(plate_mapper)
    renderer.AddActor(plate_actor) 
        
    if show_plate_edges: 
        #Add edges last so that they're not obscured by the plate actor
        renderer.AddActor(polyedgeActor)      

def _MaxLoads(model, combo_name=None, case=None):

    max_pt_load = 0
    max_moment = 0
    max_dist_load = 0
    max_area_load = 0
    
    # Find the requested load combination or load case
    if case == None:
    
        # Step through each node
        for node in model.Nodes.values():
      
            # Step through each nodal load to find the largest one
            for load in node.NodeLoads:
              
                # Find the largest loads in the load combination
                if load[2] in model.LoadCombos[combo_name].factors:
                    if load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ':
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[2]]) > max_pt_load:
                            max_pt_load = abs(load[1]*model.LoadCombos[combo_name].factors[load[2]])
                    else:
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[2]]) > max_moment:
                            max_moment = abs(load[1]*model.LoadCombos[combo_name].factors[load[2]])
    
        # Step through each member
        for member in model.Members.values():
      
            # Step through each member point load
            for load in member.PtLoads:
                
                # Find and store the largest point load and moment in the load combination
                if load[3] in model.LoadCombos[combo_name].factors:
          
                    if (load[0] == 'Fx' or load[0] == 'Fy' or load[0] == 'Fz'
                    or  load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ'):
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[3]]) > max_pt_load:
                            max_pt_load = abs(load[1]*model.LoadCombos[combo_name].factors[load[3]])
                    else:
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[3]]) > max_moment:
                            max_moment = abs(load[1]*model.LoadCombos[combo_name].factors[load[3]])
      
            # Step through each member distributed load
            for load in member.DistLoads:
        
                #Find and store the largest distributed load in the load combination
                if load[5] in model.LoadCombos[combo_name].factors:
          
                    if abs(load[1]*model.LoadCombos[combo_name].factors[load[5]]) > max_dist_load:
                        max_dist_load = abs(load[1]*model.LoadCombos[combo_name].factors[load[5]])
                    if abs(load[2]*model.LoadCombos[combo_name].factors[load[5]]) > max_dist_load:
                        max_dist_load = abs(load[2]*model.LoadCombos[combo_name].factors[load[5]])
      
        # Step through each plate
        for plate in model.Plates.values():
      
            # Step through each plate load
            for load in plate.pressures:
        
                if load[1] in model.LoadCombos[combo_name].factors:
                    if abs(load[0]*model.LoadCombos[combo_name].factors[load[1]]) > max_area_load:
                        max_area_load = abs(load[0]*model.LoadCombos[combo_name].factors[load[1]])
      
        # Step through each quad
        for quad in model.Quads.values():
      
            # Step through each plate load
            for load in quad.pressures:
        
                # Check to see if the load case is in the requested load combination
                if load[1] in model.LoadCombos[combo_name].factors:
                    if abs(load[0]*model.LoadCombos[combo_name].factors[load[1]]) > max_area_load:
                        max_area_load = abs(load[0]*model.LoadCombos[combo_name].factors[load[1]])
      
    # Behavior if case has been specified
    else:
        
        # Step through each node
        for node in model.Nodes.values():
      
            # Step through each nodal load to find the largest one
            for load in node.NodeLoads:
              
                # Find the largest loads in the load case
                if load[2] == case:
                    if load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ':
                        if abs(load[1]) > max_pt_load:
                            max_pt_load = abs(load[1])
                    else:
                        if abs(load[1]) > max_moment:
                            max_moment = abs(load[1])
      
        # Step through each member
        for member in model.Members.values():
      
            # Step through each member point load
            for load in member.PtLoads:
              
                # Find and store the largest point load and moment in the load case
                if load[3] == case:
            
                    if (load[0] == 'Fx' or load[0] == 'Fy' or load[0] == 'Fz'
                    or  load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ'):
                        if abs(load[1]) > max_pt_load:
                            max_pt_load = abs(load[1])
                    else:
                        if abs(load[1]) > max_moment:
                            max_moment = abs(load[1])
        
            # Step through each member distributed load
            for load in member.DistLoads:
          
                # Find and store the largest distributed load in the load case
                if load[5] == case:
            
                    if abs(load[1]) > max_dist_load:
                        max_dist_load = abs(load[1])
                    if abs(load[2]) > max_dist_load:
                        max_dist_load = abs(load[2])
            
            # Step through each plate
            for plate in model.Plates.values():
          
                # Step through each plate load
                for load in plate.pressures:
          
                    if load[1] == case:
            
                        if abs(load[0]) > max_area_load:
                            max_area_load = abs(load[0])
          
        # Step through each quad
        for quad in model.Quads.values():
      
            # Step through each plate load
            for load in quad.pressures:
        
                if load[1] == case:
          
                    if abs(load[0]) > max_area_load:
                        max_area_load = abs(load[0])
        
    # Return the maximum loads in the load combination or load case
    return max_pt_load, max_moment, max_dist_load, max_area_load
