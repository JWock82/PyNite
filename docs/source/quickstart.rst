==========
Quickstart
==========

Here's a simple example of how to analyze a simple beam. Many more examples are available here:
`Examples <https://github.com/JWock82/PyNite/tree/main/Examples>`_.

.. code-block:: python

    # Example of a simply supported beam with a uniform distributed load.
    # Units used in this example are inches and kips
    # This example does not use load combinations. The program will create a
    # default load combindation called 'Combo 1'

    # Import `FEModel3D` from `PyNite`
    from PyNite import FEModel3D

    # Create a new finite element model
    beam = FEModel3D()

    # Add nodes (14 ft = 168 inches apart)
    beam.add_node('N1', 0, 0, 0)
    beam.add_node('N2', 168, 0, 0)

    # Define a material
    E = 29000       # Modulus of elasticity (ksi)
    G = 11200       # Shear modulus of elasticity (ksi)
    nu = 0.3        # Poisson's ratio
    rho = 2.836e-4  # Density (kci)
    beam.add_material('Steel', E, G, nu, rho)

    # Add a section with the following properties:
    # Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 20 in^2
    my_model.add_section('MySection', 20, 100, 150, 250)

    #Add a member
    beam.add_member('M1', 'N1', 'N2', 'Steel', 'MySection')

    # Provide simple supports
    beam.def_support('N1', True, True, True, False, False, False)
    beam.def_support('N2', True, True, True, True, False, False)

    # Add a uniform load of 200 lbs/ft to the beam (from 0 in to 168 in)
    beam.add_member_dist_load('M1', 'Fy', -200/1000/12, -200/1000/12, 0, 168)

    # Alternatively the following line would do apply the load to the full
    # length of the member as well
    # beam.add_member_dist_load('M1', 'Fy', 200/1000/12, 200/1000/12)

    # Analyze the beam
    beam.analyze()

    # Print the shear, moment, and deflection diagrams
    beam.members['M1'].plot_shear('Fy')
    beam.members['M1'].plot_moment('Mz')
    beam.members['M1'].plot_deflection('dy')

    # Print reactions at each end of the beam
    print('Left Support Reaction:', beam.nodes['N1'].RxnFY, 'kip')
    print('Right Support Reacton:', beam.nodes['N2'].RxnFY, 'kip')

    # Render the deformed shape of the beam magnified 100 times, with a text
    # height of 5 inches
    from PyNite.Visualization import Renderer
    renderer = Renderer(beam)
    renderer.annotation_size = 6
    renderer.deformed_shape = True
    renderer.deformed_scale = 100
    renderer.render_loads = True
    renderer.render_model()
