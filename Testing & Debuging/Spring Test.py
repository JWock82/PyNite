# Example of a spring assembly
# Units used in this example are inches and pounds

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

system = FEModel3D()

system.AddNode('1', 0, 0, 0)
system.AddNode('2', 30, 0, 0)
system.AddNode('3', 10, 0, 0)
system.AddNode('4', 20, 0, 0)

system.AddSpring('S1', '1', '3', 1000, False, False)
system.AddSpring('S2', '3', '4', 2000, False, False)
system.AddSpring('S3', '4', '2', 3000, False, False)

system.DefineSupport('1', True, True, True, True, True, True)
system.DefineSupport('2', True, True, True, True, True, True)

system.DefineSupport('3', False, True, True, True, True, True)
system.DefineSupport('4', False, True, True, True, True, True)

system.AddNodeLoad('4', 'FX', 5000)

system.Analyze(True)

print(system.GetNode('3').DX['Combo 1'])
print(system.GetNode('4').DX['Combo 1'])