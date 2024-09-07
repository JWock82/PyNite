"""Module to import various solvers for PyNite
The GPU solver method is preferred, but requires PyTorch and a CUDA capable GPU, and is not available on all platforms.

The CPU solver is the next best option, and is available on all platforms via numpy and scipy
"""
import os
try:
    if os.environ.get('PYNITE_GPU',None) != 'True':
        raise ImportError(f'PYNITE_GPU environment variable not set to `True`')
    
    import torch
    import numpy
    
    if not torch.cuda.is_available():
        raise ImportError(f'CUDA not available')
    
    def solve(a:numpy.ndarray,b:numpy.ndarray)->numpy.ndarray:
        device = torch.device("cuda")
        
        a = torch.from_numpy(a).cfloat().to(device)
        b = torch.from_numpy(b).cfloat().to(device)

        res = torch.linalg.solve(a, b)

        return res.cpu().numpy()
    
    def spsolve(a:numpy.ndarray,b:numpy.ndarray)->numpy.ndarray:
        device = torch.device("cuda")
        a = a.todense() #FIXME: sparse solver is currently not in torch standard library, investigate other options
        a = torch.from_numpy(a).cfloat().to(device)
        b = torch.from_numpy(b).cfloat().to(device)

        res = torch.linalg.solve(a, b)

        return res.cpu().numpy()  
    
    print(f'PyNite Running with GPU solver: {torch.cuda.get_device_name()}')

except Exception as e:
    print(f'GPU solver not available: {e}')
    from numpy.linalg import solve
    from scipy.sparse.linalg import spsolve
