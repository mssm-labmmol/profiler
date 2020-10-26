import numpy as np

def fastCross(v1, v2):
    dims = len(v1.shape)
    eijk = np.zeros((3, 3, 3))
    eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
    eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1
    if dims == 2:
        return np.einsum('ijk,aj,ak->ai', eijk, v1, v2)
    elif dims == 1:
        out = np.einsum('ijk,aj,ak->ai', eijk, [v1], [v2])
        return out[0]
    else:
        raise Exception

            
