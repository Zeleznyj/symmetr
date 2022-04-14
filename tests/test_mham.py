import json
import numpy as np
from symmetr import parse, sym_mham
from symmetr.tensors import Tensor

def test_mham23():
    opt = parse('mham -s 2,3 -f IrMn3/findsym.in' )
    X = sym_mham(opt)

    with open('IrMn3/mham23.out') as f:
        res = json.load(f)
    Xr = Tensor.load(res['X'])
    X_zero = Tensor(0,3,2)

    assert X - Xr == X_zero
