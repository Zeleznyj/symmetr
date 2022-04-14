import json
import numpy as np
from symmetr import parse, sym_res
from symmetr.tensors import Tensor

def test_jE():
    opt = parse('res j E -f MnTe/findsym.in' )
    X = sym_res(opt)

    with open('MnTe/vv.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_svE():
    opt = parse('res s.v E -f MnTe/findsym.in' )
    X = sym_res(opt)

    with open('MnTe/sv.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,3)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_jE2():
    opt = parse('res j E -f MnTe/MnTe_inp1.in --exp 2 --syms 2,5,7' )
    X = sym_res(opt)

    with open('MnTe/cond2_rot.json') as f:
        res = json.load(f)
    Xr = Tensor.load(res)
    X_zero = Tensor(0,3,2)

    assert X - Xr == X_zero
