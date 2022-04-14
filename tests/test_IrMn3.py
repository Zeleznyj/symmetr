import json
import numpy as np
from symmetr import parse, sym_res
from symmetr.tensors import Tensor

def test_jE():
    opt = parse('res j E -f IrMn3/findsym.in  ' )
    X = sym_res(opt)

    with open('IrMn3/jE.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_svE():
    opt = parse('res s.j E -f IrMn3/findsym.in' )
    X = sym_res(opt)

    with open('IrMn3/svE.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,3)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_jE_equiv():

    opt = parse('res j E -f IrMn3/findsym.in -e ')
    X = sym_res(opt)

    confs = {}
    for i in X[1].confs:
        confs[i] = {}
        for j in X[1].confs[i]:
            confs[i][j] = np.array(X[1].confs[i][j],dtype=float).squeeze()
            
    with open('IrMn3/jE_equiv.out') as f:
        res = json.load(f)

    X_zero = Tensor(0,3,2)
    for i in X[1].Xs:
        for j in range(2):
            assert X[1].Xs[i][j] - Tensor.load(res['Xs'][str(i)][j]) == X_zero

    for i in confs:
        for j in confs[i]:
            resc = np.array(res['confs'][str(i)][str(j)])
            assert np.linalg.norm(confs[i][j] - resc) < 1e-10

