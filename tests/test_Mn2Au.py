import json
import numpy as np
from symmetr import parse, sym_res
from symmetr.tensors import Tensor

def test_sE_p1():
    opt = parse('res s E -f Mn2Au/findsym.in -p 1')
    X = sym_res(opt)
    with open('Mn2Au/sv0.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_sE_p1_p2():
    opt = parse('res s E -f Mn2Au/findsym.in -p 1 -p2 2')
    X = sym_res(opt)
    with open('Mn2Au/sv_0_1.out') as f:
        res = json.load(f)
    Xr00 = Tensor.load(res['X00'])
    Xr01 = Tensor.load(res['X01'])
    Xr10 = Tensor.load(res['X10'])
    Xr11 = Tensor.load(res['X11'])

    X_zero = Tensor(0,3,2)

    assert X[0][0] - Xr00 == X_zero
    assert X[0][1] - Xr01 == X_zero
    assert X[1][0] - Xr10 == X_zero
    assert X[1][1] - Xr11 == X_zero

def test_jE():
    opt = parse('res v E -f Mn2Au/findsym.in  ')
    X = sym_res(opt)
    with open('Mn2Au/vv.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_sv_1_equiv():

    opt = parse('res s E -f Mn2Au/findsym.in -p 1 -e  ')
    X = sym_res(opt)

    confs = {}
    for i in X[1].confs:
        confs[i] = {}
        for j in X[1].confs[i]:
            confs[i][j] = np.array(X[1].confs[i][j],dtype=float).squeeze()
            
    with open('Mn2Au/sv_0_equiv.out') as f:
        res = json.load(f)

    X_zero = Tensor(0,3,2)
    for i in X[1].Xs:
        for j in range(2):
            assert X[1].Xs[i][j] - Tensor.load(res['Xs'][str(i)][j]) == X_zero

    for i in confs:
        for j in confs[i]:
            resc = np.array(res['confs'][str(i)][str(j)])
            assert np.linalg.norm(confs[i][j] - resc) < 1e-10 


