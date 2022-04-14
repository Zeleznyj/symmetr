import json
from symmetr import parse, sym_res
from symmetr.tensors import Tensor

def test_jE():
    opt = parse('res j E -g P-43m --ignore-same-op-sym')
    X = sym_res(opt)
    with open('groups/jE_-43m.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_sE_exp():

    opt = parse('res s E -g P4mm --exp 1 ')
    X = sym_res(opt)
    with open('groups/sE_4mm_exp1.out') as f:
        res = json.load(f)
    Xr = Tensor.load(res['X'])
    X_zero = Tensor(0,3,2)

    assert X- Xr == X_zero

