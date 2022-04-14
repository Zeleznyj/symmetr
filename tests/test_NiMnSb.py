import json
from symmetr import parse, sym_res
from symmetr.tensors import Tensor

def test_sE_exp1():
    opt = parse('res s E --exp 1 -f NiMnSb/NiMnSb.in_nonmag ')
    X = sym_res(opt)
    with open('NiMnSb/sE_exp_1.out') as f:
        res = json.load(f)
    Xr = Tensor.load(res['X'])
    X_zero = Tensor(0,3,2)

    assert X - Xr == X_zero

def test_sE_exp2():
    opt = parse('res s E --exp 2 -f NiMnSb/NiMnSb.in_nonmag ')
    X = sym_res(opt)
    with open('NiMnSb/sE_exp_2.out') as f:
        res = json.load(f)
    Xr = Tensor.load(res['X'])
    X_zero = Tensor(0,3,2)

    assert X - Xr == X_zero
