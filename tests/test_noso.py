import json
import numpy as np
from symmetr import parse, sym_res
from symmetr.tensors import Tensor
from symmetr.symT import get_syms_noso

def test_svE_IrMn3():
    opt = parse('res s.v E -f IrMn3/findsym.in --noso' )
    X = sym_res(opt)

    with open('IrMn3/svE_noso.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,3)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_sE_Mn2Au():
    opt = parse('res s E -f Mn2Au/findsym.in --noso' )
    X = sym_res(opt)

    with open('Mn2Au/sE_noso.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_svE_NiF2():
    opt = parse('res s.v E -f NiF2/findsym.in --noso' )
    X = sym_res(opt)

    with open('NiF2/svE_noso.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,3)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_sE_080():
    opt = parse('res s E -f 0.80/findsym.in --noso' )
    X = sym_res(opt)

    with open('0.80/sE_noso.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_sE_080_p1():
    opt = parse('res s E -f 0.80/findsym.in --noso -p 1' )
    X = sym_res(opt)

    with open('0.80/sE_1_noso.out') as f:
        res = json.load(f)
    Xr0 = Tensor.load(res['X0'])
    Xr1 = Tensor.load(res['X1'])
    X_zero = Tensor(0,3,2)

    assert X[0] - Xr0 == X_zero
    assert X[1] - Xr1 == X_zero

def test_nsyms_IrMn3():
    syms = get_syms_noso(parse('res j E -f IrMn3/findsym.in'))
    assert len(syms) == 96

def test_nsyms_NiF2():
    syms = get_syms_noso(parse('res j E -f NiF2/findsym.in'))
    assert len(syms) == 36

def test_nsyms_NiF2():
    syms = get_syms_noso(parse('res j E -f 0.80/findsym.in'))
    assert len(syms) == 32
