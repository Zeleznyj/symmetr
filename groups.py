import re
import sys


def group_sym(name,dirname='.'):
    """Returns symmetry operations for a given group.

    Args:
        name: Name of the point group. Must be one of the magnetic point groups as given in the Sympy tables.
                Both number and name are accepted.
        dirname: Directory in which group tables are located.
                   Requires files syms_table.dat, syms_table_hex.dat and mag_groups.dat

    Returns:
        hex_group: Boolean. True if hexagonal or trigonal group, False otherwise. Useful for basis transformations.
        syms: The symmetrey operations in the format that is used in the code. The same format is outputted by the read.py module.
    """

    with open(dirname+'/syms_table.dat') as f:
        lines = f.readlines()
        syms_t = {}
        for line in lines:
            syms_t[line.split()[0]] = (line.split()[1],line.split()[2])

    with open(dirname+'/syms_table_hex.dat') as f:
        lines = f.readlines()
        syms_t_hex = {}
        for line in lines:
            syms_t_hex[line.split()[0]] = (line.split()[1],line.split()[2])

    hex_syms = ['6z','3z','3z-1','6z-1','21','22','23','-6z','-3z','-3z-1','-6z-1','m1','m2','m3']

    with open(dirname+'/mag_groups.txt') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.split(',',1)
        name_text = line[0].split()[2]
        name_num = line[0].split()[1]
        if name == name_text or name == name_num:
            line_f = line
            ops = [x[2:-1] for x in re.findall('\'\([12346xyzm\-]+\|',line[1])]
            ops_T = [x[2:-1] for x in re.findall('\"\([12346xyzm\-]+\|',line[1])]
            hex_group = False
            for op in ops + ops_T:
                if op in hex_syms:
                    hex_group = True

            sym_a=(name_text,name_num,hex_group,ops,ops_T)

            break

    syms=[]

    n = len(sym_a[3])
    for j in range(n):

        op = sym_a[3][j]

        if sym_a[2]:
            ops = syms_t_hex[op][0]
            opm = syms_t_hex[op][1]
        else:
            ops = syms_t[op][0]
            opm = syms_t[op][1]

        ops = ops.split(',')
        opm = opm.split(',')

        syms.append([j,ops,opm,'+1'])

    for j in range(len(ops_T)):

        op = sym_a[3][j]

        if sym_a[2]:

            ops = syms_t_hex['T'+op][0]
            opm = syms_t_hex['T'+op][1]
            ops = ops.split(',')
            opm = opm.split(',')

        else:

            ops = syms_t[op][0]
            opm = syms_t[op][1]
            ops = ops.split(',')
            opm = opm.split(',')

            for i in range(3):

                m = opm[i]

                if m == 'mx':
                    opm[i] = '-mx'
                if m == 'my':
                    opm[i] = '-my'
                if m == 'mz':
                    opm[i] = '-mz'
                if m == '-mx':
                    opm[i] = 'mx'
                if m == '-my':
                    opm[i] = 'my'
                if m == '-mz':
                    opm[i] = 'mz'

        syms.append([n+j,ops,opm,'-1'])

    return hex_group,syms


