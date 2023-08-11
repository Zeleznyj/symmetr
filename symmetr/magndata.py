# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from six.moves import urllib
from lxml import html
from bs4 import BeautifulSoup
import pandas as pd
import re
from pymatgen.core import Structure, Lattice, Element

def download_html(id):
    page = urllib.request.urlopen('http://webbdcrista1.ehu.es/magndata/index.php?index='+id)
    mybytes = page.read()
    page_txt = mybytes.decode("utf-8", errors='ignore')
    page.close()
    return page_txt

def get_magndata_structure(id):
    try:
        htmlstr = download_html(id)
        root = html.fromstring(htmlstr)
        dfs = pd.read_html(process_html(htmlstr))

    except Exception as e:
        raise Exception("Error in downloading magndata page {}".format(e))

    info = {}
    # Parse html contents
    lattice_parameter = parse_lattice_parameter(root)
    lattice = Lattice.from_parameters(*lattice_parameter)

    # Fetch tabled contents
    sites, tabled_info = parse_tables(dfs)
    species = sites['species']
    coords = sites['coords']
    moments = sites['moments']
    occupancy = sites['occupancy']
    constraints = sites['constraints']
    labels = sites['labels']

    properties = {'magmom': moments, 'magmom_basis': ['crystal_scaled'] * len(moments), 'occupancy': occupancy,
                  'label': labels, 'magnetic_constraints': constraints}
    try:
        structure = Structure(lattice, species, coords, site_properties=properties)
    except Exception as e:
        raise Exception("Error in extracting the structure from MAGNDATA: {}".format(e))
    return structure


def parse_lattice_parameter(root):
    body_xpath = root.xpath('/html/body')
    body_lines = body_xpath[0].text_content().split('\n')

    for i, line in enumerate(body_lines):
        if re.search('Lattice parameters of the magnetic unit cell', line):
            lattice_line = body_lines[i + 1]
            break

    try:
        lstr = re.findall(r'([0-9]*\.[0-9]*)', lattice_line)
        lattice_parameter = [float(i) for i in lstr]
    except Exception as e:
       raise Exception("Problem in extracting lattice parameters from MAGNDATA: {}".format(e))

    return lattice_parameter

def parse_tables(dfs):
    # Extract atomic tables
    tables = []
    for df in dfs:
        column = df.columns.map(str)
        if 'Label' in column:
            tables.append(df)
        if 'Atom' in column:
            # tables.append(df[~df['Atom'].str.contains('click')])
            delete = []
            for i, row in df.iterrows():
                rowstr = row.Atom
                #if type(rowstr) is str or type(rowstr) is unicode:
                if type(rowstr) is str:
                    if not rowstr.isdigit():
                        delete.append(i)
            tables.append(df.drop(index=delete))

    # Remove duplication
    for i, table in enumerate(tables):
        column = table.columns.map(str)
        if '|M|' in column:
            start_index = i
            break
    for i in range(start_index + 1, len(tables)):
        column = tables[i].columns.map(str)
        if '|M|' in column:
            end_index = i
            break
    atom_tables = [tables[i] for i in range(start_index,end_index)]

    # Get infomations
    labels = []
    atom_types = []
    mag_atoms = []
    multiplicities = []
    coord_tables = []
    dis_atoms = []
    dis_labels = []
    dis_occupancy = {}
    dis_coord = []
    for table in atom_tables:
        column = table.columns.map(str)
        if 'Label' in column:
            labels.extend(table.loc[:,'Label'])
            atom_types.extend(table.loc[:, 'Atom type'])
            atom_types = strip_atom(atom_types)
            multiplicities.extend(table.loc[:,'Multiplicity'])
            if '|M|' in column:
                mag_atoms.extend(table.loc[:,'Label'])
            if 'Occupancy' in column:
                alloy = True
                for j, row in table.iterrows():
                    if row.Occupancy < 1:
                        dis_atoms.append(row['Atom type'])
                        dis_labels.append(row['Label'])
                        dis_occupancy[row.Label] = row.Occupancy
                        coord = [row.x,row.y,row.z]
                        dis_coord.append(coord)
            else:
                alloy = False
        elif 'Atom' in column:
            coord_tables.append(table)
        else:
            print('Error: Failed to parse MAGNDATA tables.')
            exit(1)

    # Find the positional overlapping of disorder atoms
    lapped_pairs = []
    if alloy:
        for i, label1 in enumerate(dis_labels):
            for j in range(i+1,len(dis_labels)):
                label2 = dis_labels[j]
                if dis_coord[i] == dis_coord[j]:
                    lapped_pairs.append([label1,label2])

    species = []
    ext_labels = []
    coords = []
    moments = []
    mag_constraints =[]
    occupancy = []
    for i, table in enumerate(coord_tables):
        column = table.columns.map(str)

        for j, row in table.iterrows():
            species.append(atom_types[i])
            ext_labels.append(labels[i])
            coord = [row.x,row.y,row.z]
            coords.append(convert_float(coord))
            if labels[i] in mag_atoms:
                moment = [row.Mx,row.My,row.Mz]
                moments.append(convert_float(moment))
                mag_constraints.append(row['Symmetry constraints on M'])
            else:
                moments.append([0.0] * 3)
                mag_constraints.append('0,0,0')
            if alloy:
                if labels[i] in dis_labels:
                    try:
                        occupancy.append(float(dis_occupancy[labels[i]]))
                    except:
                        print_exc()
                        exit(3)
                else:
                    occupancy.append(1.0)
            else:
                occupancy.append(1.0)

    elements = list(set(atom_types))
    max_elem_num = max_element_number(elements)

    num = float(sum(occupancy))
    if num.is_integer():
        defect = False
    else:
        defect = True

    sites = {
        'labels': ext_labels,
        'species': species,
        'coords': coords,
        'moments': moments,
        'occupancy': occupancy,
        'constraints': mag_constraints
    }

    info = {
        'num_of_atoms': num,
        'overlapped_pairs': lapped_pairs,
        'alloy': alloy,
        'defect': defect,
        'magnetic_atoms': mag_atoms,
        'disorder_atoms': dis_atoms,
        'elements': elements,
        'max_element_number': max_elem_num,
        }

    return sites, info

def process_html(htmlstr):
    # Merge the sepalate tables in MAGNDATA html (for pandas.read_html)

    # Decompose merged cells
    new = re.sub(r'<tr><td[^>]*colspan=[^>]*>.*?<\/td><\/tr>','',htmlstr)

    # Fix HTML
    soup = BeautifulSoup(new, "lxml")
    for body in soup("tbody"):
        body.unwrap()

    return str(soup)

def convert_float(vec):
    try:
        newvec = [float(x) for x in vec]
        return newvec
    except Exception as e:
        raise Exception("Problem in convert_float {}".format(e))

def strip_atom(atoms):
    spec = []
    for atom in atoms:
        a = re.search(r'[a-z|A-Z]+', atom).group()
        # Deuteriums are converted into simple hydrogen
        if a == 'D':
            a = 'H'
        spec.append(a)
    return spec

def max_element_number(elements):
    elements_nums = [Element(e).number for e in elements]
    return max(elements_nums)
