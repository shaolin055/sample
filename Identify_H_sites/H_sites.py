import os
import csv

from numpy import pi, cos, sin, arccos, arange

from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer



from lib.lib1 import reduce_H_site

class hydrogen_site:


    def sym_site(struct, element):
        '''
        Identifies the symmetric oxygen in a oxide materials atomic structure. and return indices that are unique
        Input:
        struct: Atomic structure for symmetry analysis
        element : which element to find symmetry in the structure. e.g. "O", "C" etc.
        '''

        analyzer_init = SpacegroupAnalyzer(struct)
        symmetry_analysis = SymmetrizedStructure(struct, analyzer_init.get_space_group_operations(), analyzer_init.get_symmetry_dataset()["equivalent_atoms"], analyzer_init.get_symmetry_dataset()['wyckoffs'])
        indice=[]
        for equivalent_indices in symmetry_analysis.equivalent_indices:
            if struct[equivalent_indices[0]].specie.symbol in element:
                indice.append(equivalent_indices[0])
        return indice



    def create_h_sites(path, file_name,write_path):
        '''
        Generate initial H site locations around symmetric oxygens and generate input atomic coordinate of differnet materials files for quantum calculation

        Input :
        path = Input path of materials (without H coordinate) file location
        file_name = file name of the atomic coordinate of the materials.
        write_path = The path where the atomic coordinate of the materials with H location will be written

        return:
        Generate multiple diectory for each atomic structure file. Each directory contains input file for quantum calculation using VASP module.
        The generated atomic coordinate file contains atomic coordinate of the structure along with one possible location of hydrogen to relax to local minimum position.
        '''

        # Read structure
        struct = Poscar.from_file(path + file_name).structure

        # Identify 32 location with 1Å radius around a position


        distance = 1
        energy_tolarance = 200
        scale_matrix = [1, 1, 1]
        num_of_points = 32
        indices = arange(0, num_of_points, dtype = float) + 0.5
        phi = arccos(1 - 2 * indices / num_of_points)
        theta = pi * (1 + 5 ** 0.5) * indices
        x_coord, y_coord, z_coord = distance* cos(theta) * sin(phi), distance * sin(theta) * sin(phi), distance * cos(phi)

        X=[]
        for iter in range(0,len(x)):
            X.append([x_coord[iter],y_coord[iter],z_coord[iter]])

        #Copy structure for modification

        st=[]
        for iter in struct.sites:
            st.append(iter)
        all_struct = Structure.from_sites(st).get_sorted_structure()

        Hydrogen=[]
        # Identify symmetric oxygen and mark possible proton site around oxygen.

        symmetry_indices=sym_site(struct,'O')
        count = 0
        for sites_iter in struct.sites:
            if count in symmetry_indices:
                hydrogen_coordinate = sites_iter.coords
                for k in range(0,len(x_coord)):
                    Hydrogen.append(X[k] + hydrogen_coordinate)
                    struct.append("H", X[k] + hydrogen_coordinate, coords_are_cartesian = True )
            count = count + 1
        print(struct.composition)

        # Calculate Ewald energy and filter high energy structure

        structure_energy = reduce_H_site.print_sites( struct,write_path, write_file = False)
        sorted_energy = sorted(structure_energy, key = lambda x: x[1])
        filtered_energy = [x for x in sorted_energy if x[1] < min(sorted_energy[1]) + energy_tolarance]
        all_struct.make_supercell(scale_matrix)

        # Incorporate all H in the structure.

        for iter1 in filtered_energy:
            for iter2 in range(0, len(Hydrogen)):
                if iter1[0] == iter2:
                    all_struct.append("H",Hydrogen[iter2], coords_are_cartesian = True )

        # Write structure with possible proton sites for further quantum calculation

        reduce_H_site.print_sites( all_struct, write_path, ewald = False, write_file = True)
        Poscar_to_write = Poscar(all_struct)
        Poscar_to_write.write_file(write_path + "H_site.vasp")

    #------end-----------#


#Read atomic structures from text files

path='/'
file='structure_list.txt'
file_reader = csv.reader(open(file, 'rU'), dialect = csv.excel_tab)
folder_name =[]
count=0

#--------run for every structure------#

for iter1 in file_reader:
    split = (iter1[0].split(','))
    structure_number = split[0]
    hydrogen_site.create_h_sites(path, structure_number + '.vasp', path + structure_number)
    count = count + 1
