import errno
import os

from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.vasp.sets import MPRelaxSet

class reduce_H_site:

    def vasp_input(path, file_name):
        '''
        This function define the quantum calcualtion values for input structure.
        Input:

        path : where the atomic structure files located
        file_name : file name of the atomic structure

        return:
        generate the input files for quantum calculation in the "path" location
        '''
        job_file = '/Users/shaolin_old/code/job_mo_vc_gam'                  #INPUT
        struct = Poscar.from_file(path + file_name, check_for_POTCAR = False, read_velocities = False).structure
        user_incar_settings = {'ISYM': 0, 'ISIF': 2, 'EDIFFG': -0.02,'NELM':300, 'EDIFF': 5e-7, 'NSW': 99, 'NPAR':4,'ISMEAR':0} # Input parameter to override MP parameter values
        input_set = MPRelaxSet(struct, user_incar_settings = user_incar_settings)
        input_set.write_input(path)
        paste_file=open(path + "job_mo",'w')
        file = open(job_file, "r")
        for line in file:
            paste_file.write(line)
        file.close()
        paste_file.close()
    #------end-----------#


    def print_sites(struct,write_path,ewald=True,write_file=True):
        '''
        This function write the input files for quantum calcualtion
        Input:
        struct: Structure of atomic coordinate
        write_path : where to write the input files
        ewald = Calculate the Ewald energy of the structure
        write_file = Write the files?

        return:
        Ewald energy (if ewald = True) of the structure and write files for quantum calcultion (write_file =ture)
        '''
        energy = []
        Hydrogen = []
        selective_dynamic = []
        hydrogen_less_structure = []
        structure_with_one_hydrogen = []
        for j in struct.sites:
            if bool(j.specie.symbol in "H"):
                Hydrogen.append(j)
            else:
                hydrogen_less_structure.append(j)


        hydrogen_count=0

        for i in Hydrogen:
            for j in hydrogen_less_structure:
                structure_with_one_hydrogen.append(j)
            structure_with_one_hydrogen.append(i)
            all_struct = Structure.from_sites(structure_with_one_hydrogen).get_sorted_structure()
            all_struct.add_oxidation_state_by_element({"Cr":4, "Mo":6, "Ta":5, "Cd":2, "Nd":3, "Nb":5, "Li":1, "Bi":5, "Dy":3, "Er":3, "Ho":3, "Lu":3, "Mg":2, "Pr":3, "Pb":2, "Zn":4, "Tm":3, "Re":7, "Ru":4, "Ir":4, "Os":4, "Sb":5, "W":6, "Te":4, "U":4, "Ga":3, "Mn":4, "La":3, "Tb":3, "Th":4, "Y":3, "Ac":3, "Sc":3, "Eu":2, "V":4, "Ca":2, "Ba": 2, "Zr":4, "Yb": 3, "Ti": 4, "O": -2, "H": 1, "In": 3, "Cu":2, "Na":1, "Co":3, "Sr":2, "Ni":3, "Sm":3, "Ce":4, "Fe":3})  # dictionary holding oxidation state of different elements for Ewald energy calcualtion
            for site in all_struct.sites:
                if bool(site.specie.symbol in "H"):
                    selective_dynamic.append([True, True, True])
                else:
                    selective_dynamic.append([False, False, False])

            # Calculate Ewald energy of atomic strucure
            if ewald:
                energy.append([hydrogen_count, EwaldSummation(all_struct).total_energy])

            # Write quantum calculation input files for atomic system

            if write_file:
                Poscar_to_write = Poscar(all_struct, selective_dynamics = selective_dynamic)
                if not os.path.exists(os.path.dirname(write_path + str(hydrogen_count) + "/")):
                    try:
                        os.makedirs(os.path.dirname(write_path + str(hydrogen_count) + "/"))
                    except OSError as exc: # Guard against race condition
                        if exc.errno != errno.EEXIST:
                            raise
                Poscar_to_write.write_file(write_path + str(hydrogen_count) + "/main.vasp")  # Write the atomic structure file in the write_path directory
                vasp_input(write_path + str(hydrogen_count)+"/", "main.vasp")        # Write the input files of the atomic structure for quantum computation.
            structure_with_one_hydrogen = []
            hydrogen_count = hydrogen_count + 1
        return energy
    #------end-----------#
