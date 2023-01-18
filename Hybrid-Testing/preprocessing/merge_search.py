from .preprocessing_utils import ppm_to_da
#from .sqlite import database_file
from .objects import Database
from .gen_spectra import max_mass
import shutil
import sys
import time

def get_data(kmer, start, end, protein_num):
    data_list = []
    for ion in 'by':
        for charge in [1,2]:
            mass = max_mass(kmer, ion=ion, charge=charge)
            ion_int = 0 if ion == 'b' else 1
            input_tuple = (mass, start, end, ion_int, charge, protein_num)
            data_list.append(input_tuple)

    return data_list

def db_make_set_for_protein(i,prot,max_len, dbf, data):
    seq_len = len(prot)
    count_max = 1000000
    for size in range(2, max_len + 1):
        # size -> [2, max_len]
        for start in range(0, seq_len - size + 1):
            end = start + size
            kmer = prot[start:end]
            # last_index = seq - size 6, end = start + size - 1 = 7
            # [data.append(x) for x in get_data(kmer, start, end)]
            data_list = get_data(kmer, start, end, i)
            data.extend(data_list)
            # insertion code
            if len(data) > count_max:
                dbf.insert(data)
                data.clear()
            
    return

def db_make_database_set_for_proteins(proteins,max_len,dbf):
    plen = len(proteins)
    last_percent = 0
    data = []
    
    for i, (_, prot_entry) in enumerate(proteins):
        percent = int((i+1) * 100 / plen)
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        if percent != last_percent:
            # print(f'\rInserting {percent}%', end='')
            last_percent = percent
            free = shutil.disk_usage('/')[2]
            free = free/(1024**3)
            if free < 10:
                print("\nUsed too much space, Space available =", free, "GB" )
                sys.exit(1)
        db_make_set_for_protein(i,prot_entry,max_len, dbf, data)
        
    if len(data) != 0:
        dbf.insert(data)
        
def modified_make_database_set(proteins: list, max_len: int, dbf):    
    print("\nBeginning Insertions")
    start = time.time()
    db_make_database_set_for_proteins(proteins,max_len, dbf)
    duration = time.time() - start
    print("Insertion took: ", duration)
    # db.read() #Only for debugging
    print('\nIndexing the set of kmers based on mass, ion')
    dbf.index_ion_mass()
    print('\nIndexing the set of kmers based on protein, start position, end position')
    dbf.index_ion_mass_b()
    print('\nIndexing the set of kmers based on protein, end position, start position')
    dbf.index_ion_mass_y()
    print('Done making database')
    
    
    # print('\nSorting the set of protein masses...')
    # kmer_list = []
    # # handle_sorting_keys(db_dict_b, db_dict_y, kmer_list)
    # kmer_list = sorted(kmer_list, key=lambda x: x[0])
    # print('Sorting the set of protein masses done')
    return

def modified_match_masses(input_masses: list, db: Database, max_len: int, ppm_tolerance, dbf):
    # max_boundary = max(boundaries.keys())
    # estimated_max_len = ceil(boundaries[max_boundary][1] / 57.021464)
    # max_len = min(estimated_max_len, max_pep_len)
    
    #dbf = database_file(max_len, False)
    
    matched_masses_b, matched_masses_y = dict(), dict()
    
    for input_mass in input_masses:
        tol = ppm_to_da(input_mass, ppm_tolerance)
        # ion_int = 0 if ion == 'b' else 1
        matched_masses_b[input_mass], matched_masses_y[input_mass] = dbf.query_mass(input_mass, tol) #same place: location start, protein_num
        # print(input_mass, 'b', matched_masses_b[input_mass])
        # print(input_mass, 'y', matched_masses_y[input_mass])
        
        # if debug:
        # write_matched_masses(write_path, matched_masses_b, matched_masses_y, kmer_set, debug)

    return matched_masses_b, matched_masses_y