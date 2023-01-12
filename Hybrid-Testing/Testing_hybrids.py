from preprocessing import database, preprocessing_utils, clustering, merge_search, sqlite
from preprocessing.merge_search import modified_match_masses
from constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS
# import matplotlib.pyplot as plt

#These are user inputted parameters and their typical values
ppm_tolerance = 20
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 25
make_new_db = True

#Set your filepaths to the database and the spectra folder
prot_path = '/home/    /data/database/sample_database.fasta'
proteins = database.build(prot_path)

spectra_path = '/home /data/spectra/NOD2_E3'
spectra_files = preprocessing_utils.get_spectra_files(spectra_path)

#Loads in the spectra as a list of spectrum objects
spectra, _ = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)
    
for spectrum in spectra:
    #This matches every mz in a spectrum with a list of kmers it can match to. Format is (m/z, location_start, location_end, ion, charge, parent_protein)
    matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, make_new_db)
    
    #Your code here
