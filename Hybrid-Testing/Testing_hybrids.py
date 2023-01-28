from preprocessing import database, preprocessing_utils, clustering, merge_search, scoring
from preprocessing.constants import WATER_MASS, PROTON_MASS, AMINO_ACIDS
from preprocessing.sqlite import database_file
import matplotlib.pyplot as plt

#These are user inputted parameters and their typical values
ppm_tolerance = 20 
peak_filter = 25
relative_abundance_filter = .1
prec_tol = 10
max_pep_len = 10 #want at 25 eventually, but save time/resources to build at 10
make_new_db = False #change to false after first time running 

#Set your filepaths to the database and the spectra folder
prot_path = 'C:/Users/kayle/OneDrive/Documents/Hybrid-Testing/Hybrid-Testing/data/database/sample_database.fasta'
proteins = database.build(prot_path)

dbf = database_file(max_pep_len, make_new_db)

if make_new_db:
    kv_prots = [(k, v) for k, v in proteins.proteins]    
    merge_search.modified_make_database_set(kv_prots, max_pep_len, dbf)

spectra_path = 'C:/Users/kayle/OneDrive/Documents/Hybrid-Testing/Hybrid-Testing/data/spectra/NOD2_E3'
spectra_files = preprocessing_utils.get_spectra_files(spectra_path)

#Loads in the spectra as a list of spectrum objects
spectra = preprocessing_utils.load_spectra(spectra_files, ppm_tolerance, peak_filter, relative_abundance_filter)
    
#container stores result of if a spectrum is hybrid or not 
is_hybird = list

#This loop checks each spectrum to determine if it is a hybrid peptide
for spectrum in spectra:
    #This matches every mz in a spectrum with a list of kmers it can match to. Format is (m/z, location_start, location_end, ion, charge, parent_protein)
    matched_masses_b, matched_masses_y = merge_search.modified_match_masses(spectrum.mz_values, proteins, max_pep_len, ppm_tolerance, dbf)

    #Get precursor mass to check matches 
    b_precursor, y_precursor = list(matched_masses_b.keys())[-2], list(matched_masses_b.keys())[-1]

    #Getting everything that matched the precursor weights (there's a lot?)
    all_b_hits, all_y_hits = matched_masses_b[b_precursor], matched_masses_b[y_precursor]

    #Investigating the first hit for the b and y precursor weights 
    b_hit, y_hit = all_b_hits[0], all_y_hits[0]

    #getting the sequence associated with the hits 
    b_seq, y_seq = preprocessing_utils.find_sequence(b_hit, proteins.proteins), preprocessing_utils.find_sequence(y_hit, proteins.proteins)

    #scoring the hit 
    b_score, y_score = scoring.overlap_scoring(b_seq, ppm_tolerance, spectrum.mz_values), scoring.overlap_scoring(y_seq, ppm_tolerance, spectrum.mz_values)

    #use the score to guess if it's a hybrid? 
