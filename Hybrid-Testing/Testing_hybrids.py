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

    #FIRST ROUTE TO CHECK IF HYBRID - CHECK PRECURSOR WEIGHT MATCHES 
    #Get precursor mass
    b_precursor, y_precursor = list(matched_masses_b.keys())[-2], list(matched_masses_y.keys())[-1]

    #Getting everything that matched the precursor weights
    all_b_hits, all_y_hits = matched_masses_b[b_precursor], matched_masses_y[y_precursor]

    #Loop over all precursor weight hits and score each one 
    b_scores, y_scores = list(), list()
    for hit in all_b_hits: 
        #Obtain sequence from hit 
        seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
        #Score the hit - score = number of peaks matched (out of 25)
        b_scores.append(scoring.overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    for hit in all_y_hits: 
        seq = preprocessing_utils.find_sequence(hit, proteins.proteins)
        y_scores.append(scoring.overlap_scoring(seq, ppm_tolerance, spectrum.mz_values))
    
    #Use the scores to guess if it's a hybrid?
    #Note that sometimes this divides by zero, meaning there are 0 hits for the precursor match 
    if(len(b_scores != 0)):
        b_mean = sum(b_scores)/len(b_scores)
    else: 
        b_mean = 0
    if(len(y_scores != 0)):
        y_mean = sum(y_scores)/len(y_scores)
    else:
        b_mean = 0

    #Perhaps if the mean scores are below some threshold, they are hybrid? 
        #Need to truth this with some instances 

    #SECOND ROUTE TO CHECK IF HYBRID - CHECK HIGHEST MASS MATCH 
    #Best scoring, highest mass match (with charge accounted for...)
    #How does it compare to the precursor weight? 
    #If there is a mass close to the precursor with a 'high' score, it may not be hybrid 

