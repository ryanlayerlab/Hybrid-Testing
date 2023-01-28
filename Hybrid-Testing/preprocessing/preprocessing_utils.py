import os
from .spectra import load
from .gen_spectra import convert_precursor_to_ion

def get_spectra_files(spectra_folder):
    spectra_files = []
    for (root, _, filenames) in os.walk(spectra_folder):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))
    return spectra_files

def ppm_to_da(mass: float, ppm_tolerance: float) -> float:
    return abs((ppm_tolerance / 1000000)*mass)

def load_spectra(
    spectra_files: list, ppm_tol: int, peak_filter: int = 0, relative_abundance_filter: float = 0.0):
    linear_spectra = []
    all_spectra = []
    for spectra_file in spectra_files:
        these_spectra = load(
            spectra_file, 
            peak_filter=peak_filter, 
            relative_abundance_filter=relative_abundance_filter
        )
        for object in these_spectra:
            b_prec, y_prec = convert_precursor_to_ion(object.precursor_mass, object.precursor_charge)
            object.mz_values.append(b_prec)
            object.mz_values.append(y_prec)
            object.abundance.append(100000000) #I gave it a crazy high abundance to represent precursor. Still a hack
            object.abundance.append(100000000)
        all_spectra += these_spectra
        # leave commented; uncomment only to test just specific indices
        # index_list = [0] #For using a condensed database
        # these_spectra, all_spectra = reduce_database(all_spectra, these_spectra, index_list) 
        linear_spectra += list(set([
            x for spectrum in these_spectra for x in spectrum.mz_values
        ]))
    linear_spectra.sort()
    return (all_spectra)

def find_sequence(hit, protein_list):
    pid = hit[5]
    start, end = hit[1], hit[2]
    prot_seq = protein_list[pid][1]
    return prot_seq[start:end]