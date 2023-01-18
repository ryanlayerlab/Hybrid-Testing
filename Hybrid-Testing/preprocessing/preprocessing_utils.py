import os
from file_io import spectra
from gen_spectra import convert_precursor_to_ion

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
        these_spectra = spectra.load(
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

def peak_filtering(masses: list, abundances: list, num_peaks: int):
    '''Take the most abundant peaks and return the sorted masses with the abundances.
    It is assumed that the masses and abundances lists share ordering


    :param masses: m/z values 
    :type masses: list
    :param abundances: abundance value for the m/z values. Abundance at entry 
        *i* corresponds to m/z valuat entry *i*
    :type abundances: list
    :param num_peaks: the top X most abundant peaks 
    :type num_peask: int

    :returns: filtered masses, filtered abundaces
    :rtype: (list, list)
    '''

    # zip the abundance and the m/z values together
    mass_abundances = zip(masses, abundances)
    
    # sort by key 1, the abundance, and take the top peak filter results
    mass_abundances = sorted(mass_abundances, key=lambda x: x[1], reverse=True)[:num_peaks]

    # sort them now by the value of m/z
    mass_abundances.sort(key=lambda x: x[0])

    # seperate them
    masses = [float(x) for x, _ in mass_abundances]
    abundances = [float(x) for _, x in mass_abundances]

    return (masses, abundances)