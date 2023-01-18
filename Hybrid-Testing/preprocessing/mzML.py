from .utils import file_exists
from .objects import Spectrum
from pyteomics import mzml

def relative_abundance_filtering(
    masses: list, 
    abundances: list, 
    percentage: float
    ):
    '''Take all peaks from the spectrum who's abundance is at least *percentage* 
    of the total abundances. It is assumed that the masses and abundances lists 
    share ordering

    :param masses: m/z values 
    :type masses: list
    :param abundances: abundance value for the m/z values. Abundance at entry 
        *i* corresponds to m/z valuat entry *i*
    :type abundances: list
    :param percentage: the minimum percentage of the total abundance a peak must
        have to pass the filter. Values are in the range [0, 1). A relatively 
        realistic value is .005 (.5%)
    :type percentage: float

    :returns: filtered masses, filtered abundaces
    :rtype: (list, list)
    '''

    # total intensity
    ti = sum(abundances)

    # find the filter value
    min_value = ti * percentage

    # zip them up and take values that pass the filter
    filtered_mass_abundances = [x for x in zip(masses, abundances) if x[1] >= min_value]

    # split them off and return
    masses = [float(x) for x, _ in filtered_mass_abundances]
    abundances = [float(x) for _, x in filtered_mass_abundances]

    return (masses, abundances)

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

def read(filename: str, peak_filter=0, relative_abundance_filter=0) -> list:
    '''
    read an .mzML file into memory. Filter out the peaks by the type specified.
    If both filter types are set to 0, all peaks are returned, otherwise filtering 
    will happen. If both filters are given a value, peak_filter is given preference.

    Inputs:
        filename:       (str) path to the file to import
    kwargs:
        peak_filter:                (int) the top number of peaks to keep. Default=0
        relative_abundance_filter:  (float) thee percentage of abundance a peak must have to be
                                    considered. Must be in range(0, 1) or the integer is converted
                                    to a decimal in that range. Default=0
    Outputs:
        (list) Spectrum namedtuple instances
    '''
    if not file_exists(filename):
        print('File {} not found. Please make sure that this file exists'.format(filename))
        return

    spectra = []
    
    filecontents = mzml.read(filename)

    content: dict

    for spec_num, content in enumerate(filecontents):

        masses = list(content['m/z array'])
        abundances = list(content['intensity array'])

        # peak filter if the number is > 0
        if peak_filter > 0:
            masses, abundances = peak_filtering(masses, abundances, peak_filter)

        # if peak filter is not set and relative abundances is, filter by that
        elif relative_abundance_filter > 0:

            # if an integer is given, make it a float in the range (0, 1)
            while relative_abundance_filter > 1:
                relative_abundance_filter /= 100

            masses, abundances = relative_abundance_filtering(masses, abundances, relative_abundance_filter)

        # get the total intensity
        ti = sum(abundances)

        # get the precursor and its charge
        # we will assume its the first entry in the list
        precursor = None
        precursor_charge = 0

        if not len(content['precursorList']['precursor']) or not len(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon']):
            precursor = max(masses)
            precursor_charge = 1

        else:
            precursor = float(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
            precursor_charge = int(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])

        # get the id
        id = content.get('id', '')

        spectra.append(
            Spectrum(
            spec_num,
            masses,
            abundances,
            precursor,
            precursor_charge,
            filename, 
            id
        ))

    return spectra