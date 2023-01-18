from .objects import Database
from pyteomics import fasta

def build(fasta_file: str) -> Database:
    '''Create a Database namedtuple from a fasta file

    :param fasta_file: the full path to a fasta database file 
    :type fasta_file: str

    :returns: a Database object with the fasta file and protein fields filled in
    :rtype: Database
    '''

    db = Database(fasta_file)
    prots = []
    # prots = defaultdict(list)
    # get_name = lambda x: x.split('|')[-1].split()[0]
    for entry in fasta.read(fasta_file):
        # p_name = get_name(entry.description)
        # prots[p_name].append(entry)
        prots.append(entry)
    db = db._replace(proteins=prots)
    return db