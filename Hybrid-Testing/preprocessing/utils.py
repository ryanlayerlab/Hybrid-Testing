import os

def file_exists(file_name: str) -> bool:
    '''Determine if a file exists

    :param file_name: Path to the file in question
    :type file_name: str
    
    :returns: True if the file exists
    :rtype: bool
    '''
    
    return os.path.isfile(file_name)
