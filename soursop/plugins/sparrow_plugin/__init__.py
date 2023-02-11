#
# In general, we recommend placing user-facing pluggin code within the
# __init__.py file either as a stateless function or as a class which
# enables statefull functionality. Examples of both (achieving the same
# outcome) are shown below
#
#
#


# NOTE - we do not want to make plugins hard dependencies for SOURSOP,
# so if you are providing a plugin that requires additional packages
# which may or may not be avaialble please ensure that any import statements
# are wrapped as a try/except block as is done below with an appropriate
# message printed if the import fails
#
# We also recommend importing using the as _<name> format - this means that
# the imported module cannot be imported from the calling Python code, which
# in general helps maintain modularity
#
try:
    from sparrow import Protein as _Protein
except ModuleNotFoundError:
    print('The sparrow_plugin reqires sparrow. This can be installed from [https://github.com/idptools/sparrow]')
    raise Exception



# Stateless stand-alone function that takes in an SSProtein object and
# returns the NCPR 
#
#
def get_protein_net_charge(ssprotein_obj):
    """
    Returns the net charge of the protein

    Parameters
    ----------------
    soursop.ssprotein.SSProtein ssprotein_obj
        A passed SSProtein object

    Returns
    ----------
    float
        Returns the net charge per residue associated with the
        passed sequence (only works for the 20 standard amino 
        acids).
    """

    # get sequence
    tmp = ssprotein_obj.get_amino_acid_sequence(oneletter=True)

    # remove caps
    tmp = tmp.replace('>','')
    tmp = tmp.replace('<','')

    # calculate NCPR
    return _Protein(tmp).NCPR



# Statefull class that provides an interface into a sparrow.protein.Protein
# object
#
#
class SparrowProtein:

    def __init__(self, ssprotein_obj):
        """
        Returns an interface into the 

        Parameters
        ----------------
        soursop.ssprotein.SSProtein ssprotein_obj
            A passed SSProtein object

        Returns
        ----------
        SparrowProtein object
        

        """
        
        # get sequence
        tmp = ssprotein_obj.get_amino_acid_sequence(oneletter=True)

        # remove caps
        tmp = tmp.replace('>','')
        tmp = tmp.replace('<','')


        
        self.protein = _Protein(tmp)

    @property
    def NCPR(self):
        return self.protein.NCPR
                               
