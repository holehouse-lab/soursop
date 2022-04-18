import pytest
import numpy as np
from soursop import ssmutualinformation


def test_shan_entropy():

    assert 0.693 == np.round(ssmutualinformation.shan_entropy(np.repeat(1,2)),3)
    assert 1.099 == np.round(ssmutualinformation.shan_entropy(np.repeat(1,3)),3)
    assert 1.386 == np.round(ssmutualinformation.shan_entropy(np.repeat(1,4)),3)



def test_calc_MI():

    X = np.random.normal(100,0.1,size=10000)
    Y = np.random.normal(100,0.1,size=10000)
    bins = np.arange(0,300,2)

    assert 0.693 == np.round(ssmutualinformation.calc_MI(X,X,bins),3)
    

