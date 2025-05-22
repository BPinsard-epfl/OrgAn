from OrgAn import get_elution_order, calculate_polarity_index, estimate_retention_factor, generate_chromatogram, gives_data_frame


solutes = gives_data_frame("tests/test_smiles.csv")
print(get_elution_order(solutes, True)["name"].tolist())


def test_polarity_index_value_error(value1 = 0, value2 = 0, value3 = 0, value4 = 0):
    try:
        calculate_polarity_index(value1, value2, value3, value4)
    except ValueError:
        return True
    except:
        return False
    else:
        return False
    
def test_polarity_index_type_error(value1):
    try:
        calculate_polarity_index(value1)
    except TypeError:
        return True
    except:
        return False
    else:
        return False
    
def test_generate_chromatogram():
    try:
        generate_chromatogram(solutes, 5, 3)
    except:
        return False
    else:
        return True


assert get_elution_order(solutes)["name"].tolist() == ['butane', 'furan', 'acetylene', 'Propionic Acid', 'Ethanol', 'Acetone', 'Acetic Acid', 'DL-alanine'], "get_elution_order did not give the correct elution order"
assert get_elution_order(solutes, True)["name"].tolist() == ['DL-alanine', 'Acetic Acid', 'Ethanol', 'Acetone', 'Propionic Acid', 'acetylene', 'furan', 'butane'], "get_elution_order does not give the correct reverse phase elution order"
print ("get_elution_order passed all tests")

assert calculate_polarity_index(water=1) == 10.2, "calculate_polarity_index did not calculate the right polarity index"
assert round(calculate_polarity_index(etoac=0.4, et2o=0.6), 2) == 3.44, "calculate_polarity_index did not calculate the right polarity index"
assert round(calculate_polarity_index(n_hex=0.2, ccl4=0.3, etoh=0.5), 2) == 2.65, "calculate_polarity_index did not calculate the right polarity index"
assert test_polarity_index_type_error("not a number"), "calculate_polarity_index did not throw the right exception"
assert test_polarity_index_value_error(value1=0.2, value2= 0.3), "calculate_polarity_index did not throw the right exception"
assert test_polarity_index_value_error(value1=0.4, value2=0.2, value3=0.7, value4=0.1), "calculate_polarity_index did not throw the right exception"
print("calculate_polarity_index passed all tests")

assert round(estimate_retention_factor(1.34, 10.2), 2) == 21.88, "estimate_retention_factor did not yield the correct result"
assert round(estimate_retention_factor(1.34, 8.2), 2) == 2.19, "estimate_retention_factor did not yield the correct result"
assert round(estimate_retention_factor(1.14, 10.2), 2) == 13.80, "estimate_retention_factor did not yield the correct result"
print("estimate_retention_factor passed all the tests")

assert test_generate_chromatogram, "generate_chromatogram raised an error"
print("generate_chromatogram passed all the tests")