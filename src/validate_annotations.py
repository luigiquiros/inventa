import pandas as pd
import numpy as np
def MF_checker_on_two_MF_columns(table, column1, column2, check_column_name):
        # We are matching molecular formulas only by the number of carbons as the SIRIUS adducts can be quite variable and deals strangely with adduct/ionisationin the MF. Like loss of water etc. Ammonium adduct.
    table[str(column1)+'_partial'] = table[column1].str[:3]
    table[str(column2)+'_partial'] = table[column2].str[:3]
​
    # Make functions to check the partial match
    def mf_check(x):    
        return 'yes' if x[str(column1)+'_partial'] == x[str(column2)+'_partial'] else 'no'
​
    # Apply the functions for all valid MF pairs above zodiac score
    table[check_column_name] = table.apply(mf_check, axis=1)