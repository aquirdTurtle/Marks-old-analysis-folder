
def import_HP4395A_Data(filename):
    import pandas as pd
    return pd.read_csv(filename, delimiter='\t', header=11)


def import_SRSSR780_Data(filename):
    import pandas as pd
    return pd.read_csv(filename, delimiter=',', header=None)

