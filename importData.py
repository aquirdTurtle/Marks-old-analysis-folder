def importData(date, name, runNumber):
    """
    This function imports the data created by my first-step python code.
    :param date: the date in YYMMDD format so that the code knows where to look for the data.
    :param name: the name of the compiled data file, minus _run?? (e.g. it's standardized for carrier calibrations
    :param runNumber: The run number.
    :return: returns an OrderedDictionary object from the collections library that contains all the information about a
             run.
    """
    from collections import OrderedDict as dic
    data = dic()
    with open('//andor/share/Data and documents/Data repository/' + str(date) + '/' + str(name) + '_run'
              + str(runNumber) + '.csv') as myfile:
        fileText = myfile.read().replace('\n', '')
        rows = []
        for row in fileText.split(':'):
            rows.append(row.replace('\t', ''))
        del rows[0]
        for itemInc in range(0, len(rows), 2):
            if (rows[itemInc + 1].find(';')) != -1:
                itemData = dic()
                subrows = []
                for subrow in rows[itemInc + 1].split(';'):
                    subrows.append(subrow)
                del subrows[0]
                for subitem in range(0, len(subrows), 2):
                    itemData[subrows[subitem]] = subrows[subitem + 1]
                data[rows[itemInc]] = itemData
            else:
                data[rows[itemInc]] = rows[itemInc + 1]
    return data
