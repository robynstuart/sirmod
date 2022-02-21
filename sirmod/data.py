'''
Function to create and manipulate parameters
'''

# Imports
import numpy as np
import pandas as pd
import sciris as sc
import sciris as sc
from . import utils as smu
from xlrd import open_workbook, colname


#%% Helper functions
def force_bool(entry, location=''):
    ''' Convert an entry to be Boolean '''
    if entry in [1, 'TRUE', 'true', 'True', 't', 'T']:
        return 1
    elif entry in [0, 'FALSE', 'false', 'False', 'f', 'F']:
        return 0
    else:
        errormsg = f'Boolean data {entry} not understood in spreadsheet location {location}'
        raise Exception(errormsg)


def validate_data(thesedata, sheetname, thispar, row, checkupper=False, checklower=True, checkblank=True, startcol=0):
    ''' Do basic validation on the data: at least one point entered, between 0 and 1 or just above 0 if checkupper=False '''

    # Check that only numeric data have been entered
    for column, datum in enumerate(thesedata):
        if not sc.isnumber(datum):
            errormsg = f'Invalid entry in sheet {sheetname}, parameter {thispar}:\n'
            errormsg += f'row={row + 1}, column={colname(column + startcol)}, value={datum}\n'
            errormsg += 'Be sure all entries are numeric'
            if ' ' or '\t' in datum: errormsg += ' (there seems to be a space or tab)'
            raise Exception(errormsg)

    # Now check integrity of data itself
    validdata = np.array(thesedata)[~np.isnan(thesedata)]
    if len(validdata):
        valid = np.array([True] * len(validdata))  # By default, set everything to valid
        if checklower: valid *= np.array(validdata) >= 0
        if checkupper: valid *= np.array(validdata) <= 1
        if not valid.all():
            invalid = validdata[valid == False]
            errormsg = 'Invalid entry in sheet "%s", parameter "%s":\n' % (sheetname, thispar)
            errormsg += 'row=%i, invalid="%s", values="%s"\n' % (row + 1, invalid, validdata)
            errormsg += 'Be sure that all values are >=0 (and <=1 if a probability)'
            errormsg += versioncheck
            raise Exception(errormsg)
    elif checkblank:  # No data entered
        errormsg = 'No data or assumption entered for sheet "%s", parameter "%s", row=%i' % (sheetname, thispar, row)
        errormsg += versioncheck
        raise Exception(errormsg)
    else:
        return None


def blank2nan(thesedata):
    ''' Convert a blank entry to a nan '''
    return list(map(lambda val: np.nan if val == '' else val, thesedata))


def get_years(sheetdata):
    ''' Get years from a worksheet'''
    years = []  # Initialize epidemiology data years
    for col in range(sheetdata.ncols):
        thiscell = sheetdata.cell_value(1, col)  # 1 is the 2nd row which is where the year data should be
        if thiscell == '' and len(years) > 0:  # We've gotten to the end
            lastdatacol = col  # Store this column number
            break  # Quit
        elif thiscell != '':  # Nope, more years, keep going
            years.append(float(thiscell))  # Add this year

    return lastdatacol, years


#%% Data specs
def load_data_specs(filename=None, folder=None, verbose=2):
    ''' Function to parse the data parameter definitions '''
    inputsheets = ['Data inputs']
    fullpath = sc.makefilepath(filename=filename, folder=folder)
    workbook = open_workbook(fullpath)

    data_specs = sc.odict()
    for inputsheet in inputsheets:
        sheet = workbook.sheet_by_name(inputsheet)
        rawpars = []
        for rownum in range(sheet.nrows - 1):
            rawpars.append({})
            for colnum in range(sheet.ncols):
                attr = str(sheet.cell_value(0, colnum))
                cellval = sheet.cell_value(rownum + 1, colnum)
                if cellval == 'None': cellval = None
                rawpars[rownum][attr] = cellval
        data_specs[inputsheet] = rawpars

    sheets = sc.odict()  # Lists of parameters in each sheet
    sheettypes = sc.odict()  # The type of each sheet -- e.g. time parameters or matrices
    checkupper = sc.odict()  # Whether or not the upper limit of the parameter should be checked
    sheetcontent = sc.odict()
    for par in data_specs['Data inputs']:
        if par['sheet'] not in sheets.keys():  # Create new list if sheet not encountered yet
            sheets[par['sheet']] = []  # Simple structure for storing a list of parameter names, used in loadspreadsheet
            sheetcontent[par['sheet']] = []  # Complex structure for storing all information, used in makespreadsheet
        sheets[par['sheet']].append(par['short'])  # All-important: append the parameter name
        sheetcontent[par['sheet']].append(par)  # Append entire dictionary
        sheettypes[par['sheet']] = par['type']  # Figure out why kind of sheet this is
        checkupper[par['short']] = True if par['rowformat'] in ['decimal',
                                                                'percentage'] else False  # Whether or not to check the upper limit

    # Store useful derivative information
    data_specs['sheets'] = sheets
    data_specs['sheetcontent'] = sheetcontent
    data_specs['sheettypes'] = sheettypes
    data_specs['checkupper'] = checkupper

    return data_specs


#%% Load data
def load_data(data_sheet='simple.xlsx', verbose=1):
    '''
    Loads the data spreadsheet
    '''

    fullpath = sc.makefilepath(filename=data_sheet, ext='xlsx', default='simple.xlsx')
    sc.printv('Loading data from %s...' % fullpath, 1, verbose)

    # Create dictionary of parameters to load
    data_specs = load_data_specs(filename='model-inputs.xlsx', folder='../sirmod', verbose=verbose)
    sheets = data_specs['sheets']
    sheettypes = data_specs['sheettypes']
    checkupper = data_specs['checkupper']

    ## Initialize dictionaries
    data = sc.odict()  # Create sheetsure for holding data
    data['meta'] = sc.odict()
    data['meta']['date'] = sc.now()
    data['meta']['sheets'] = sc.dcp(sheets)  # Store parameter names

    ## Initialize populations
    data['pops'] = sc.odict()  # Initialize to empty list
    data['pops']['short'] = []  # Store short population/program names
    data['pops']['long'] = []  # Store long population/program names
    data['pops']['male'] = []  # Store whether or not population is male
    data['pops']['female'] = []  # Store whether or not population is female
    data['pops']['age'] = []  # Store the age range for this population

    ## Initialize partnerships
    data['pships'] = []  # Store partnerships

    ## Initialize other quantities
    blhindices = {'best': 0, 'low': 1, 'high': 2}  # Define best-low-high indices
    skipblanksheets = ['Optional indicators']  # Don't check optional indicators, check everything else
    skipblankpars = []

    ## Actually open workbook
    try:
        workbook = open_workbook(fullpath)  # Open workbook
    except Exception as E:
        errormsg = 'Failed to load spreadsheet "%s": %s' % (fullpath, repr(E))
        raise Exception(errormsg)

    ## Open workbook and calculate columns for which data are entered, and store the year ranges
    sheetdata = workbook.sheet_by_name('Population size')  # Load this workbook
    lastdatacol, data['years'] = get_years(sheetdata)
    assumptioncol = lastdatacol + 1  # Figure out which column the assumptions are in; the "OR" space is in between

    ##################################################################
    ## Now, actually load the data
    ##################################################################

    ## Load population data
    popssheet = workbook.sheet_by_name('Populations')
    for row in range(popssheet.nrows):
        thesedata = popssheet.row_values(row, start_colx=2,
                                         end_colx=11)  # Data starts in 3rd column, finishes in 11th column
        subparam = popssheet.cell_value(row, 1)  # Figure out which population it is
        if subparam != '':  # Warning -- ugly to duplicate this, but doesn't follow the format of data sheets, either
            sc.printv('Loading population "%s"...' % subparam, 4, verbose)
            data['pops']['short'].append(str(thesedata[0]))
            data['pops']['long'].append(str(thesedata[1]))
            data['pops']['male'].append(force_bool(thesedata[2], 'male, row %i' % row))
            data['pops']['female'].append(force_bool(thesedata[3], 'female, row %i' % row))
            data['pops']['age'].append([int(thesedata[4]), int(thesedata[5])])

    ##  Loop over each group of data sheets
    for sheetname in sheets.keys():  # Loop over each type of data
        subparlist = sheets[sheetname]  # List of subparameters
        sheetdata = workbook.sheet_by_name(sheetname)  # Load this workbook
        sheettype = sheettypes[sheetname]  # Get the type of this sheet -- e.g., is it a time parameter or a matrix?
        sc.printv('Loading sheet "%s"...' % sheetname, 3, verbose)

        # Loop over each row in the workbook, starting from the top
        for row in range(sheetdata.nrows):
            paramcategory = sheetdata.cell_value(row, 0)  # See what's in the first column for this row
            subparam = sheetdata.cell_value(row,
                                            1)  # Get the name of a subparameter

            if paramcategory != '':  # It's not blank
                sc.printv('Loading parameter "%s"...' % paramcategory, 3, verbose)

                # It's anything other than the constants sheet: create an empty list
                if sheettype != 'constant':
                    try:
                        thispar = subparlist.pop(0)  # Get the name of this parameter, e.g. 'popsize'
                    except:
                        errormsg = 'Incorrect number of headings found for sheet "%s"\n' % sheetname
                        errormsg += 'Check that there is no extra text in the first two columns'

                        raise Exception(errormsg)
                    data[thispar] = []  # Initialize to empty list

            elif subparam != '':  # The second column isn't blank: it's time for the data
                sc.printv('Parameter: %s' % subparam, 4, verbose)

                # It's key data, save both the values and uncertainties
                if sheettype == 'key':
                    startcol = 3  # Extra column for high/best/low
                    if len(data[thispar]) == 0:
                        data[thispar] = [[] for z in range(3)]  # Create new variable for best, low, high
                    thesedata = blank2nan(sheetdata.row_values(row, start_colx=startcol,
                                                               end_colx=lastdatacol))  # Data starts in 4th column -- need room for high/best/low
                    assumptiondata = sheetdata.cell_value(row, assumptioncol)
                    if assumptiondata != '': thesedata = [
                        assumptiondata]  # Replace the (presumably blank) data if a non-blank assumption has been entered
                    blh = sheetdata.cell_value(row, 2)  # Read in whether indicator is best, low, or high
                    data[thispar][blhindices[blh]].append(thesedata)  # Actually append the data
                    validate_data(thesedata, sheetname, thispar, row, checkblank=(blh == 'best'),
                                 checkupper=checkupper[thispar],
                                 startcol=startcol)  # Make sure at least the best estimate isn't blank

                # It's basic data, append the data and check for programs
                elif sheettype == 'time':
                    startcol = 2
                    thesedata = blank2nan(sheetdata.row_values(row, start_colx=startcol,
                                                               end_colx=lastdatacol - 1))  # Data starts in 3rd column, and ends lastdatacol-1
                    assumptiondata = sheetdata.cell_value(row, assumptioncol - 1)
                    if assumptiondata != '':  # There's an assumption entered
                        thesedata = [
                            assumptiondata]  # Replace the (presumably blank) data if a non-blank assumption has been entered
                    data[thispar].append(thesedata)  # Store data
                    checkblank = False if sheetname in skipblanksheets or thispar in skipblankpars else True  # Don't check optional indicators, check everything else
                    validate_data(thesedata, sheetname, thispar, row, checkblank=checkblank,
                                 checkupper=checkupper[thispar], startcol=startcol)

                # It's a matrix, append the data
                elif sheettype == 'matrix':
                    startcol = 2
                    thesedata = sheetdata.row_values(row, start_colx=startcol,
                                                     end_colx=sheetdata.ncols)  # Data starts in 3rd column
                    thesedata = list(map(lambda val: 0 if val == '' else val, thesedata))  # Replace blanks with 0
                    data[thispar].append(thesedata)  # Store data
                    validate_data(thesedata, sheetname, thispar, row, startcol=startcol)

                # It's a constant, create a new dictionary entry
                elif sheettype == 'constant':
                    startcol = 2
                    endcol = 5
                    try:
                        thispar = subparlist.pop(0)  # Get the first item in this list
                    except Exception as E:
                        errormsg = 'Error reading constants sheet: perhaps too many rows?\n'
                        errormsg += '%s' % repr(E)
                        errormsg += versioncheck
                        raise Exception(errormsg)
                    thesedata = blank2nan(sheetdata.row_values(row, start_colx=startcol,
                                                               end_colx=endcol))  # Data starts in 3rd column, finishes in 5th column
                    validatedata(thesedata, sheetname, thispar, row)
                    data[thispar] = thesedata  # Store data

                # It's not recognized: throw an error
                else:
                    errormsg = 'Sheet type "%s" not recognized: please do not change the names of the sheets!' % sheettype
                    raise Exception(errormsg)

    # Check that matrices have correct shape
    data['npops'] = len(data['pops']['short'])
    for key in sheets['Partnerships & transitions']:
        thesedata = data[key]
        matrixshape = shape(array(thesedata))
        correctfirstdim = data['npops'] if key != 'birthtransit' else sum(data['pops']['female'])
        correctseconddim = data['npops']
        if matrixshape[0] != correctfirstdim or matrixshape[1] != correctseconddim:
            errormsg = f'Matrix {key} in partnerships & transitions sheet is not the correct shape'
            errormsg += f'(rows = {matrixshape[0]}, columns = {matrixshape[1]}, should be {correctfirstdim} and {correctseconddim})\n'
            errormsg += 'Check for missing rows or added text'
            raise Exception(errormsg)

    # Store tuples of partnerships
    popkeys = data['pops']['short']
    for row in range(data['npops']):
        for col in range(data['npops']):
            if data['part'][row][col]: data['pships'].append((popkeys[row], popkeys[col]))

    # Do a final check
    failedtopopulate = sc.odict()
    for sheetname, sheetpars in sheets.items():
        if len(sheetpars) > 0:
            failedtopopulate[sheetname] = sheetpars
    if len(failedtopopulate):
        errormsg = 'Not all parameters were successfully populated; missing parameters were:'
        errormsg += f'\n{failedtopopulate}'
        raise Exception(errormsg)

    return data



