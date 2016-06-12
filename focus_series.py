import os
import numpy as np
import math
import datetime as dt
from astropy.io import ascii
from sklearn.metrics import mutual_info_score

__project__  = 'focus-series'
__author__   = 'Jean Elsner'
__version__  = '1.0.0'
__email__    = 'jean.elsner@googlemail.com'
__homepage__ = 'https://github.com/JeanElsner/focus-series'

ELEV = 'CURRENT.OBJECT.INSTRUMENTAL.ZD'
TEMP = 'AUXILIARY.SENSOR[8].VALUE'
FILTER_BESTFOC = True

def _calc_MI(x, y, bins=None):
    if bins is None:
        bins = 25
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi

def mutual_info(X, Y, bins=None):
    r"""Calculcates the mutual information between two discrete distributions.
    
    Parameters
    ----------
    X, Y : array_like
        Discrete distributions to compare.
    bins : int, optional
        Number of bins for the distributions.
        
    Returns
    -------
    float
        Mutual information score between the two distributions.
    """
    return _calc_MI(X, Y, bins=bins) / np.sqrt(_calc_MI(X, X, bins=bins) * _calc_MI(Y, Y, bins=bins))

def get_grid_position(x, y, dim_x, dim_y, bounds):
    r""" Gives the position of point on a grid
    
    Devides the domain as given by the boundaries into
    a grid of dimension :math:`\mathtt{dim_x}\times\mathtt{dim_y}`
    and gives the grid-coordinates the point falls into.
    
    Parameters
    ----------
    x, y : float
        Coordinates of a point within the domain.
    dim_x, dim_y : int
        Dimensions of the 2d grid.
    bounds : array_like
        An array_like object of dimension (2,2)
        describing the boundaries of the domain,
        e.g. ((min_x, max_x), (min_y, max_x))
    
    Returns
    -------
    grid_x, grid_y : int
        Coordinates of the tile the point falls into.
    """
    grid_x = min(math.floor(abs(x - bounds[0][0]) / 
                        abs(bounds[0][0] - bounds[0][1]) * dim_x), dim_x-1)
    grid_y = min(math.floor(abs(y - bounds[1][0]) / 
                        abs(bounds[1][0] - bounds[1][1]) * dim_y), dim_y-1)
    return (grid_x, grid_y)

def reduce_to_grid(x, y, z, dim_x=4, dim_y=4, fallback=lambda z,x,y: min(z)):
    r"""Reduces a surface to a grid comprised of median values.
    
    Reduces a surface to a grid of knots. The knots are calculated
    using the median of the surface points in a rectangle surrounding
    the knot point. The rectangles are defined by dividing the domain
    as given by :math:`x\otimes y` into :math:`\mathtt{dim_x}\times\mathtt{dim_y}`.
    
    Parameters
    ----------
    x, y, z : array_like
        A surface defined by :math:`z=f(x, y)`.
    dim_x, dim_y: int, optional
        The number of knots to calculate for a given axis.
    fallback : callable, optional
        A fallback function that returns a value for a knot,
        if there is no data in the corresponding rectangle.
        Receives the array_like `z` as well as the `x` and `y`
        position of the knot as arguments.
        
    Returns
    -------
    grid_means : list
        Two dimensional list containing the knot values.
    knots : tuple of list
        A tuple containing two lists of the knot `x` and `y` coordinates.
    """
    grid = [[[] for _y in range(dim_y)] for _x in range(dim_x)]
    knots = ([], [])
    bounds = ((min(x), max(x)), (min(y), max(y)))
    for i in range(len(z)):
        grid_x, grid_y = get_grid_position(x[i], y[i], dim_x, dim_y, bounds)
        grid[grid_x][grid_y].append(z[i])
    grid_mean = [[0 for _y in range(dim_y)] for _x in range(dim_x)]
    for _y in range(dim_y):
        knots[1].append(bounds[1][0] + abs(bounds[1][0]-bounds[1][1])*((_y+.5)/dim_y))
    for _x in range(dim_x):
        knots[0].append(bounds[0][0] + abs(bounds[0][0]-bounds[0][1])*((_x+.5)/dim_x))
        for _y in range(dim_y):
            mean = np.mean(grid[_x][_y])
            grid_mean[_x][_y] = mean if not np.isnan(mean) else fallback(z, knots[0][_x], knots[1][_y])
    return (grid_mean, knots)

def get_tsi_filenames(tsi_path):
    r"""Get all the tsi meta-data filenames in the directory.
    
    Parameters
    ----------
    tsi_path : str
        Path to search for tsi files.
    
    Returns
    -------
    tsi : list
        A list of all the filenames.
    """
    tsi =[]
    for _file in os.listdir(tsi_path):
        if FILTER_BESTFOC and not "bestfoc" in _file:
            continue
        if _file.endswith(".tab"):
            tsi.append(_file)
    return tsi

def get_psf_filenames(psf_path):
    r"""Get all the psf filenames in the directory.
    
    Parameters
    ----------
    psf_path : str
        Path to search for psf files.
    
    Returns
    -------
    psf : list
        A list of all the filenames.
    """
    psf =[]
    for _file in os.listdir(psf_path):
        if FILTER_BESTFOC and not "bestfoc" in _file:
            continue
        if _file.endswith(".psf"):
            psf.append(_file)
    return psf

def generate_grid_list(x, y, headers):
    r"""Generate a two dimensional list object populated with dictionaries.
    
    The two dimensional list contains dictionaries with
    list object as an entries for the keys given in `headers`.
    
    Parameters
    ----------
    x, y : int
        The length of the x and y dimensions of the list.
    headers : array_like
        A list of headers to use as keys in the dictionaries.
    
    Returns
    -------
    lst : list
        Two dimensional as described above.
    """
    lst = []
    for _x in range(x):
        lst.append([])
        for _y in range(y):
            lst[_x].append({})
            for h in headers:
                lst[_x][_y][h] = []
    return lst

def generate_temp_elev_header_dict(headers):
    r"""Generate a dictionary with dictionary item for each header.
    
    The items are dictionaries as generated by `generate_temp_elev_dict`.
    
    Parameters
    ----------
    headers : array_like
        The headers to be used as keys in the dictionary.
    
    Returns
    -------
    dict
        Dictionary containing a dictionary for each header.
    
    See Also
    --------
    generate_temp_elev_dict : Entries are generated using this routine.
    """
    lst = {}
    for h in headers:
        lst[h] = generate_temp_elev_dict()
    return lst

def generate_temp_elev_dict():
    r"""Generate a dictionary with list items for statistical data.
    
    Returns
    -------
    dict
        A dictionary with fields for statistical values as
        well as temperature and elevation.
    """
    return {'temp' : [], 'elev' : [], 'mean' : [], 'med' : [], 
            'std' : [], 'var' : [], 'err' : []}

def generate_temp_elev_val_dict():
    r"""Generate a dictionary with list items for values.
    
    Returns
    -------
    dict
        Dictionary with fields `temp`, `elev` and `val`.
        Each entry has a an empty list object as value.
    """
    return {'temp' : [], 'elev' : [], 'val' : []}

def read_tsi(f):
    r"""Reads a single tsi meta-data file.
    
    Parameters
    ----------
    f : str
        Filename to read.
        
    Returns
    -------
    tab : dict
        A dictionary of key, value pairs.
    """
    tab = {}
    with open(f, 'r') as _f:
        for r in _f.readlines():
            k, v = r.split("=")
            v = v.split("#")[0]
            if is_number(v):
                tab[k] = float(v)
    return tab

def is_number(s):
    r"""Check whether a given value is a number.
    
    Checks whether a given value can be interpreted as a float.
    
    Parameters
    ----------
    s
        Value to be checked.
        
    Returns
    -------
    bool
        True if value is a number.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def check_if_const(dst, threshold=.9):
    r"""Check whether the given distribution is (nearly) constant.
    
    The discrete distribution has to be binned. However it need
    not be normalised.
    
    Parameters
    ----------
    dst : array_like
        The distribution as a propability density function.
    threshold : float, optional
        If at least the fraction given by `threshold` is contained
        within a single bin, the distribution is considered constant.
    
    Returns
    -------
    bool
        True if distribution is constant.
    """
    if len(dst) <= 1:
        return True
    return max(dst)/float(sum(dst)) > threshold
    
def read_all_tsi(tsi_path, limit=None):
    r"""Read all the tsi files in the directory."""
    tsi = {}
    i = 0
    for f in get_tsi_filenames(tsi_path):
        for k, v in read_tsi(tsi_path + f).items():
            if k in tsi.keys():
                tsi[k].append(v)
            else:
                tsi[k] = [v]
        i += 1
        if (limit is not None and i >= limit):
            break
    return tsi
    
def read_all_psf_grids(med_path, headers, grid_x=3, grid_y=3):
    r"""Read the the gridded psf files in the directory."""
    psf = []
    for h in headers:
        for x in range(grid_x):
            for y in range(grid_y):
                f = ascii.read(med_path + "grid_" + str(x) + "_" + str(y) + "_" + h + ".psf")
                psf.append([f[h], h, x, y])
    return psf
    
def read_all_psf_temp_elev(headers, psf_path, tsi_path, limit = None):
    r"""Read all the psf values and relate them to elevation and temperature."""
    lst = generate_temp_elev_header_dict(headers)
    i = 0
    for f in get_psf_filenames(psf_path):
        i += 1
        t = ascii.read(psf_path + f)
        try:
            tsi = read_tsi(tsi_path + f[:-3] + "tab")

            for h in headers:
                a = np.array(t[h], dtype=float)
                lst[h]['temp'].append(tsi[TEMP])
                lst[h]['elev'].append(tsi[ELEV])
                lst[h]['mean'].append(np.nanmean(a))
                lst[h]['std'].append(np.nanstd(a))
                lst[h]['med'].append(np.nanmedian(a))
                lst[h]['var'].append(np.nanvar(a))
        except FileNotFoundError:
            print("Warning: file", f, "has no tsi data")
        if (limit is not None and i >= limit):
            break
    return lst
    
def read_all_tsi_val_temp_elev(val, tsi_path, limit = None):
    r"""Read all the values for the given tsi header 
    in relation to temperature and elevation.
    """
    lst = generate_temp_elev_val_dict()
    i = 0
    for f in get_tsi_filenames(tsi_path):
        i += 1
        tsi = read_tsi(tsi_path + f)

        lst['temp'].append(tsi[TEMP])
        lst['elev'].append(tsi[ELEV])
        lst['val'].append(tsi[val])
        if (limit is not None and i >= limit):
            break
    return lst

def populate_psf_grid(psf_path, grid_x=3, grid_y=3, max_x=9000, max_y=9000, headers=None, _max=None):
    r"""Populates a grid with median values of all the psf files."""
    files = get_psf_filenames(psf_path)
    if headers is None:
        headers = get_headers(psf_path + files[0])
    med = generate_grid_list(grid_x, grid_y, headers)
    i = 0
    for f in files:
        i = i+1
        psf = generate_grid_list(grid_x, grid_y, headers)
        t = ascii.read(psf_path + f)
        for r in t:
            x_ind = int(math.floor(r['x']/max_x*grid_x))
            y_ind = int(math.floor(r['y']/max_y*grid_y))
            for c in headers:
                if is_number(r[c]):
                    psf[x_ind][y_ind][c].append(r[c])
        for x in range(grid_x):
            for y in range(grid_y):
                for c in psf[x][y]:
                    if c in headers:
                        _med = np.median(psf[x][y][c])
                        if is_number(_med) and not math.isinf(_med):
                            med[x][y][c].append(_med)
        if (_max is not None and i >= _max):
            break
    return med

def _get_temp_elev_dict_ext(dic):
    return (max(dic['temp']), min(dic['temp']),
           max(dic['elev']), min(dic['elev']))

def _get_temp_elev_bound(ext, grid_temp, grid_elev, temp, elev):
    return (ext[1] + (ext[0] - ext[1]) / grid_temp * temp,
            ext[1] + (ext[0] - ext[1]) / grid_temp * (temp + 1),
            ext[3] + (ext[2] - ext[3]) / grid_elev * elev,
            ext[3] + (ext[2] - ext[3]) / grid_elev * (elev + 1))

def populate_temp_elev_dict(dic, grid_temp, grid_elev, single_value=True):
    r"""Populate a grid with average values from the dictionary."""
    ext = _get_temp_elev_dict_ext(dic)
    g = generate_temp_elev_dict()
    for __temp in range(grid_temp):
        for __elev in range(grid_elev):
            bound = _get_temp_elev_bound(ext, grid_temp, grid_elev, __temp, __elev)
            std, mean, med, var, val = [[],[],[],[],[]]
            for k, v in enumerate(dic['temp']):
                if dic['temp'][k] >= bound[0] and dic['temp'][k] <= bound[1]:
                    if dic['elev'][k] >= bound[2] and dic['temp'][k] <= bound[3]:
                        if single_value:
                            val.append(dic['val'][k])
                        else:
                            std.append(dic['std'][k])
                            mean.append(dic['mean'][k])
                            med.append(dic['med'][k])
                            var.append(dic['var'][k])
            g['temp'].append((bound[0] + bound[1])/2)
            g['elev'].append((bound[2] + bound[3])/2)
            if single_value:
                g['std'].append(np.nanstd(val))
                g['mean'].append(np.nanmean(val))
                g['med'].append(np.nanmedian(val))
                g['var'].append(np.nanvar(val))
                g['err'].append(np.nanstd(val)/math.sqrt(len(val)))
            else:
                g['std'].append(np.nanmedian(std))
                g['mean'].append(np.nanmedian(mean))
                g['med'].append(np.nanmedian(med))
                g['var'].append(np.nanmedian(var))
                g['err'].append(np.nanstd(mean)/math.sqrt(len(std)))
    return g

def search_meteo(dt_obj, meteo_path, time_delta = 5):
    r"""Searches meteological log files.
    
    Given a date and time, tries to find the corresponding
    log file and searches it for data within the timeframe.
    
    Parameters
    ----------
    dt_obj : datetime.datetime
        The date and time to search for.
    meteo_path : str
        Path to the directory containing the log files.
    time_delta : float, optional
        Maximum time difference in seconds. Only data
        within the difference is used. Default is the
        interval measurements are taken at (5)
    
    Returns
    -------
    dict
        A dictionary containing all the measurements.
    """
    f = meteo_path + dt_obj.strftime('%y%m%d') + '.log'
    keys, values, __h = ([], [], None)
    try :
        with open(f, 'r') as _f:
            for r in _f:
                #2015-07-01T12:50:58.714407
                try:
                    tpl = r[:26].split('T')
                    t = tpl[1].split(':')
                    d = tpl[0].split('-')
                except IndexError:
                    print("Couldn't parse line:", r)
                    continue
                _dt = dt.datetime(int(d[0]), int(d[1]), int(d[2]), 
                                  int(t[0]), int(t[1]), int(float(t[2])))
                if not dt_obj.tzinfo is None:
                    _dt = dt_obj.tzinfo.localize(_dt)
                if not __h is None and not _dt == __h:
                        break
                delta = (_dt - dt_obj).seconds
                if delta  <= time_delta:
                    __h = _dt
                    tpl = r[27:].split('=')
                    if not tpl[0].strip() in keys:
                        keys.append(tpl[0].strip())
                        values.append(tpl[1].strip())
    except FileNotFoundError:
        print('No meteological data for', dt_obj)
    return dict(zip(keys, values))