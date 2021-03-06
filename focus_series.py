import os
import numpy as np
import math
import pandas
import datetime as dt
import matplotlib.pyplot as plt
from astropy.io import ascii
from sklearn.metrics import mutual_info_score
import seaborn as sns
sns.set(context='paper', font='monospace')

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
    r"""Read all the tsi files in the directory.
    
    Reads all the tsi files and stores lists of their
    values in a dictionary.
    
    Parameters
    ----------
    tsi_path : str
        Path to directory to search for tsi files.
    limit : int, optional
        A limit to the number of files to load.
    
    Returns
    -------
    tsi : dict
        A dictionary containing lists of all the tsi
        values read.
    """
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
    r"""Read the the gridded psf files in the directory.
    
    Reads all the median values per grid as generated
    by `populate_psf_grid`.
    
    Parameters
    ----------
    med_path : str
        Path to the directory containing the median grids.
    headers : array_like
        A list of the psf headers to look for.
    grid_x, grid_y : int, optional
        Dimension of the grid.
    
    Returns
    -------
    psf : list
        List object containing the median values and their
        position within the grid.
    
    See Also
    --------
    populate_psf_grid : Calculates median values per grid
    """
    psf = []
    for h in headers:
        for x in range(grid_x):
            for y in range(grid_y):
                f = ascii.read(med_path + "grid_" + str(x) + "_" + str(y) + "_" + h + ".psf")
                psf.append([f[h], h, x, y])
    return psf

def read_all_psf_temp_elev(headers, psf_path, tsi_path, limit = None):
    r"""Read all the psf values and relate them to elevation and temperature.
    
    The averages for the psf values of each file as well as the
    corresponding temperature and elevation are stored for each
    header in a dictionary.
    
    Parameters
    ----------
    headers : array_like
        List of the headers to read.
    psf_path : str
        Path to the directory containing the psf files.
    tsi_path : str
        Path to the directory containing the tsi meta-data files.
    limit : int, optional
        A limit to the number of files to load. Default is `None`
        (no limit).
    
    Returns
    -------
    dict
        Each header is used as a key to store a dictionary
        object. These dictionaries contain lists, where each
        row corresponds to a psf file.
    """
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
    
    Parameters
    ----------
    val : str
        Name of the tsi field to read.
    tsi_path : str
        Path to the directory containing the tsi meta-data files.
    limit : int, optional
        Limit to the number of files to load. Default is `None`
        (no limit).
    
    Returns
    -------
    dict
        Dictionary object with three fields, each a list of
        values (temperature, elevation and the tsi value).
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

def get_headers(filename):
    """Get the headers of the given file.
    
    Reads the headers of the given file, using the default
    astropy.io.ascii behavior.
    
    Parameters
    ----------
    filename : str
        Filename to be read.
    
    Returns
    -------
    array_like
        The names of the headers found in the file.
    """
    f = ascii.read(filename)
    return f.colnames

def populate_psf_grid(psf_path, grid_x=3, grid_y=3, max_x=9000, max_y=9000, 
                      headers=None, limit=None):
    r"""Populates a grid with median values of all the psf files.
    
    Reads all the psf files in the given directory and assigns
    each psf parameter as specified by `headers` to a grid.
    For each of these grids, a median value is calculated,
    based on the psf values within the grid.
    
    Parameters
    ----------
    psf_path : str
        Path to the directory containing the psf files.
    grid_x, grid_y : int, optional
        Number of grids. Default is a :math:`3\times 3` grid.
    max_x, max_y : int
        The maximum `x` and `y` pixel coordinates of the psf.
        Default is 9000.
    headers : array_like, optional
        A list of the psf parameters to read. If set to `None`,
        all the parameters will be read.
    limit : int, optional
        Maximum number of psf files to read. Default is `None`.
        If set to `None` all the files in the directory will be
        read.
    
    Returns
    -------
    list
        A two dimensional list containing median values
        for each grid and header.
    """
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
        if (limit is not None and i >= limit):
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
    r"""Populate a grid with average values from the dictionary.

    Takes values with temperature and elevation information
    and assigns them to a grid position. For each grid averages
    of the values within the grid are calculated.
    
    Parameters
    ----------
    dic : dict
        Dictionary containing lists of the values and their
        temperature and elevation.
    grid_temp, grid_elev : int
        Dimension of the grid.
    single_value : bool
        Used to determine, whether there is a single value, as
        is the case with tsi data or an average of data, as is
        the case for psf data. Default is `True`.
    
    Returns
    -------
    dict
        A dictionary of lists as defined by `generate_temp_elev_dict`.
    """
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

def get_sparse_tsi(params, n = 2):
    """ Returns a list of shortened tsi parameters.
    
    Takes a list of tsi parameters and shortens
    their name.
    
    Parameters
    ----------
    params : array_like
        List of tsi parameters.
    n : int, optional
        The number of words to include in the shortened
        name.
    
    Returns
    -------
    list
        List of the same parameters with shortened names.
    """
    sparse = []
    for p in params:
        spl = p.split('.')
        pp = ''
        for i in range(1, min(len(spl), n)+1):
            pp = spl[-i] + ('.' if i > 1 else '') + pp
        sparse.append(pp)
    return sparse

def mutual_info_heatmap(file_src, plot_path, f_num, size=None,  
                        key1='key1', key2='key2', keymi='mutual_info', 
                        xlabel_callback=lambda l: l,
                        ylabel_callback=lambda l: l):
    """ Plots a heatmap of the mutual information.
    
    Creates a correlation matrix from the mutual information
    in a file and plots a heatmap. Will also automatically
    divide the matrix plot into a grid and save each grid
    separately.
    
    Parameters
    ----------
    file_src : str
        Filename to read.
    plot_path : str
        Path to directory where the plots will be saved.
    f_num : int
        The maximum number of parameters that will be
        displayed on the axes. The plot will be divided
        into a grid accordingly.
    size : float, optional
        Size of the plot in inches. If `None`, will be
        set automatically.
    key1, key2 : str, optional
        Name of the columns in the file containing the
        parameter names which where compared.
    keymi : str, optional
        Column name containing the mutual information value.
    xlabel_callback, ylabel_callback : callable, optional
        These callback functions will be called with the list
        of parameters. The return value will be used for the
        axis labels.
    """
    mi = ascii.read(file_src)
    fields = sorted(set(mi[key1]))
    fields2 = sorted(set(mi[key2]))
    f_len  = len(fields)
    f_len2 = len(fields2)
    print(len(set(fields+fields2)), 'parameters total')
    for i in range(math.ceil(f_len/f_num)):
        start = i*f_num
        end   = min((i+1)*f_num,f_len)
        sub_fields = fields[start:end]
    
        for j in range(math.ceil(f_len2/f_num)):
            start2 = j*f_num
            end2   = min((j+1)*f_num,f_len2)
            sub_fields2 = fields2[start2:end2]
            mat = {k : {k : 0 for k in sub_fields2} for k in sub_fields}   
            for r in mi:
                if r[key1] in sub_fields and r[key2] in sub_fields2:
                    mat[r[key1]][r[key2]] = float(r[keymi])
        
            plt.ioff()
            f, ax = plt.subplots(figsize=size)
            df = pandas.DataFrame(mat)
            sns.heatmap(
                df, square=False,vmin=0,vmax=1,
                xticklabels=xlabel_callback(df.axes[1]),
                yticklabels=ylabel_callback(df.axes[0])
            )
            plt.yticks(rotation='horizontal')
            plt.xticks(rotation='vertical')
            f.tight_layout()
            f.savefig(
                plot_path + 'heatmap_' + str(i) + '-' + str(j) + '.pdf'
            )
            plt.close()