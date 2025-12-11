import sys, os, subprocess
from os import system
from os import path

import cctk
import autode as ade

import numpy as np
import scipy as sp
import rdkit
import copy

import multiprocessing, multiprocessing.pool
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import MolDrawing
from rdkit.Chem.Draw import IPythonConsole

import datamol as dm

# import useful_rdkit_utils as uru
import pandas as pd
from tqdm.auto import tqdm
import seaborn as sns
import mols2grid
# import pingouin

import matplotlib.pyplot as plt

from urllib.request import urlopen
from urllib.parse import quote

import ipywidgets
from ipywidgets import interact, fixed
import warnings
warnings.filterwarnings('ignore')
import py3Dmol

# import pygwalker as pyg

import pyscf
# from pyscf.semiempirical import MINDO3  # Module not available in standard PySCF
# Note: For semi-empirical calculations, use xtb (GFN2-xTB) via command-line or AutodE
# xtb binary is available at: /opt/homebrew/Caskroom/miniforge/base/envs/mml_comp_chem/bin/xtb
from rdkit.Chem.Draw.IPythonConsole import drawMol3D

sys.path.append("/etc/opt/orca-5.0.3/")
sys.path.append("/etc/opt/openmpi-4.1.1/bin")

ORCA_PATH = '/etc/opt/orca-5.0.3'
OPEN_MPI_PATH = '/etc/opt/openmpi-4.1.1/'
OPEN_MPI_BIN_PATH = '/etc/opt/openmpi-4.1.1/bin'
OPEN_MPI_LIB_PATH = '/etc/opt/openmpi-4.1.1/lib'

os.environ['PATH'] += os.pathsep + ORCA_PATH
os.environ['PATH'] += os.pathsep + OPEN_MPI_PATH
os.environ['PATH'] += os.pathsep + OPEN_MPI_BIN_PATH

#os.environ['NBOBIN'] = '/opt/nbo7/bin/'
os.environ['NBOEXE'] = '/opt/nbo7/bin/nbo7.i4.exe'
os.environ['GENEXE'] = '/opt/nbo7/bin/gennbo.i4.exe'

old = os.environ.get("LD_LIBRARY_PATH")
if old:
    os.environ["LD_LIBRARY_PATH"] = old + ":" + 'PATH'
else:
    os.environ["LD_LIBRARY_PATH"] = 'PATH'

os.environ['LD_LIBRARY_PATH'] += os.pathsep + ORCA_PATH
os.environ['LD_LIBRARY_PATH'] += os.pathsep + OPEN_MPI_LIB_PATH

os.environ['OMPI_ALLOW_RUN_AS_ROOT'] = '1'
os.environ['OMPI_ALLOW_RUN_AS_ROOT_CONFIRM'] = '1'

ade.Config.n_cores = 6
ade.Config.max_core = 2000

def plotter(x, y, data_label, title, x_label="X", y_label="Y", color="black", marker=None, multiple=False, scatter=False, nothing=False):
    # make figure and axes
    fig, ax = plt.subplots(figsize=(10, 10))
    # aesthetic settings for plot
    font = 'Arial'
    SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 16, 20, 24
    ax.grid(True, linewidth=1.0, color='0.95')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_axisbelow(True)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3.0)
    for tick in ax.get_yticklabels():
        tick.set_fontname(font)
        tick.set_fontsize(SMALL_SIZE)
    for tick in ax.get_xticklabels():
        tick.set_fontname(font)
        tick.set_fontsize(SMALL_SIZE)

    if multiple:
        for i in range(len(y)):
            plt.scatter(x[i], y[i], label = data_label[i], color=color[i], marker=marker, linewidth=2.0)
    elif scatter:
        plt.scatter(x, y, label=data_label, color=color, marker=".", linewidth=2.0)
    elif nothing:
        return fig, ax
    else:
        plt.plot(x, y, label=data_label, color=color, marker=marker, linewidth=2.0)
    ax.legend(fontsize=SMALL_SIZE, loc=0)
    plt.title(title, fontsize=SMALL_SIZE, fontname=font)
    plt.xlabel(f"{x_label}", fontsize=SMALL_SIZE, fontname=font)
    plt.ylabel(f"{y_label}", fontsize=SMALL_SIZE, fontname=font)
    
    return fig, ax

def CIRconvert(ids):
    #try:
    url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
    ans = urlopen(url).read().decode('utf8')
    return ans
    #except:
    #    return 'Did not work'

class HeaderNotFoundError(Exception):
    pass


def jump_to_header(lines, header):
    """
    Given a list of lines, truncate the start of the list so that the first line
    of the new list contains the header.

    Args:
            lines: List of lines.
            header: Substring to match.

    Returns:
            Truncated lines.

    Raises:
            HeaderNotFoundError
    """

    # Search for the header
    for i, line in enumerate(lines):
        if header in line.strip():
            return lines[i:]

    # Search failed
    raise HeaderNotFoundError(f"Header {header} could not be found in the lines.")


# Orbitals
orbitals = ["s", "p", "d", "f"]


def get_percentage(line, orbital):
    """
    Retrieve the percent character of an orbital.

    Args:
            line: Line containing orbital and percentage.
            orbital: Type of orbital (s, p, d, f).

    Returns:
            Percentage of character.

    Raises:
            n/a
    """

    # Locate orbital in line
    index = line.find(orbital)
    line = line[index:]

    # Locate the first open bracket
    index = line.find("(")
    line = line[index:]

    # Isolate the percentage
    return line[1:7].strip()


def z_int(string):
    """
    Convert string to integer.
    If string empty, return -1.

    Args:
            string: Input to be cast to int.

    Returns:
            Int representation.

    Raises:
            n/a
    """
    try:
        return int(string)
    except ValueError:
        return -1


def parse_natural_populations(lines) -> pd.DataFrame:
    """
    Parse the natural populations section of NBO output.

    Args:
            lines: Gaussian output lines.

    Returns:
            Data frame of formatted output.

    Raises:
            HeaderNotFoundError
    """

    # Natural populations
    lines = jump_to_header(lines, "Summary of Natural Population Analysis:")

    # Jump to column names
    lines = lines[4:]
    columns = lines[0].split()

    # Jump to values
    lines = lines[2:]
    data = []
    for line in lines:

        # Termination condition
        if "=" in line:
            break

        # Extract the values
        values = line.split()
        data.append(
            [
                str(values[0]),
                int(values[1]),
                float(values[2]),
                float(values[3]),
                float(values[4]),
                float(values[5]),
                float(values[6]),
            ]
        )

    # Store values in a dataframe
    pop_df = pd.DataFrame(data=data, columns=columns)
    return pop_df


def parse_hybridization_character(lines):
    """
    Parse the hybridization character section of NBO output.

    Args:
            lines: Gaussian output lines.

    Returns:
            Data frames of formatted output.

    Raises:
            HeaderNotFoundError
    """

    # NBO Analysis
    # import ipdb; ipdb.set_trace()
    lines = jump_to_header(lines, "(Occupancy)   Bond orbital / Coefficients / Hybrids")

    # Jump to values
    lines = lines[2:]

    # Save the data for different types of orbitals
    lp_data = []
    bd_data = []

    # Iterate over the lines
    i = -1
    while True:
        i += 1
        line = lines[i]

        # Termination condition
        if "NHO DIRECTIONALITY AND BOND BENDING" in line:
            break

        # Lone pair
        if "LP" in line or "LV" in line:
            entry = {orbital: 0.0 for orbital in orbitals}
            entry["bond index"] = line[0:4].strip()
            entry["occupancy"] = line[7:14].strip()
            entry["type"] = line[16:19].strip()
            entry["orbital index"] = line[20:22].strip()
            entry["atom symbol"] = line[23:25].strip()
            entry["atom number"] = line[25:28].strip()

            # Populate the orbital percentages
            for orbital in orbitals:
                if orbital in line:
                    entry[orbital] = get_percentage(line, orbital)

            # Move one line down
            i += 1
            line = lines[i]

            # Populate the orbital percentages
            for orbital in orbitals:
                if orbital in line:
                    entry[orbital] = get_percentage(line, orbital)

            # Save the entry
            lp_data.append(entry)

        # Bonding
        if "BD" in line:
            entry = {
                f"atom {i} {orbital}": 0.0 for orbital in orbitals for i in range(1, 3)
            }
            entry["bond index"] = line[0:4].strip()
            entry["occupancy"] = line[7:14].strip()
            entry["type"] = line[16:19].strip()
            entry["orbital index"] = line[20:22].strip()
            entry["atom 1 symbol"] = line[23:25].strip()
            entry["atom 1 number"] = line[25:28].strip()
            entry["atom 2 symbol"] = line[29:31].strip()
            entry["atom 2 number"] = line[31:34].strip()

            # Move one line down
            i += 1
            line = lines[i]

            entry["atom 1 polarization"] = line[16:22].strip()
            entry["atom 1 pol coeff"] = line[24:33].strip()

            # Populate the orbital percentages
            for orbital in orbitals:
                if orbital in line:
                    entry[f"atom 1 {orbital}"] = get_percentage(line, orbital)

            # Move one line down
            i += 1
            line = lines[i]

            # Populate the orbital percentages
            for orbital in orbitals:
                if orbital in line:
                    entry[f"atom 1 {orbital}"] = get_percentage(line, orbital)

            # Move down until you see an orbital
            while "s" not in line:
                i += 1
                line = lines[i]

            entry["atom 2 polarization"] = line[16:22].strip()
            entry["atom 2 pol coeff"] = line[24:33].strip()

            # Populate the orbital percentages
            for orbital in orbitals:
                if orbital in line:
                    entry[f"atom 2 {orbital}"] = get_percentage(line, orbital)

            # Move one line down
            i += 1
            line = lines[i]

            # Populate the orbital percentages
            for orbital in orbitals:
                if orbital in line:
                    entry[f"atom 2 {orbital}"] = get_percentage(line, orbital)

            # Save the entry
            bd_data.append(entry)

    # Store values in a dataframe
    lp_df = pd.DataFrame(data=lp_data)
    bd_df = pd.DataFrame(data=bd_data)

    return [lp_df, bd_df]


def parse_perturbation_energy(lines):
    """
    Parse the perturbation energy section of NBO output.

    Args:
            lines: Gaussian output lines.

    Returns:
            Data frame of formatted output.

    Raises:
            HeaderNotFoundError
    """

    # 2nd order perturbation theory analysis
    lines = jump_to_header(
        lines, "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS"
    )

    # Jump to values
    i = -1
    while True:
        i += 1
        line = lines[i]
        if "within" in line:
            lines = lines[i:]
            break

    # Extract 2nd order data
    e2_data = []
    for line in lines:

        # Termination condition
        if "NATURAL BOND ORBITALS" in line:
            break

        # Skip conditions
        if "" == line.strip():
            continue
        if "unit" in line:
            continue
        if "None" in line:
            continue
        if "RY" in line:
            continue

        # Extract the values
        entry = {}
        line = line[3:-2]
        chars = "-(). "

        entry["donor bond index"] = int(line[0:3].strip(chars))
        entry["donor type"] = str(line[4:8].strip(chars))
        entry["donor orbital index"] = int(line[7:11].strip(chars))
        entry["donor atom 1 symbol"] = str(line[12].strip(chars))
        entry["donor atom 1 number"] = int(line[13:16].strip(chars))
        entry["donor atom 2 symbol"] = str(line[17:19].strip(chars))
        entry["donor atom 2 number"] = z_int(line[19:22].strip(chars))
        entry["acceptor bond index"] = int(line[24:30].strip(chars))
        entry["acceptor type"] = str(line[31:35].strip(chars))
        entry["acceptor orbital index"] = int(line[36:38].strip(chars))
        entry["acceptor atom 1 symbol"] = str(line[39:41].strip(chars))
        entry["acceptor atom 1 number"] = int(line[41:44].strip(chars))
        entry["acceptor atom 2 symbol"] = str(line[44:47].strip(chars))
        entry["acceptor atom 2 number"] = z_int(line[48:50].strip(chars))
        pe = line[49:51].strip(chars).split()
        try: 
            float(pe[0])
            entry["perturbation energy"] = float(line[49:51].strip(chars))
        except:
            entry["perturbation energy"] = np.nan
        entry["energy difference"] = float(line[61:69].strip(chars))
        entry["fock matrix element"] = float(line[69:78].strip(chars))
        e2_data.append(entry)

    # Store values in a dataframe
    e2_df = pd.DataFrame(data=e2_data)
    return e2_df


def nbo_parser(fname):
    """
    Parse all the important sections of NBO output.

    Args:
            fname: Path to gaussian NBO output.

    Returns:
            Data frames of formatted output.

    Raises:
            HeaderNotFoundError
    """

    # Open the lines
    with open(fname, "r") as f:
        lines = f.readlines()

    # Compile the dataframes
    dfs = {}
    dfs["natural_populations"] = parse_natural_populations(lines)
    dfs["hybridization_character"] = parse_hybridization_character(lines)
    dfs["perturbation_energy"] = parse_perturbation_energy(lines)
    return dfs
    
def MolTo3DView(xyz, size=(300, 300), style="stick", surface=False, opacity=0.5):
    """Draw molecule in 3D
    
    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'carton'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    assert style in ('line', 'stick', 'sphere', 'carton')
    # mblock = Chem.MolToMolBlock(mol)
    # xyz = mblock
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(open(xyz, 'r').read(), 'xyz')
    viewer.setStyle({style:{}})
    viewer.setBackgroundColor('0xeeeeee')

    viewer.setHoverable({},True,'''function(atom,viewer,event,container) {
                   if(!atom.label) {
                    atom.label = viewer.addLabel(atom.atom+":"+atom.serial,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                   }}''',
               '''function(atom,viewer) { 
                   if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                   }
                }''')


    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    return viewer

#################################################################
#                                                               #
#                      06-310 Utilities                         #
#                                                               #
#################################################################
"""
def vibration(xyz, frequency, size=(300, 300), style="stick", surface=False, opacity=0.5):
    assert style in ('line', 'stick', 'sphere', 'carton')
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(open(xyz, 'r').read(), 'xyz', {'vibrate': {'frames':10,'amplitude':frequency}})
    viewer.setStyle({style:{}})
    viewer.setBackgroundColor('0xeeeeee')

    viewer.setHoverable({},True,'''function(atom,viewer,event,container) {
                   if(!atom.label) {
                    atom.label = viewer.addLabel(atom.atom+":"+atom.serial,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                   }}''',
               '''function(atom,viewer) { 
                   if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                   }
                }''')


    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    return viewer

def view_scan(directory):
    confs = [os.path.join(directory, xyz) for xyz in os.listdir(directory) if xyz.endswith("xyz")]
    numbers = [float(i.split("/")[1].split("_")[2]) for i in confs]
    zipped = zip(confs, numbers)
    sorted_result = sorted(zipped, key = lambda x: x[1])
    confs, numbers = zip(*sorted_result)
    return confs, numbers

def structure_viewer(idx):
    mol = structures[idx]
    return MolTo3DView(mol).show()

def get_normal_mode(molecule, normal_mode):
    elements = molecule.atomic_symbols
    coords = np.array(molecule.coordinates) # To transform from au to A
    natm = molecule.n_atoms
    vib_xyz = "%d\n\n" % natm
    nm = normal_mode.reshape(natm, 3)
    for i in range(natm):
        # add coordinates:
        vib_xyz += elements[i] + " %15.7f %15.7f %15.7f " % (coords[i,0], coords[i,1], coords[i,2])
        # add displacements:
        vib_xyz += "%15.7f %15.7f %15.7f\n" % (nm[i,0], nm[i,1], nm[i,2])
    return vib_xyz

def plotter(x, y, data_label, title, x_label="X", y_label="Y", color="black", marker=None, multiple=False, scatter=False):
    # make figure and axes
    fig, ax = plt.subplots(figsize=(12, 8))
    # aesthetic settings for plot
    font = 'Arial'
    SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 16, 20, 24
    ax.grid(True, linewidth=1.0, color='0.95')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_axisbelow(True)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3.0)
    for tick in ax.get_yticklabels():
        tick.set_fontname(font)
        tick.set_fontsize(SMALL_SIZE)
    for tick in ax.get_xticklabels():
        tick.set_fontname(font)
        tick.set_fontsize(SMALL_SIZE)

    if multiple:
        for i in range(len(y)):
            plt.plot(x, y[i], label = data_label[i], color=color[i], marker=marker, linewidth=2.0)
    elif scatter:
        plt.scatter(x, y, label=data_label, color=color, marker=".", linewidth=2.0)
    else:
        plt.plot(x, y, label=data_label, color=color, marker=marker, linewidth=2.0)
    ax.legend(fontsize=SMALL_SIZE, loc=0)
    plt.title(title, fontsize=SMALL_SIZE, fontname=font)
    plt.xlabel(f"{x_label}", fontsize=SMALL_SIZE, fontname=font)
    plt.ylabel(f"{y_label}", fontsize=SMALL_SIZE, fontname=font)
    
    return fig, ax
"""