# Loopo Toolkit
A Python toolbox for parsing, manipulating, and analyzing Molecular Dynamics simulation data.  

It provides:
- A **core library** (`loopo/`) with classes for configurations, chains, molecules, correlators, and radial profiles.  
- A set of **command-line utilities** (`utils/`) for configurations analysis (mean squared displacement, looping over trajectories, etc.).  

---

## 🚀 Installation

Clone the repository and install it, if you prefer in a virtual environment:
```bash
git clone git@github.com:l09rin/Loopo.git
cd Loopo
python3 -m venv .loopo-venv
source .loopo-venv/bin/activate
pip install -U .
```

To include full analysis features, install further optional dependencies `freud-analysis` and `ovito`:
```bash
python3 -m pip install -U ovito freud-analysis
```

---

## 📂 Project Structure

```
src/loopo/              # Core library
    ├── configuration.py
    ├── chains.py
    ├── bonds.py
    ├── correlation.py
    ├── molecule.py
    ├── profiles.py
src/loopo_actions/      # action classes used by loopOnTraj
utils/                  # Command-line scripts
    ├── loopOnTraj.py
    ├── msd.py
    ├── msd_parallel.py
    └── generate_tsteps_cycle.py
```

---

## 🔧 Command-line Utilities

Once installed, the following commands are available in your shell:

- `loopOnTraj` — versatile tool to parse, manipulation and analyse particle configurations (compatible with various formats).  
- `msd` — compute mean squared displacement (serial).
- `msd-parallel` — compute mean squared displacement (thread-parallelized).
- `generate-tsteps-cycle` — generate a file containing custom sequence of timesteps, also cycles with logarithmically spaced intervals.

### Example:

```bash
loopOnTraj -confs "dump.traj.*" lmp -gyration mol 3 files rg_vs_time.dat average_rg.dat
msd-parallel -confs list confs_list.dat -format xyz -parallel 4
```

---

## 📘 Core Library

The `loopo` package defines the following Python classes:

- **`configuration`** — main container for a single configuration (simulation box, particles with their attributes, etc.).
- **`chains`** — to store and manipulate polymer chains and bonded structures.
- **`bonds`** — bonds parsing and manipulation.
- **`molecule`** — molecule-based wrapper.
- **`correlation`** — defines functor classes to represent and calculate averaged correlation functions and time-dependent quantities.
- **`profiles`** — density and other radial profiles.

These classes can be imported and used directly in Python:

```python
from loopo.configuration import CONFIGURATION
conf = CONFIGURATION.smart_auto_read("traj.dump", "lmp")
```

---

## 🛠️ `loopOnTraj` Utility

This program selects and manipulates particles over configurations dumped by Molecular Dynamics simulation softwares, such as [**LAMMPS**](https://github.com/lammps/lammps), [**MoDyCa**](https://github.com/l09rin/MoDyCa), [**EDMD**](https://github.com/l09rin/EDMD) code.
It supports many type of **actions** such as changing attributes, make selections, computing physical and geometrical quantities and more.

### Usage

```bash
loopOnTraj -confs single_conf.dat { -parallel <nthreads> } { <action> <options> }
loopOnTraj -confs "dump.traj.*" { -parallel <nthreads> } { <action> <options> }
loopOnTraj -confs list file_names.dat { -parallel <nthreads> } { <action> <options> }
```

Supported formats: `lmp`, `xyz`, `sph`.
Configurations can be passed as:
- a single file,
- a glob pattern (`"dump.traj.*"`),
- or a file containing a list of configuration files.

Run `loop_on_configs.py --help` for a comprehensive list of all functionalities.

### Examples

- **Compute gyration radius**:
```bash
loopOnTraj -confs "mgel.dump.*" -gyration mol 3 files rg_vs_time.dat average_rg.dat
```

- **Recenter configuration on molecule's center of mass**:
```bash
loopOnTraj -confs "mgel.dump.*" -recenter mol 1
```

- **Radial, cylindrical, slice density profile**:
```bash
loopOnTraj -confs "mgel.dump.*" -profile sphere 0.5 specie mol 3 center mol 3 file density_profile_mol3.dat
loopOnTraj -confs "mgel.dump.*" -profile cylinder z 0.5 0.5 specie type 2 center fixed 0,0,0 file density_profile_type2.dat
loopOnTraj -confs "mgel.dump.*" -profile linear x 0.5 specie type 2 center fixed 0,0,0 file density_profile_type2.dat
```

- **Convex hull analysis**:
```bash
loopOnTraj -confs "mgel.dump.*" -convex_hull mol 3 file convex_hulls.dat
```

- **Surface mesh analysis**:
```bash
loopOnTraj -confs "mgel.dump.*" -surface_mesh mol 3 radius 5.0 file surface_meshes.dat
```

- **Select/remove particles**:
```bash
loopOnTraj -confs "mgel.dump.*" -sel mol 3
loopOnTraj -confs "mgel.dump.*" -rm type 1
```

- **Change particles attributes**:
Several attributes of particles can be modified, such as charge, mass, type or molecule indices:
```bash
loop_on_configs.py -confs "mgel.dump.*" -change type 3 q -1.0
loop_on_configs.py -confs "mgel.dump.*" -change IDfile filename.dat q none
```

- **Dump and convert formats**:
```bash
loop_on_configs.py -confs "mgel.dump.*" -print output.dat lmp pos:id:type:q mol 3
```


For the full set of actions and options, run:
```bash
loopOnTraj --help
```

---

## 📊 `msd` / `msd-parallel` Utilities

Computes the **mean squared displacement (MSD)** of a system of particles.  
Accepts input files in `lmp`, `sph`, `xyz` formats (default = xyz).  

### Format (xyz example)
```
# N <number_of_particles>
# step <time_step>
x1 y1 z1
x2 y2 z2
...
(blank line)
(repeat for next step)
```

### Usage

```bash
msd-parallel -confs list confs_list.dat -format xyz -parallel 24 -out msd.dat
```

Optional filters:
- `-select` / `-remove` (by molecule, type, charge, crosslinkers)
- `-recenter mol <id>`
- `-sort_id` (if IDs are unsorted)

---

## 📂 Input Formats

Loopo supports multiple configuration formats:
- **LAMMPS dumps (`lmp`)**
- **XYZ trajectories (`xyz`)**
- **SPH (`sph`)**
- **Patchy particles format (`patch`)** containing information on the (2D) angle and attractive patches

---

## 📤 Outputs

Depending on the action, Loopo produces:
- **Trajectory dumps** in `xyz`, `lmp`, `sph` formats.
- **Analysis results** such as:
  - Gyration radius vs time
  - Center of mass trajectories
  - Convex hull / surface mesh geometrical descriptors
  - Density profiles
  - Mean squared displacement

All files are plain-text and suitable for post-processing or visualization.

---

## ⚡ Optional Dependencies

Some features require additional libraries:  

- [`freud-analysis`](https://freud.readthedocs.io/) — order parameters, local environments.  
- [`ovito`](https://www.ovito.org/python/) — surface meshes, advanced visualization/analysis.  

Install them via:
```bash
pip install -U freud-analysis ovito
```

If missing, Loopo will warn you, and the corresponding methods will be disabled.

---

## 📜 License
GPLv3.0 License.

---
