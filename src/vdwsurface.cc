//  Copyright (c) 2014, Robert McGibbon
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public License
//  as published by the Free Software Foundation; either version 2.1
//  of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////////
// Calculation of the fused-sphere van der waals surface of a molecule.
//////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <map>
#include <set>

#include "dotsphere.h"
#include "vdwsurface.h"
#include "Vec3.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace OpenMM;


// A. Bondi (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68: 441. doi:10.1021/j100785a001
static map<string, double> create_bondi_radii() {
    map<string, double> m;
    m["H"] = 1.1;
    m["Li"] = 2.2;
    m["Be"] = 1.9;
    m["B"] = 1.8;
    m["C"] = 1.7;
    m["N"] = 1.55;
    m["O"] = 1.52;
    m["F"] = 1.47;
    m["Na"] = 2.4;
    m["Mg"] = 2.2;
    m["Al"] = 2.1;
    m["Si"] = 2.1;
    m["P"] = 1.8;
    m["S"] = 1.8;
    m["Cl"] = 1.75;
    m["K"] = 2.8;
    m["Ca"] = 2.31;
    m["Sc"] = 2.3;
    m["Ti"] = 2.15;
    m["V"] = 2.05;
    m["Cr"] = 2.05;
    m["Mn"] = 2.05;
    m["Fe"] = 2.05;
    m["Co"] = 2.0;
    m["Ni"] = 2.0;
    m["Cu"] = 2.0;
    m["Zn"] = 2.1;
    m["Ga"] = 2.1;
    m["As"] = 2.05;
    m["Se"] = 1.9;
    m["Br"] = 1.83;
    m["Rb"] = 2.9;
    m["Sr"] = 2.55;
    m["Y"] = 2.4;
    m["Zr"] = 2.3;
    m["Nb"] = 2.15;
    m["Mo"] = 2.1;
    m["Tc"] = 2.05;
    m["Ru"] = 2.05;
    m["Rh"] = 2.0;
    m["Pd"] = 2.05;
    m["Ag"] = 2.1;
    m["Cd"] = 2.2;
    m["In"] = 2.2;
    m["I"] = 1.98;
    m["W"] = 2.1;
    m["Os"] = 2.1;
    m["Ir"] = 2.1;
    return m;
}
static map<string, double> BONDI_RADII = create_bondi_radii();


// #############################################################################
// # Functions
// #############################################################################
//
vector<Vec3> vdw_surface(vector<Vec3> coordinates, vector<string> elements,
                            double scale_factor, double density) {
    // Compute points on the VDW surface of a molecule
    //
    // Parameters
    // ----------
    // coordinates : np.ndarray, shape=(n_atoms, 3)
    //     The cartesian coordinates of the nuclei, in units of ANGSTROMS
    // elements : list, shape=(n_atoms)
    //     The element symbols (C, H, etc) for each of the nuceli
    // scale_factor : float
    //     The points on the molecular surface are set at a distance of
    //     scale_factor * vdw_radius away from each of the atoms.
    // density : float
    //     The (approximate) number of points to generate per square angstrom
    //     of surface area. 1.0 is the default recommended by Kollman & Singh.
    if (coordinates.size() != elements.size()) {
        fprintf(stderr, "coordinate.size doesnot match elements.size");
        exit(1);
    }

    vector<Vec3> surfacepoints;
    vector<double> radii;

    for (size_t i = 0; i < elements.size(); i++) {
        // todo: check for error if element not in BONDI_RADII table
        if (BONDI_RADII.find(elements[i]) != BONDI_RADII.end()) {
            radii.push_back(BONDI_RADII[elements[i]] * scale_factor);
        } else {
            fprintf(stderr, "%s is not a supported element", elements[i].c_str());
            exit(1);
        }
    }

    for (size_t i = 0; i < coordinates.size(); i++) {
        // this could be optimized -- we don't need to compute the dotsphere.
        // at maximum we only need to compute it once per unique element / radius
        vector<Vec3> dots = dotsphere(density * ((4.0/3.0) * M_PI * radii[i]*radii[i]*radii[i]));

        for (size_t j = 0; j < dots.size(); j++)
            dots[j] = coordinates[i] + dots[j]*radii[i];

        // all of the atoms that i is close to
        typedef std::pair<double, size_t> Neighbor;
        vector<Neighbor> neighbors;
        for (size_t j = 0; j < coordinates.size(); j++) {
            if (i == j)
                continue;
            double d = (coordinates[i] - coordinates[j]).norm();
            if (d < (radii[i] + radii[j])) {
                neighbors.push_back(make_pair(d, j));
            }
        }
        sort(neighbors.begin(), neighbors.end());

        for (size_t k = 0; k < dots.size(); k++) {
            int accessible = 1;
            for (vector<Neighbor>::iterator it = neighbors.begin() ; it != neighbors.end(); ++it) {
                if ((coordinates[(*it).second] - dots[k]).norm() < radii[(*it).second]) {
                    accessible = 0;
                    break;
                }
            }
            if (accessible)
                surfacepoints.push_back(dots[k]);
        }
    }

    return surfacepoints;
}