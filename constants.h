#ifndef _constants
#define _constants

const int Dim = 2; // dimension of construction, since layers are added one by one
const double bond_length = 1.396469226;// bond length 1.42A
const double inter_layer_distance=3.35; // inter-layer distance 3.35A
const double atom_min_dist2 = 1;//square of minimum distance acceptable between two atoms (actual value = 1A)
const double nucl_min_dist2 = 25;//square of minimum distance acceptable between two nucleation points (actual value = 5A)
const int max_neighbors = 3;//maximum number of neighbors inclusive of bonded and non-bonded atoms within a radius of (lattice const - 1) or (1.732050808*bond_length - 1)
const int nucl_neigh_count = 50; //maximum number of neighboring nucleation points considered while generating angles for the grain, two in each of the four quadrants
const int quadrants = 4; //considering 4 quadrant system
const double min_angle = 0.0; //minimum angle between adjacent grains


#endif
