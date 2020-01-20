#ifndef _initialization
#define _initialization

class nucleation_type
{
public:
	double *xy;
	double angle;
	long **neighbors;
};

class configuration_type
{
public:
	double *r;
	long *bonded_atom;
	double *bonded_angle;
	int neigh_count;
	bool saturated;
	long grain;
	long mirror_no;
};

class bin_type
{
public:
	long **bin_list;
	long *next;
	long *tot_bins;
	double *bin_length;
	long total_atoms;
};

void graphene_initialization(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long max_atoms, long bins_x, long bins_y, int layers, int procedure, bool construct);

#endif
