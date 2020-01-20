#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "constants.h"
#include "graphene_initialization.h"


using namespace std;


void graphene_initialization(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long max_atoms, long bins_x, long bins_y, int distinct_layers, int procedure, bool construct)
{
    //procedure = 1 - for initialization of nucleation_type
    //procedure = 2 - for initialization of configuration_type and bin_type
    if (procedure == 1)
    {
        if (construct == true)
        {
            for (long i=0; i<N; i++)
            {
                nucleation[i].xy = new double [Dim];
                nucleation[i].xy[0] = -1;//not really required
                nucleation[i].xy[1] = -1;//not really required
                nucleation[i].angle = -1;
                nucleation[i].neighbors = new long *[quadrants];

                for (int j=0; j<quadrants; j++)
                {
                    nucleation[i].neighbors[j] = new long [nucl_neigh_count];
                    for (int k=0; k<nucl_neigh_count; k++)
                    {
                    nucleation[i].neighbors[j][k] = -1;
                    }
                }
            }
        }
        if (construct == false)
        {
            for (long i=0; i<N; i++)
            {
                for (int j=0; j<quadrants; j++)
                {
                    delete [] nucleation[i].neighbors[j];
                }
                delete [] nucleation[i].xy;
                delete [] nucleation[i].neighbors;
            }
            delete [] nucleation;
        }
    }

    if (procedure == 2)
    {
        if (construct == true)
        {
            //initializing bin
            for (int i=0; i<distinct_layers; i++)//each layer has its own bin
            {
                bin[i].next = new long [max_atoms];
                bin[i].tot_bins = new long [Dim];
                bin[i].tot_bins[0] = bins_x;
                bin[i].tot_bins[1] = bins_y;
                bin[i].bin_length = new double [Dim];
                bin[i].total_atoms = 0;

                for (int j=0; j<max_atoms; j++)
                {
                    bin[i].next[j] = -1;
                }

                for (int j=0; j<Dim; j++)
                {
                    bin[i].bin_length[j] = L[j]/bin->tot_bins[j];
                }

                bin[i].bin_list = new long *[bins_x];

                for (long j=0; j<bins_x; j++)
                {
                    bin[i].bin_list[j] = new long [bins_y];
                    for (long k=0; k<bins_y; k++)
                    {
                        bin[i].bin_list[j][k] = -1;
                    }
                }
            }


            //initializing configuration
            for (long i=0; i<max_atoms*distinct_layers; i++)//there is only configuration construct, i.e. the next layers data is input after multiples of max_atoms
            {
                configuration[i].r = new double [Dim];
                configuration[i].bonded_atom = new long [3];
                configuration[i].bonded_angle = new double [3];
                configuration[i].neigh_count = 0;
                configuration[i].saturated = false;
                configuration[i].grain = -1;//check whether -1 is required
                configuration[i].mirror_no = -1;//indicates whether it has a mirror in the higher layer, changed to true when its mirror in the next atom is accepted
                for (int j=0; j<3; j++)
                {
                    configuration[i].bonded_atom[j] = -1;
                    configuration[i].bonded_angle[j] = -1;//check whether -1 is required
                }
            }
        }

        if (construct == false)
        {

            for (long i=0; i<max_atoms*distinct_layers; i++)
            {
                delete [] configuration[i].r;
                delete [] configuration[i].bonded_atom;
                delete [] configuration[i].bonded_angle;
            }
            delete [] configuration;

            for (int i=0; i<distinct_layers; i++)
            {
                for (long j=0; j<bins_x; j++)
                {
                    delete [] bin[i].bin_list[j];
                }
                delete [] bin[i].bin_list;
                delete [] bin[i].tot_bins;
                delete [] bin[i].next;
                delete [] bin[i].bin_length;
            }
            delete [] bin;

        }
    }

}
