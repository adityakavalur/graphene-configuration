#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "constants.h"
#include "graphene_initialization.h"
#include "graphene_grain_area_calc.h"

using namespace std;

void graphene_grain_area_calc(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long max_atoms, int distinct_layers, const char *grainarea_file)
{

    double *grain_area;
    grain_area = new double [N];

    for (long i=0; i<N; i++)
    {
        grain_area[i] = 0;
    }

    for (int i=0; i<distinct_layers; i++)
    {
        for (long j=0; j<bin[i].tot_bins[0]; j++)
        {
            for (long k=0; k<bin[i].tot_bins[1]; k++)
            {
                long *temp_bin;
                temp_bin = new long [max_atoms];//maximum possible atoms in a bin, if there is only one bin

                for (long l=0; l<max_atoms; l++)
                {
                    temp_bin[l] = -1;
                }

                long current_atom, current_atom_reduced, current_grain, atom_counter=0;//atom counter of the atoms in the current bin
                current_atom = bin[i].bin_list[j][k];

                /* Important */
                current_atom_reduced = current_atom%max_atoms;//used to access memory in bin array

                if (current_atom == -1) {cout << "empty bin, increase bin size for better grain area calculation " << endl; exit(1);}
                temp_bin[atom_counter] = configuration[current_atom].grain;//saving grain number
                while (bin[i].next[current_atom_reduced] != -1)
                {
                    atom_counter++;
                    current_atom = bin[i].next[current_atom_reduced];
                    current_atom_reduced = current_atom%max_atoms;
                    temp_bin[atom_counter] = configuration[current_atom].grain;
                }
                atom_counter++;//one pending for real value

                for (long l=0; l<atom_counter; l++)
                {
                    current_grain = temp_bin[l];
                    grain_area[current_grain] = grain_area[current_grain] + (bin[i].bin_length[0]*bin[i].bin_length[1]/atom_counter);//all atoms in the bin are given equal weightage for area
                }

                delete [] temp_bin;
            }

        }

    }



    //making a histogram of the grain area distribution
    double max_area, min_area;
    int bin_no = 25;//number of bins
    max_area = 0;
    min_area = L[0]*L[1];
    for (long i=0; i<N; i++)
    {
        if (grain_area[i] > max_area) {max_area = grain_area[i];}
        if (grain_area[i] < min_area) {min_area = grain_area[i];}
    }

    long *histogram_area;
    double range = max_area-min_area;
    double bin_width = range/bin_no;

    histogram_area = new long [bin_no];
    for (long i=0; i<bin_no; i++)
    {
        histogram_area[i] = 0;
    }



    long temp;
    for (long i=0; i<N; i++)
    {
        temp = floor(grain_area[i]-min_area)/bin_width;
        if(temp==bin_no) {temp = bin_no;}
        histogram_area[temp] ++;
    }

    //dumping grain area distribution
    ofstream myfile;
    myfile.open (grainarea_file);
    for (long i=0; i<bin_no; i++)
    {
        myfile << min_area+(i+0.5)*bin_width << " " << histogram_area[i] << endl;
    }



    delete [] grain_area;
    delete [] histogram_area;
}
