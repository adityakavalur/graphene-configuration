#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "constants.h"
#include "graphene_initialization.h"
#include "graphene_grain_growth.h"
#include "graphene_add_layer.h"
#define PI 3.14159265
using namespace std;

bool graphene_atom_layer_check(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long max_atoms, int layer_number, long generated_atom, double *temp_xy, long original_atom)
{

    //bringing the lattice point back to the original cell
    for (int i=0; i<Dim; i++)
    {
        while (temp_xy[i] > L[i]) {temp_xy[i] -= L[i];}
        while (temp_xy[i] < 0) {temp_xy[i] += L[i];}
    }


    double *dist_xy;// distance values
    long *temp_bin;
    long *actual_bin;
    dist_xy = new double [Dim];
    temp_bin = new long [Dim];
    actual_bin = new long [Dim];
    long temp_neigh_count = 0;//re-initialized for every new atom

    //calculating the bin of the point under consideration
    for (int i=0; i<Dim; i++)
    {
        temp_bin[i] = floor(temp_xy[i]/bin[layer_number-1].bin_length[i]);
        if (temp_bin[i]<0 || temp_bin[i]>=bin[layer_number-1].tot_bins[i]) {
           cout <<"Bin numbering error in graphene_atom_layer_check!"<<endl;
           exit(1);
        }
        //if (temp_bin[i]==bin[layer_number-1].tot_bins[i]) {temp_bin[i] --;}
    }


    for (long i=temp_bin[0]-1; i<temp_bin[0]+2; i++)
    {
        actual_bin[0] = i%bin[layer_number-1].tot_bins[0];
        if (actual_bin[0] < 0) {actual_bin[0] += bin[layer_number-1].tot_bins[0];}

        for (long j=temp_bin[1]-1; j<temp_bin[1]+2; j++)
        {
            actual_bin[1] = j%bin[layer_number-1].tot_bins[1];
            if (actual_bin[1] < 0) {actual_bin[1] += bin[layer_number-1].tot_bins[1];}
            long atom_s = bin[layer_number-1].bin_list[actual_bin[0]][actual_bin[1]];


            while (atom_s != -1)
            {

                /* Important */
                long atom_s_reduced = atom_s%max_atoms;//we will use atom_s_reduced to access arrays of bin.next as they reduced to the equivalent value below max_atoms


                for (int k=0; k<Dim; k++)
                {
                    dist_xy[k] = abs(temp_xy[k] - configuration[atom_s].r[k]);
                    if (dist_xy[k] > L[k]/2) {dist_xy[k] = L[k]-dist_xy[k];}
                }
                if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < atom_min_dist2) { delete [] temp_bin; delete [] dist_xy; delete [] actual_bin; return false;}//atom not accepted
                if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < pow(((bond_length*1.732050808)-1),2))
                {
                    temp_neigh_count ++;
                    if (configuration[atom_s].neigh_count == max_neighbors) {
                       delete [] actual_bin; delete [] temp_bin; delete [] dist_xy; return false;
                    }
                    //if (configuration[atom_s].saturated == true){delete [] actual_bin; delete [] temp_bin; delete [] dist_xy; return false;} //atom not accepted
                }
                if (temp_neigh_count  > max_neighbors) { delete [] temp_bin; delete [] dist_xy; delete [] actual_bin; return false;}//atom not accepted
                atom_s = bin[layer_number-1].next[atom_s_reduced];//using atom_s_reduced to access the next atom in the bin
            }
        }
    }

    //atom accepted
    //updating configuration

    configuration[original_atom].mirror_no = generated_atom;
    configuration[generated_atom].grain = configuration[original_atom].grain;
    configuration[generated_atom].neigh_count = temp_neigh_count;
    if (temp_neigh_count == max_neighbors) {configuration[generated_atom].saturated = true;}

    for (int i=0; i<Dim; i++)
    {
        configuration[generated_atom].r[i] = temp_xy[i];
    }


    //updating bin

    int atom_check = bin[layer_number-1].bin_list[temp_bin[0]][temp_bin[1]];


    /* Important */
    int atom_check_reduced = atom_check%max_atoms;//we will use atom_check reduced to access memory in the array as it has been reduced to an equivalent value below max_atoms



    if (atom_check == -1)
    {
        bin[layer_number-1].bin_list[temp_bin[0]][temp_bin[1]] = generated_atom;
    }
    else
    {
        while (bin[layer_number-1].next[atom_check_reduced] != -1)
        {
            atom_check = bin[layer_number-1].next[atom_check_reduced];
            atom_check_reduced = atom_check%max_atoms;
        }
        bin[layer_number-1].next[atom_check_reduced] = generated_atom;
    }


     //updating configuration.neigh_count for non-generated atoms
    for (long i=temp_bin[0]-1; i<temp_bin[0]+2; i++)
    {
        actual_bin[0] = i%bin[layer_number-1].tot_bins[0];
        if (actual_bin[0] < 0) {actual_bin[0] += bin[layer_number-1].tot_bins[0];}

        for (long j=temp_bin[1]-1; j<temp_bin[1]+2; j++)
        {
            actual_bin[1] = j%bin[layer_number-1].tot_bins[1];
            if (actual_bin[1] < 0) {actual_bin[1] += bin[layer_number-1].tot_bins[1];}


            long atom_s = bin[layer_number-1].bin_list[actual_bin[0]][actual_bin[1]];

            while (atom_s != -1)
            {
                /* Important */
                long atom_s_reduced = atom_s%max_atoms;//we will use atom_s_reduced to access arrays of bin.next as they reduced to the equivalent value below max_atoms, this is for atoms beyond base layer only

                if (atom_s != generated_atom)
                {

                    for (int k=0; k<2; k++)
                    {
                        dist_xy[k] = abs(temp_xy[k] - configuration[atom_s].r[k]);
                        if (dist_xy[k] > L[k]/2) {dist_xy[k] = L[k]-dist_xy[k];}
                    }

                    if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < pow(((bond_length*1.732050808)-1),2))
                    {
                        configuration[atom_s].neigh_count ++;
                    }
                    if ( configuration[atom_s].neigh_count == max_neighbors) {configuration[atom_s].saturated = true;}

                }
                atom_s = bin[layer_number-1].next[atom_s_reduced];//using atom_s_reduced to access the next atom in the bin
            }
        }
    }





    delete [] temp_bin;
    delete [] dist_xy;
    delete [] actual_bin;
    return true;


}





long graphene_add_layer(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long max_atoms, int layer_number, long time_step, int distinct_layers, const char *output_file)
{

    long atom_count=0;
    //adding atoms to the layer

    for (long i=(max_atoms*(layer_number-2)) ; i<((max_atoms*(layer_number-2))+(bin[layer_number-2].total_atoms)); i++)//driven by the atoms in the below layer, to generate a copy of each of them
    {
        bool atom_check;//all these could have been put outside putting only temp_xy inside the layer_number if condition
        long original_atom = i;
        long generated_atom = atom_count+max_atoms*(layer_number-1);
        long original_grain = configuration[original_atom].grain;
        double grain_angle = nucleation[original_grain].angle;

        double *temp_xy;// generated co-ordinates
        temp_xy = new double [Dim];

        if (layer_number == 2)
        {
            temp_xy[0] = ((1.732050808*bond_length/2)*cos(grain_angle*PI/180)) - ((bond_length/2)*sin(grain_angle*PI/180)) + configuration[original_atom].r[0];
            temp_xy[1] = ((1.732050808*bond_length/2)*sin(grain_angle*PI/180)) + ((bond_length/2)*cos(grain_angle*PI/180)) + configuration[original_atom].r[1];
            //cout << temp_xy[0] << " " << configuration[original_atom].r[0] << endl;
            //cout << temp_xy[1] << " " << configuration[original_atom].r[1] << endl;
            //exit(1);
        }

        if (layer_number == 3)
        {
            temp_xy[0] = - ((1.732050808*bond_length/2)*cos(grain_angle*PI/180)) - ((bond_length/2)*sin(grain_angle*PI/180)) + configuration[original_atom].r[0];
            temp_xy[1] = - ((1.732050808*bond_length/2)*sin(grain_angle*PI/180)) + ((bond_length/2)*cos(grain_angle*PI/180)) + configuration[original_atom].r[1];
        }

        //check whether generated atom is acceptable
        atom_check = graphene_atom_layer_check(nucleation, configuration, bin, N, L, max_atoms, layer_number, generated_atom, temp_xy, original_atom);

        if (atom_check == true)
        {
            atom_count ++;
        }
        delete [] temp_xy;
    }

    cout << "maximum possible mirrors added to layer - " << layer_number << " current layer atom count : "<< atom_count << endl;

    //exit(1);

    //updating bonded atoms' information
    for (long i=(max_atoms*(layer_number-2)) ; i<((max_atoms*(layer_number-2))+(bin[layer_number-2].total_atoms)); i++)//driven by the atoms in the below layer
    {
        long original_atom = i;
        if (configuration[original_atom].mirror_no != -1)
        {
            int bonded_atom_count=0;
            long generated_atom = configuration[original_atom].mirror_no;
            for (int j=0; j<3; j++)
            {

                long bonded_atom = configuration[original_atom].bonded_atom[j];
                if (configuration[bonded_atom].mirror_no != -1 && bonded_atom != -1)
                {
                    configuration[generated_atom].bonded_atom[bonded_atom_count] = configuration[bonded_atom].mirror_no;
                    configuration[generated_atom].bonded_angle[bonded_atom_count] = configuration[original_atom].bonded_angle[j];
                    bonded_atom_count ++;
                }
            }

        }

    }



    time_step ++;//time step added for the creation of the new layer

    //grain-growth wherever possible
    time_step = graphene_grain_growth(nucleation, configuration, bin, N, L, max_atoms, layer_number, atom_count, time_step, distinct_layers, output_file);

    return(time_step);
}
