#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <time.h>
#include "constants.h"
#include "graphene_initialization.h"
#include "graphene_grain_growth.h"
#include "graphene_atom_generation.h"
#define PI 3.14159265
using namespace std;


//used to check whether the atom generated in the function below is acceptable or not
bool graphene_atom_check(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long atom_count, long current_atom, long *no_atom_current_step, double *temp_xy, double angle_generated, int layer_number, long max_atoms)
{

    double *dist_xy;// distance values
    long *temp_bin;
    long *actual_bin;
    int temp_neigh_count = 0;
    dist_xy = new double [Dim];
    temp_bin = new long [Dim];
    actual_bin = new long [Dim];


    //bringing the lattice point back to the original cell
    for (int i=0; i<Dim; i++)
    {
        while (temp_xy[i] > L[i]) {temp_xy[i] -= L[i];}
        while (temp_xy[i] < 0) {temp_xy[i] += L[i];}
    }

    //calculating the bin of the point under consideration
    for (int i=0; i<Dim; i++)
    {
        temp_bin[i] = floor(temp_xy[i]/bin[layer_number-1].bin_length[i]);

        if (temp_bin[i]<0 || temp_bin[i]>=bin[layer_number-1].tot_bins[i]) {
           cout <<"Bin numbering error in graphene_atom_check!"<<endl;
           exit(1);
        }
        //if (temp_bin[i]==bin[layer_number-1].tot_bins[i]) {temp_bin[i] --;}
    }

    //checking whether another atom already exists within 1Angstrom of point under-consideration, need to check only 9 bins since bin size is greater than 2Angstroms
    // also if addition of generated atom comes within lattice constant - 1 of already saturated atom

    for (long i=temp_bin[0]-1; i<temp_bin[0]+2; i++)
    {
        actual_bin[0] = i%bin[layer_number-1].tot_bins[0];
        if (actual_bin[0] < 0) {actual_bin[0] += bin[layer_number-1].tot_bins[0];}

        if (actual_bin[0]<0 || actual_bin[0]>=bin[layer_number-1].tot_bins[0]) {
           cout <<"Bin indexing error in x!"<<endl;
           exit(1);
        }

        for (long j=temp_bin[1]-1; j<temp_bin[1]+2; j++)
        {
            actual_bin[1] = j%bin[layer_number-1].tot_bins[1];
            if (actual_bin[1] < 0) {actual_bin[1] += bin[layer_number-1].tot_bins[1];}

            if (actual_bin[1]<0 || actual_bin[1]>=bin[layer_number-1].tot_bins[1]) {
               cout <<"Bin indexing error in y!"<<endl;
               exit(1);
            }


            long atom_s = bin[layer_number-1].bin_list[actual_bin[0]][actual_bin[1]];

            while (atom_s != -1)
            {

                /* Important */
                long atom_s_reduced = atom_s%max_atoms;//we will use atom_s_reduced to access arrays of bin.next as they reduced to the equivalent value below max_atoms, this is for atoms beyond base layer only

                for (int k=0; k<2; k++)
                {
                    dist_xy[k] = abs(temp_xy[k] - configuration[atom_s].r[k]);
                    if (dist_xy[k] > L[k]/2) {dist_xy[k] = L[k]-dist_xy[k];}
                }

                //(1) atom not accepted: there is an atom within 1 angstrom.
                if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < atom_min_dist2) {
                   delete [] actual_bin; delete [] temp_bin; delete [] dist_xy; return false;//atom not accepted
                }
                if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < pow(((bond_length*sqrt(3.0))-1),2))
                {
                    temp_neigh_count ++;

                    //(2) atom not accepted: there is a neighbor that has already its own 3 neighbors.
                    //if (configuration[atom_s].saturated == true){delete [] actual_bin; delete [] temp_bin; delete [] dist_xy; return false;} //atom not accepted
                    if (configuration[atom_s].neigh_count == max_neighbors) {
                       delete [] actual_bin; delete [] temp_bin; delete [] dist_xy; return false;
                    }
                }
                //(3) atom not accepted: this attempted atom has more than 3 neighbors.
                if (temp_neigh_count > max_neighbors) {delete [] actual_bin; delete [] temp_bin; delete [] dist_xy; return false;}

                atom_s = bin[layer_number-1].next[atom_s_reduced];//using atom_s_reduced to access the next atom in the bin
            }
        }
    }

    //atom accepted
    //updating count
    (*no_atom_current_step) ++;

    //updating bin
    long generated_atom = (atom_count - 1  + *no_atom_current_step)+ (max_atoms*(layer_number-1));
    int atom_check = bin[layer_number-1].bin_list[temp_bin[0]][temp_bin[1]];

    /* Important */
    int atom_check_reduced=atom_check%max_atoms;//we will use atom_check reduced to access memory in the array as it has been reduced to an equivalent value below max_atoms


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

    //updating configuration
    for (int i=0; i<Dim; i++)
    {
        configuration[generated_atom].r[i] = temp_xy[i];
    }
    configuration[generated_atom].bonded_atom[0] = current_atom;
    configuration[generated_atom].neigh_count = temp_neigh_count;
    if (temp_neigh_count == max_neighbors) {configuration[generated_atom].saturated = true;}

    double angle_secondary = angle_generated-180+360;//angle of generated atom with the current atom with the horizontal direction (i.e. zig-zag direction)
    while (angle_secondary > 360) {angle_secondary = angle_secondary-360;}
    while (angle_secondary < 0 ) {angle_secondary = angle_secondary+360;}
    configuration[generated_atom].bonded_angle[0] = angle_secondary;

    int i=0;
    while (configuration[current_atom].bonded_atom[i] != -1)
    {
        i++;
    }
    configuration[current_atom].bonded_atom[i] = generated_atom;
    while (angle_generated > 360) {angle_generated = angle_generated -360;}
    while (angle_generated < 0) {angle_generated = angle_generated + 360;}
    configuration[current_atom].bonded_angle[i] = angle_generated;
    configuration[generated_atom].grain = configuration[current_atom].grain;


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

                    if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < pow(((bond_length*sqrt(3.0))-1),2))
                    {
                        configuration[atom_s].neigh_count ++;
                    }
                    if ( configuration[atom_s].neigh_count == max_neighbors) {configuration[atom_s].saturated = true;}

                }
                atom_s = bin[layer_number-1].next[atom_s_reduced];//using atom_s_reduced to access the next atom in the bin
            }
        }
    }

    delete [] actual_bin;
    delete [] dist_xy;// distance values
    delete [] temp_bin;

    return true;

}

//generating a new atom
void graphene_atom_generation(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long atom_count, long current_atom, long *no_atom_current_step, int layer_number, long max_atoms)
{
    //srand (time(NULL));

    bool atom_added = false;
    int temp_rand;
    int random_value, random_value2, random_value3;

    double *temp_xy;// generated co-ordinates
    temp_xy = new double [Dim];



    if (configuration[current_atom].bonded_atom[0] == -1) //no neighbor, so 3 possible addition sites, only used in timestep=0 and no generated atom can be rejected
    {

        if (configuration[current_atom].bonded_atom[1] != -1 || configuration[current_atom].bonded_atom[2] != -1) {
           cout <<"Error0!!!"<<endl;
           exit(1);
        }

        double angle_generated;

        //first option of neighbor
        temp_rand = rand();
        random_value = temp_rand % 3 ; // generated integer [0,2]

        angle_generated = (random_value*120.0 + 30.0 + nucleation[configuration[current_atom].grain].angle);//add 30deg so that 0 deg is zig-zag direction for a 0deg misorientation grain

        //cout <<current_atom<<": random_value0="<<random_value<<" "<<angle_generated-nucleation[configuration[current_atom].grain].angle<<" "<<nucleation[configuration[current_atom].grain].angle<<" "<<configuration[current_atom].grain<<endl;

        temp_xy[0] = configuration[current_atom].r[0] + bond_length*cos(angle_generated*PI/180.0);
        temp_xy[1] = configuration[current_atom].r[1] + bond_length*sin(angle_generated*PI/180.0);

        atom_added = graphene_atom_check(nucleation, configuration, bin, N, L, atom_count, current_atom, no_atom_current_step, temp_xy, angle_generated, layer_number, max_atoms);

        if (atom_added == true)
        {
            delete [] temp_xy;
            return;
            //return atom_added;//no significance of value returned, it must return however
        }
        //second option of neighbor
        else
        {
            do
            {
                random_value2 = rand() % 3 ; // generated integer [0,2]
            }
            while (random_value2 == random_value);
        }
        angle_generated = (random_value2*120.0 + 30.0 + nucleation[configuration[current_atom].grain].angle);
        temp_xy[0] = configuration[current_atom].r[0] + bond_length*cos(angle_generated*PI/180.0);
        temp_xy[1] = configuration[current_atom].r[1] + bond_length*sin(angle_generated*PI/180.0);
        atom_added = graphene_atom_check(nucleation, configuration, bin, N, L, atom_count, current_atom, no_atom_current_step, temp_xy, angle_generated, layer_number, max_atoms);

        if (atom_added == true)
        {
            delete [] temp_xy;
            return;
            //return atom_added;//no significance of value returned, it must return however
        }
        //third option of neighbor
        else
        {

            do
            {
                random_value3 = rand() % 3 ; // generated integer [0,2]
            }
            while (random_value3 == random_value || random_value3 == random_value2);
        }
        angle_generated = (random_value3*120.0 + 30.0 + nucleation[configuration[current_atom].grain].angle);
        temp_xy[0] = configuration[current_atom].r[0] + bond_length*cos(angle_generated*PI/180.0);
        temp_xy[1] = configuration[current_atom].r[1] + bond_length*sin(angle_generated*PI/180.0);
        atom_added = graphene_atom_check(nucleation, configuration, bin, N, L, atom_count, current_atom, no_atom_current_step, temp_xy, angle_generated, layer_number, max_atoms);

        configuration[current_atom].saturated = true;//regardless of any atom added the current atom is saturated

        delete [] temp_xy;
        return;
        //return atom_added;//no significance of value returned, it must return however


    }



    if (configuration[current_atom].bonded_atom[1] == -1)
    {

        if (configuration[current_atom].bonded_atom[0] == -1 || configuration[current_atom].bonded_atom[2] != -1) {
           cout <<"Error1!!!"<<endl;
           exit(1);
        }

        //retrieving information of the 1st bonded atom
        double angle_generated, angle_retrieved;

        angle_retrieved = configuration[current_atom].bonded_angle[0];

        //first option of neighbor
        random_value = rand() % 2 ; // generated integer [0,1]
        random_value ++;//changing it to [1,2]
        angle_generated = (random_value*120.0 + angle_retrieved);


        temp_xy[0] = configuration[current_atom].r[0] + bond_length*cos(angle_generated*PI/180.0);
        temp_xy[1] = configuration[current_atom].r[1] + bond_length*sin(angle_generated*PI/180.0);

        atom_added = graphene_atom_check(nucleation, configuration, bin, N, L, atom_count, current_atom, no_atom_current_step, temp_xy, angle_generated, layer_number, max_atoms);

        if (atom_added == true)
        {
            delete [] temp_xy;
            return;
            //return atom_added;//no significance of value returned, it must return however
        }
        //second option of neighbor
        else
        {

            do
            {
                random_value2 = rand() % 2 ; // generated integer [0,1]
                random_value2 ++;//changing the range to [1,2]
            }
            while (random_value2 == random_value);
        }

        angle_generated = (random_value2*120 + angle_retrieved);
        temp_xy[0] = configuration[current_atom].r[0] + bond_length*cos(angle_generated*PI/180);
        temp_xy[1] = configuration[current_atom].r[1] + bond_length*sin(angle_generated*PI/180);
        atom_added = graphene_atom_check(nucleation, configuration, bin, N, L, atom_count, current_atom, no_atom_current_step, temp_xy, angle_generated, layer_number, max_atoms);

        configuration[current_atom].saturated = true;//regardless of any atom added the current atom is saturated

        delete [] temp_xy;
        return;
        //return atom_added;//no significance of value returned, it must return however

    }



    if (configuration[current_atom].bonded_atom[2] == -1)
    {

        if (configuration[current_atom].bonded_atom[0] == -1 || configuration[current_atom].bonded_atom[1] == -1) {
           cout <<"Error2!!!"<<endl;
           exit(1);
        }


        double angle_generated, angle_retrieved, angle_retrieved2;

        //retrieving information of the 1st bonded atom
        angle_retrieved = configuration[current_atom].bonded_angle[0];


        //retrieving information of the 2nd bonded atom
        angle_retrieved2 = configuration[current_atom].bonded_angle[1];


        if (abs(angle_retrieved2-angle_retrieved)>180.0)//assuming each angle may waver a little from 120 this case is the one where the angle is obtuse between the two bonds in the anti-clockwise measurement direction
        {
            angle_generated = (angle_retrieved2+angle_retrieved)/2.0;
            temp_xy[0] = configuration[current_atom].r[0] + bond_length*cos((angle_generated)*PI/180.0);
            temp_xy[1] = configuration[current_atom].r[1] + bond_length*sin((angle_generated)*PI/180.0);
            atom_added = graphene_atom_check(nucleation, configuration, bin, N, L, atom_count, current_atom, no_atom_current_step, temp_xy, angle_generated, layer_number, max_atoms);

            configuration[current_atom].saturated = true;//regardless of any atom added the current atom is saturated

            delete [] temp_xy;
            return;
            //return atom_added;//no significance of value returned, it must return however

        }
        if (abs(angle_retrieved-angle_retrieved2)<180.0)//assuming each angle may waver a little from 120 this case is the one where the angle is acute between the two bonds in the anti-clockwise measurement direction
        {

            if (angle_retrieved > angle_retrieved2) {angle_generated = angle_retrieved+120.0;}
            if (angle_retrieved2 > angle_retrieved) {angle_generated = angle_retrieved2+120.0;}

            /*if (angle_generated>360.0) {
               cout <<"angle_generated="<<angle_generated<<endl;
            }*/

            temp_xy[0] = configuration[current_atom].r[0] + bond_length*cos((angle_generated)*PI/180.0);
            temp_xy[1] = configuration[current_atom].r[1] + bond_length*sin((angle_generated)*PI/180.0);
            atom_added = graphene_atom_check(nucleation, configuration, bin, N, L, atom_count, current_atom, no_atom_current_step, temp_xy, angle_generated, layer_number, max_atoms);

            configuration[current_atom].saturated = true;//regardless of any atom added the current atom is saturated

            delete [] temp_xy;
            return;
            //return atom_added;//no significance of value returned, it must return however

        }
    }

    //shouldn't reach here on first layer, but will reach here for subsequent layers as saturation boolean of individual atoms is not updated in the the graphene_add_layer
    delete [] temp_xy;
    configuration[current_atom].saturated = true;//no atom added but atom saturated
    return;
    //return atom_added;//no significance of value returned, it must return however
}
