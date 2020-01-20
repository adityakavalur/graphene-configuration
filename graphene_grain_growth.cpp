#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "constants.h"
#include "graphene_initialization.h"
#include "graphene_nucleation.h"
#include "graphene_atom_generation.h"
#include "graphene_grain_growth.h"

using namespace std;

void graphene_dump_configuration(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long max_atoms, int layer_number, long atom_count, long time_step, int distinct_layers, const char *output_file)
{

    ofstream myfile;
    myfile.open (output_file, ios::app);

    long total_atoms = atom_count;//system total count, initialized with the atom count of the layer under construction not yet updated in the bin
    for (int i=0; i<distinct_layers; i++)
    {
        total_atoms = total_atoms + bin[i].total_atoms;
    }

    myfile << total_atoms << endl;
    myfile << "Atoms. Timestep: " << time_step << endl;
    for (int i=0; i<layer_number; i++)
    {
        for (long j=0; j<bin[i].total_atoms; j++) //printing atoms from already completed layers
        {
            long current_atom = i*max_atoms+j;
            myfile << "1 " << configuration[current_atom].r[0] << " " << configuration[current_atom].r[1] << " " << ((i+0.5)*inter_layer_distance) << endl;
        }
    }
        long j=0;//for below loop
        while (j < atom_count)//printing the layer currently under growth, either the nucleation sites or add_layer
        {
            myfile << "1 " << configuration[j].r[0] << " " << configuration[j].r[1] << " " << ((layer_number-0.5)*inter_layer_distance) << endl;
            j++;
        }

    myfile.close();

}

long graphene_grain_growth(nucleation_type *nucleation, configuration_type *configuration, bin_type *bin, long N, double *L, long max_atoms, int layer_number, long atom_count, long time_step, int distinct_layers, const char *output_file)
{

    cout <<"LAYER NUMBER="<<layer_number<<endl;


    bool system_saturation = false; // whether more atoms can be added to the system, re-initialized for current layer


    long *no_atom_current_step ;//atoms generated in the present time step
    no_atom_current_step = new long;

    if (time_step == 0)
    {

        int *temp_bin;
        temp_bin = new int [Dim];



        //time-step zero needs to transfer the nucleation site information to the configuration class
        for (long i=0; i<N; i++)
        {
            configuration[i].grain = i;
            for (int j=0; j<Dim; j++)
            {
                configuration[i].r[j] = nucleation[i].xy[j];
            }
        }



        //time-step zero binning of current atoms only for the first layer since all nucleation points are in the first layer
        for (long i=0; i<N; i++)
        {
            for (int j=0; j<Dim; j++)
            {
                temp_bin[j] = floor(configuration[i].r[j]/bin[layer_number-1].bin_length[j]);

                if (temp_bin[j]<0 || temp_bin[j]>=bin[layer_number-1].tot_bins[j]) {
                   cout <<"Bin numbering error!"<<endl;
                   exit(1);
                }
                //if (temp_bin[j]==bin[layer_number-1].tot_bins[j]) {temp_bin[j] -- ;}
            }
            long bin_atom = bin[layer_number-1].bin_list[temp_bin[0]][temp_bin[1]];
            if (bin_atom == -1) {bin[layer_number-1].bin_list[temp_bin[0]][temp_bin[1]] = i;}
            else
            {
                while (bin[layer_number-1].next[bin_atom] != -1)
                {
                    bin_atom = bin[layer_number-1].next[bin_atom];
                }
                bin[layer_number-1].next[bin_atom] = i;
            }
        }

        delete [] temp_bin;
    }


    //dumping configuration of an initial time-step when the function is called either the nucleation output or graphene_add_layer
    graphene_dump_configuration(nucleation, configuration, bin, N, L, max_atoms, layer_number, atom_count, time_step, distinct_layers, output_file);



    //evolving time steps
    //bool dummy;//to allow the graphene_atom_generation to return as soon as an atom is added
    while (system_saturation == false)
    {
        (*no_atom_current_step)=0;//reset for current time step

        //cout <<(max_atoms*(layer_number-1))<<" "<<(max_atoms*(layer_number-1)+atom_count)-1<<endl;

        for (long i=(max_atoms*(layer_number-1)); i<(max_atoms*(layer_number-1)+atom_count); i++)
        {
            if (configuration[i].saturated == false)
            {
                graphene_atom_generation(nucleation, configuration, bin, N, L, atom_count, i, no_atom_current_step, layer_number, max_atoms);
            }
        }
        atom_count += (*no_atom_current_step);
        if (layer_number == 2 || layer_number == 3) {cout << "one time step of grain growth in layer - " << layer_number << endl;}
        if (*no_atom_current_step == 0) {system_saturation = true;}
        if (atom_count > max_atoms) {cout << "atoms in a layer greater than max_atoms, recalculate max_atoms" << endl; exit(1);}
        time_step ++;

        //dump files
        graphene_dump_configuration(nucleation, configuration, bin, N, L, max_atoms, layer_number, atom_count, time_step, distinct_layers, output_file);

        cout <<time_step<<" "<<atom_count<<" "<<system_saturation<<endl;

    }

    for (long i=(max_atoms*(layer_number-1)); i<(max_atoms*(layer_number-1)+atom_count); i++) {
       //cout<<i<<" "<<configuration[i].saturated<<endl;
       if (!configuration[i].saturated) {
          cout <<"Pre-mature termination!"<<endl;
          exit(1);
       }
    }


    cout << "atom_count of current layer - " << layer_number << " : " << atom_count << "; single_grain/max_count allowable of a layer -" << max_atoms << endl;
    bin[layer_number-1].total_atoms = atom_count;//grain growth on current layer done, updating the total atom count to the bin of the corresponding layer

    delete [] no_atom_current_step;
    return(time_step);


}
