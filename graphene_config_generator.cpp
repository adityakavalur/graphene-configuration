#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <time.h>
#include "constants.h"
#include "graphene_initialization.h"
#include "graphene_nucleation.h"
#include "graphene_grain_growth.h"
#include "graphene_grain_area_calc.h"
#include "graphene_add_layer.h"
using namespace std;

int main ()
{
    cout << "graphene_config_generator_start" << endl;

    srand (time(NULL));

    //miscellaneous declarations
    double *L; //model cell size
    long N; //number of nucleation sites
    long atom_count=0;//atom_count of the layer under construction
    long time_step=0;//time_step
    int layers, distinct_layers=0; //number of layers and distinct layers, distinct layers set to zero to check if the input file has all the required parameters since distinct layers is the last entry
    int procedure; // procedure number for initialization function to determine what is to be initialized
    bool construct; // to determine whether the initialization function is to construct(true) or de-construct(false)
    long max_atoms; //maximum atoms possible, considering a single grain graphene system
    long bins_x, bins_y; //number of bins in X,Y direction
    double min_z;//height of box base in Z direction
    long start_atom_count;// starting number of atom count + 1 ****, these many atoms will be provided by some other file, in case it needs to be merged with some other lammps file
    int base_layer_number;//number of base layer + 1 ****, these many layers will be provided by some other file

    L = new double [Dim];


    //txt file declarations
    string initialfile;//simulation run parameters initial input file
    string outputfile;//dump file of XYZ format
    string grainareafile;//grain area calculation
    string finalconfigfile;//final configuration file
    string angledistrifile;//angle distribution file
    string orientationfile;//grain angle distribution file


    initialfile="graphene_multi-grain_start_file.txt";//only place this is defined
    const char *initial_file = initialfile.c_str();

    /* Initial file - order of parameters
    L[0] - length in X direction
    L[1] - length in Y direction
    N - number of nucleation sites
    layers - number of layers
    distinct_layers - number of distinct layers
    min_z - height of box base in Z direction, in case it needs to be merged with some other lammps file
    start_atom_count - starting number of atom count + 1 ****, these many atoms will be provided by some other file, in case it needs to be merged with some other lammps file
    base_layer_number - number of base layer + 1 ****, these many layers will be provided by some other file
    */


    outputfile = "dump_graphene_multi-grain.xyz";//only place this is defined
    const char *output_file = outputfile.c_str();


    grainareafile = "dump_graphene_multi-grain-area.txt";//only place this is defined
    const char *grainarea_file = grainareafile.c_str();


    finalconfigfile = "dump_graphene_multi-grain_final_config.txt";//only place this is defined
    const char *finalconfig_file = finalconfigfile.c_str();


    angledistrifile = "dump_graphene_multi-grain_angle_distributions.txt";//only place this is defined
    const char *angledistri_file = angledistrifile.c_str();

    orientationfile = "dump_graphene_multi-grain_orientation.txt";//only place this is defined
    const char *orientation_file = orientationfile.c_str();



    //reading input file
    ifstream File;
    File.open(initial_file);
    if (File==NULL)
    {
        cout << "Check input file - " << initialfile << ", file not found"<< endl; abort();
    }

    File >> L[0];
    File >> L[1];
    File >> N;
    File >> layers;
    File >> distinct_layers;
    File >> min_z;
    File >> start_atom_count;
    File >> base_layer_number;
    File.close();

    if (distinct_layers == 0) {cout << "check input file - " << initialfile << ", all required parameters are not listed" << endl; abort();}

    if (distinct_layers > 3) {cout << "check input file - " << initialfile << ", only 3 distinct layers possible in graphene" << endl; abort();}

    if (distinct_layers > layers) {cout << "number of distinct layers is more than number of layers, check input file" << endl; abort();}

    cout << "graphene_config_generator_file_read" << endl;



    //deleting the file from previous run if still present, as all other places have appending flag on
    ofstream myfile;
    myfile.open (output_file);
    myfile.close();



    //declarations in graphene_initialization.h
    nucleation_type *nucleation;
    configuration_type *configuration;
    bin_type *bin;


    //calculating the theoretical maximum number of atoms in given model cell size i.e. assuming single grain for one layer
    //from simplest model of graphene (4 atoms) we get a value slightly lower than 0.4 atoms/A^2 (using a lower value of 1.39A bond length for highest packing)
    //doesn't increment max_atoms dynamically, if atoms exceed this value in one layer the code will exit
    max_atoms = ceil(L[0]*L[1]*0.4);//for one layer
    max_atoms = 1.25*max_atoms; //extra space, since small model cell sizes have a higer ratio of half unit cells to complete ones, which leads to higher number of max atoms for small model cell


    //calculating the number of bins required in X,Y directions, assuming 1.42A bond length taking a starting point of 4A (a√3 - max distance between two graphene atoms + 1A min distance between two atoms) as the bin length in both directions
    int bin_length = ceil((1.732050808*bond_length) + sqrt(atom_min_dist2));// higher value results in more computations while distance calculations, but lower value leads to empty bins which is difficult to treat for grain area distribution
    bins_x = floor(L[0]/bin_length);
    bins_y = floor(L[1]/bin_length);
    //cout << bins_x << " " << bins_y << endl;



    //initializing the nucleation_type class
    nucleation = new nucleation_type[N];
    procedure = 1; //for initialization of nulceation_type
    construct = true;
    graphene_initialization(nucleation, NULL, NULL, N, L, 0, 0, 0, 0, procedure, construct);

    cout << "graphene_config_generator_initialized_nucleation" << endl;



    //carrying out nucleation
    graphene_nucleation(nucleation, N, L, angledistri_file);
    atom_count = N;//updating the atom count

    cout << "graphene_config_generator_after_nucleation" << endl;



    //initialising the configuration_type and bin_type classes
    configuration = new configuration_type [max_atoms*distinct_layers];//there is only one configuration variable, i.e. the next layers data is input after multiples of max_atoms
    bin = new bin_type [distinct_layers];//each distinct pattern layer has its own bin i.e. A, B, C (max of 3)
    procedure = 2;
    construct = true;
    graphene_initialization(nucleation, configuration, bin, N, L, max_atoms, bins_x, bins_y, distinct_layers, procedure, construct);

    cout << "graphene_config_generator_initialized_configuration&bin" << endl;



    //carrying out grain-growth
    int layer_number = 1;//grain growth is to be carried out only for the first layer from the nucleation sites

    time_step = graphene_grain_growth(nucleation,configuration,bin, N, L, max_atoms, layer_number, atom_count, time_step, distinct_layers, output_file);

    cout << "graphene_config_generator_grain_growth " <<bin[layer_number-1].total_atoms<< endl;

    for (long i=0;i!=bin[layer_number-1].total_atoms;++i) {
       if (configuration[i].neigh_count>max_neighbors) {
          cout <<"Too many neighbors: "<<layer_number<<" "<<i<<" "<<configuration[i].neigh_count<<endl;
          exit(1);
       }
       if (configuration[i].neigh_count==0) {
          cout <<"No neighbor: "<<layer_number<<" "<<i<<" "<<configuration[i].neigh_count<<endl;
          exit(1);
       }
       if (configuration[i].neigh_count<max_neighbors) {
          cout <<"WARNING! Too few neighbors: "<<layer_number<<" "<<i<<" "<<configuration[i].neigh_count<<endl;
          cout << i << " " << configuration[i].bonded_atom[0] << " " << configuration[i].bonded_atom[1] << " " << configuration[i].bonded_atom[2] << " " << endl;
       }
    }

    //exit(1);

    //adding layers if required

    for (int i=2; i<(distinct_layers+1); i++)// i must take layer number so it is initiated with 2 and run till distinct_layer
    {
        cout << "layer number initiated - " << i << endl;
        //shouldn't enter this loop if distinct_layers=1 i.e. only one layer pattern is to be built
        time_step = graphene_add_layer(nucleation, configuration, bin, N, L, max_atoms, i, time_step, distinct_layers, output_file);

       for (long j=max_atoms*(i-1);j!=max_atoms*(i-1)+bin[i-1].total_atoms;++j) {
          if (configuration[j].neigh_count>max_neighbors) {
             cout <<"Too many neighbors: "<<layer_number<<" "<<j<<" "<<configuration[j].neigh_count<<endl;
             exit(1);
          }
          if (configuration[j].neigh_count==0) {
             cout <<"No neighbor: "<<layer_number<<" "<<j<<" "<<configuration[j].neigh_count<<endl;
             exit(1);
          }
          if (configuration[j].neigh_count==1) {
             cout <<"WARNING! Too few neighbors: "<<layer_number<<" "<<j<<" "<<configuration[j].neigh_count<<endl;
          }
       }

    }

    //exit(1);

    //dumping the final configuration
    long final_atom_count=0;
    for (int i=0; i<layers; i++)
    {
        int layer_pattern = i%distinct_layers;
        final_atom_count = final_atom_count + (bin[layer_pattern].total_atoms);
    }

    //ofstream myfile;
    myfile.open (finalconfig_file);//to create a new file as this is the first output and delete a file if present from the previous run
    myfile << "LAMMPS input file" << endl;
    myfile << "#" << endl;
    myfile << final_atom_count << " atoms"<< endl;
    myfile << "0 bonds"<< endl;
    myfile << "0 angles"<< endl;
    myfile << "0 dihedrals"<< endl;
    myfile << "0 impropers"<< endl;
    myfile << "#" << endl;
    myfile << layers << " atom types"<< endl;
    myfile << "#" << endl;
    myfile << std::setprecision(15)<< "0 " << L[0] << " xlo xhi" << endl;
    myfile << "0 " << L[1] << " ylo yhi" << endl;
    myfile << "0 " << layers*inter_layer_distance << " zlo zhi" << endl;
    myfile << "#" << endl;
    myfile << "Masses" << endl;
    myfile << "#" << endl;
    myfile << "1 12" << endl;
    myfile << "#" << endl;
    myfile << "Atoms" << endl;
    myfile << "#" << endl;

    final_atom_count = 1+start_atom_count;//re-initialized
    for (long i=0; i<layers; i++)
    {
        int layer_pattern = i%distinct_layers;
        for (long j=0; j<bin[layer_pattern].total_atoms; j++)
        {
            long current_atom = j + (max_atoms*i);
            myfile << std::setprecision(15) << final_atom_count << " 1 " << i+1+base_layer_number << " " << configuration[current_atom].r[0] << " " << configuration[current_atom].r[1] << " " << ((i+0.5)*inter_layer_distance)+min_z  << endl;
            final_atom_count++;
            // 1  is molecule tag
        }
    }
    myfile.close();

    cout << "final_configuration_output_completed" << endl;



    //dumping nucleation and its bonded atoms configuration for comparison with lammps output
    myfile.open (orientation_file);
    myfile << 4*N << " atoms" << endl;
    for (int i=0; i<N; i++)
    {
        myfile << i << " " << configuration[i].r[0] << " " << configuration[i].r[1] << " " << 0.5*inter_layer_distance << endl;
        for (int j=0; j<3; j++)
        {
            long atom_p = configuration[i].bonded_atom[j];
            myfile << atom_p << " " << configuration[atom_p].r[0] << " " << configuration[atom_p].r[1] << " " << 0.5*inter_layer_distance << " " << configuration[i].bonded_angle[j] << endl;
        }
    }
    myfile.close();


    //calculation grain area distribution
    if (N > 1)
    {
        graphene_grain_area_calc(nucleation, configuration, bin, N, L, max_atoms, distinct_layers, grainarea_file);
        cout << "graphene_area_calculation_completed" << endl;
    }

    //checking which grains border each other and if their angle of mis-match meet the minimum criteria set in graphene_nucleation.cpp
    //only for first layer
    long counter_5=0;
    long counter_10 =0;
    long counter_30 = 0;
    long counter =0;
    for (long i=0; i<bin[0].tot_bins[0]; i++)
    {
        for (long j=0; j<bin[0].tot_bins[1]; j++)
        {
            long atom_p = bin[0].bin_list[i][j];
            while (atom_p != -1)
            {
                long *actual_bin;
                actual_bin = new long [Dim];
                for (long k=i-1; k<i+2; k++)
                {
                    actual_bin[0] = k%(bin[0].tot_bins[0]);
                    if (actual_bin[0]<0) {actual_bin[0] += bin[0].tot_bins[0];}
                    {
                        for (long l=j-1; l<j+2; l++)
                        {
                            actual_bin[1] = l%(bin[0].tot_bins[1]);
                            if (actual_bin[1]<0) {actual_bin[1] += bin[0].tot_bins[1];}
                            {
                                long atom_s = bin[0].bin_list[actual_bin[0]][actual_bin[1]];
                                while (atom_s != -1)
                                {
                                    double *temp_dist;
                                    temp_dist = new double [Dim];

                                    for (int m=0; m<Dim; m++)
                                    {
                                        temp_dist[m] = abs (configuration[atom_p].r[m] - configuration[atom_s].r[m]);
                                    }
                                    if ( pow (temp_dist[0],2.0)+ pow(temp_dist[1],2.0) < pow((sqrt(3.0)*bond_length -1.0),2.0))
                                    {
                                        if (configuration[atom_p].grain != configuration[atom_s].grain) //i.e. they don't belong to the same grain
                                        {
                                            double angle_diff = abs (nucleation[configuration[atom_p].grain].angle - nucleation[configuration[atom_s].grain].angle);
                                            if (angle_diff >= 60.0) {angle_diff = angle_diff - 60.0;}
                                            if (angle_diff >= 30.0) {angle_diff = 60.0 - angle_diff;}

                                            if (min_angle > angle_diff) {counter ++; cout << "grains less than border each other, with angle difference of " << angle_diff << endl;}
                                            if (angle_diff < 5.0) {counter_5++;}
                                            if (angle_diff < 10.0) {counter_10++;}
                                            if (angle_diff < 30.0) {counter_30++;}

                                        }
                                    }

                                    atom_s = bin[0].next[atom_s];
                                }
                            }
                        }
                    }

                }

                delete [] actual_bin;
                atom_p = bin[0].next[atom_p];
            }
        }
    }

    cout << "Number of atoms at grain boundary wrt min angle set " << min_angle << " deg : " << counter/2 << endl;
    cout << "Number of atoms at grain boundary with angles less than 5 deg : " << counter_5/2 << endl;
    cout << "Number of atoms at grain boundary with angles less than 10 deg :  " << counter_10/2 << endl;
    cout << "Number of atoms at grain boundary : " << counter_30/2 << endl;

    //checking the two conditions in atom generation: minimum distance and maximum no of neighbors
    for (int i=0; i<distinct_layers; i++)
    {
        for (long j=(max_atoms*i); j<((max_atoms*i)+bin[i].total_atoms); j++)
        {
            long atom_p = j;
            int neighbor_count = 0; //re-initialized for ever primary atom
            long *neighbor_list;
            neighbor_list = new long [4];
            long *atom_p_bin;//bin number of primary atom;
            atom_p_bin = new long [Dim];

            double *dist_xy;
            dist_xy = new double [Dim];

            for (int k=0; k<Dim; k++)
            {
                atom_p_bin[k] = floor(configuration[atom_p].r[k]/bin[i].bin_length[k]);
                if (atom_p_bin[k] == bin[i].tot_bins[k]) { atom_p_bin[k] -- ;}
            }

            long *atom_s_bin;
            atom_s_bin = new long [Dim];

            for (long k=atom_p_bin[0]-1; k<atom_p_bin[0]+2; k ++)
            {
                atom_s_bin[0] = k%bin[i].tot_bins[0];
                if (atom_s_bin[0] < 0) {atom_s_bin[0] += bin[i].tot_bins[0];}
                for (long l=atom_p_bin[1]-1; l<atom_p_bin[1]+2; l++)
                {
                    atom_s_bin[1] = l%bin[i].tot_bins[1];
                    if (atom_s_bin[1] < 0) {atom_s_bin[1] += bin[i].tot_bins[1];}

                    long atom_s = bin[i].bin_list[atom_s_bin[0]][atom_s_bin[1]];

                    while (atom_s != -1)
                    {
                        long atom_s_reduced = atom_s%max_atoms; //this is used to access bin.next

                        if (atom_s != atom_p)
                        {
                            for (int m=0; m<Dim; m++)
                            {
                                dist_xy[m] = abs(configuration[atom_p].r[m] - configuration[atom_s].r[m]);

                                if (dist_xy[m] > L[m])
                                {
                                    cout << "atom lattice site beyond cell parameters " << dist_xy[m] << " " << L[m] << " " << atom_p << " " << atom_s << endl;

                                    //deleting memory allocations before exiting
                                    procedure = 2;
                                    construct = false;
                                    graphene_initialization(nucleation, configuration, bin, N, L, max_atoms, bins_x, bins_y, layers, procedure, construct);
                                    procedure = 1;
                                    construct = false;
                                    graphene_initialization(nucleation, NULL, NULL, N, L, 0, 0, 0, layers, procedure, construct);
                                    delete [] L;
                                    delete [] atom_p_bin;
                                    delete [] atom_s_bin;

                                    exit(1);
                                }


                                if (dist_xy[m] > L[m]/2) {dist_xy[m] = L[m]-dist_xy[m];}

                            }

                            if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < atom_min_dist2)
                            {
                                cout << "Two atoms closer than 1A, configuration not acceptable " << (pow(dist_xy[0],2)+pow(dist_xy[1],2)) << " " << atom_p << " " << atom_s << endl;

                                //deleting memory allocations before exiting
                                procedure = 2;
                                construct = false;
                                graphene_initialization(nucleation, configuration, bin, N, L, max_atoms, bins_x, bins_y, layers, procedure, construct);
                                procedure = 1;
                                construct = false;
                                graphene_initialization(nucleation, NULL, NULL, N, L, 0, 0, 0, layers, procedure, construct);
                                delete [] L;
                                delete [] atom_p_bin;
                                delete [] atom_s_bin;

                                exit(1);
                            }

                            if ((pow(dist_xy[0],2)+pow(dist_xy[1],2)) < pow((bond_length*1.732050808)-1,2)) {neighbor_list[neighbor_count] = atom_s; neighbor_count ++;}
                            if (neighbor_count > max_neighbors)
                            {
                                cout << "More than 3 atoms in the neighborhood (lattice constant - 1) of an atom, configuration not acceptable " << atom_p << " " << neighbor_list[0] << " " << neighbor_list[1] << " " << neighbor_list[2] << " " << neighbor_list[3] << endl;
                                cout << atom_p << " " << configuration[atom_p].r[0] << " " << configuration[atom_p].r[1] << " " << configuration[atom_p].bonded_atom[0] << " " << configuration[atom_p].bonded_atom[1] << " " << configuration[atom_p].bonded_atom[2] << " " << configuration[atom_p].grain << endl;
                                cout << neighbor_list[0] << " " << configuration[neighbor_list[0]].r[0] << " " << configuration[neighbor_list[0]].r[1] << " " << configuration[neighbor_list[0]].bonded_atom[0] << " " << configuration[neighbor_list[0]].bonded_atom[1] << " " << configuration[neighbor_list[0]].bonded_atom[2] << " " << configuration[neighbor_list[0]].grain  << endl;
                                cout << neighbor_list[1] << " " << configuration[neighbor_list[1]].r[0] << " " << configuration[neighbor_list[1]].r[1] << " " << configuration[neighbor_list[1]].bonded_atom[0] << " " << configuration[neighbor_list[1]].bonded_atom[1] << " " << configuration[neighbor_list[1]].bonded_atom[2] << " " << configuration[neighbor_list[1]].grain  << endl;
                                cout << neighbor_list[2] << " " << configuration[neighbor_list[2]].r[0] << " " << configuration[neighbor_list[2]].r[1] << " " << configuration[neighbor_list[2]].bonded_atom[0] << " " << configuration[neighbor_list[2]].bonded_atom[1] << " " << configuration[neighbor_list[2]].bonded_atom[2] << " " << configuration[neighbor_list[2]].grain  << endl;
                                cout << neighbor_list[3] << " " << configuration[neighbor_list[3]].r[0] << " " << configuration[neighbor_list[3]].r[1] << " " << configuration[neighbor_list[3]].bonded_atom[0] << " " << configuration[neighbor_list[3]].bonded_atom[1] << " " << configuration[neighbor_list[3]].bonded_atom[2] << " " << configuration[neighbor_list[3]].grain  << endl;


                                //deleting memory allocations before exiting
                                procedure = 2;
                                construct = false;
                                graphene_initialization(nucleation, configuration, bin, N, L, max_atoms, bins_x, bins_y, layers, procedure, construct);
                                procedure = 1;
                                construct = false;
                                graphene_initialization(nucleation, NULL, NULL, N, L, 0, 0, 0, layers, procedure, construct);
                                delete [] L;
                                delete [] atom_p_bin;
                                delete [] atom_s_bin;

                                exit(1);
                            }
                        }
                        atom_s = bin[i].next[atom_s_reduced];
                    }

                }
            }

        delete [] neighbor_list;
        delete [] atom_p_bin;
        delete [] dist_xy;
        delete [] atom_s_bin;
        }
    }



    //deleting the configuration_type and bin_type classes
    procedure = 2;
    construct = false;

    graphene_initialization(nucleation, configuration, bin, N, L, max_atoms, bins_x, bins_y, layers, procedure, construct);


    cout << "graphene_config_generator_deleted_configuration" << endl;

    //deleting memory  of nucleation_type class
    procedure = 1;
    construct = false;
    graphene_initialization(nucleation, NULL, NULL, N, L, 0, 0, 0, layers, procedure, construct);

    cout << "graphene_config_generator_deleted_nucleation" << endl;

    delete [] L;


    return 1;

}
