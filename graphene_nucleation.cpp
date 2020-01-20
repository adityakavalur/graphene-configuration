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

using namespace std;


void graphene_nucleation(nucleation_type *nucleation, long N, double *L, const char *angledistri_file)
{
    /*long stime=time(NULL);
    long ctime=stime;
    while (stime==ctime) {
       ctime=time(NULL);
       cout <<stime<<" "<<ctime<<endl;
    }*/


    srand (time(NULL));//check best place to put this
    long site_count = 0;//initializing the nucleation site count
    double temp_x, temp_y;
    double dist_x, dist_y;
    bool accept_site;
    double temp_distance2;
    int temp_rand1,temp_rand2;


    //generating nucleation site co-ordinates
    while (site_count < N)
    {
        //generating random x, y co-ordinates between 0 - Lx, 0 - Ly respectively
        temp_rand1 = rand();
        temp_rand2 = rand();

        /*
        //
        //
        //

        Change back to randomization
        */
        temp_x = L[0] * temp_rand1 / double(RAND_MAX);//check accuracy should be more than 10^-6 also check range
        temp_y = L[1] * temp_rand2 / double(RAND_MAX);

        //if (site_count==0) {temp_x = 100.0;temp_y = 18.75;}
        //else {temp_x=101.2; temp_y = 56.25;}
        //if (site_count==0 || site_count==1) {temp_y = 18.75;}
        //else {temp_y=56.25;}


        //cout <<"Random numbers = "<<temp_rand1<<" "<<temp_rand2<<endl;
        cout << "Nucleation sites = " << temp_x << "," << temp_y << endl;

        if (temp_x<=0.0 || temp_x>=L[0]) {
           cout <<"Out of range in the x direction: "<<temp_x<<endl;
           exit(1);
        }
        if (temp_y<=0.0 || temp_y>=L[1]) {
           cout <<"Out of range in the y direction: "<<temp_y<<endl;
           exit(1);
        }

        //cout <<site_count<<" "<<temp_x<<" "<<temp_y<<endl;

        accept_site = true;//re-initiating for a new random site
        {
            long j=0;
            while (j != site_count)//to check with already created nucleation sites
            {
                dist_x = abs(temp_x - nucleation[j].xy[0]);
                dist_y = abs(temp_y - nucleation[j].xy[1]);

                if (dist_x > L[0]/2) {dist_x = L[0] - dist_x;}
                if (dist_y > L[1]/2) {dist_y = L[1] - dist_y;}
                temp_distance2 = pow(dist_x, 2.0) + pow(dist_y, 2.0);
                if (temp_distance2 < nucl_min_dist2)
                {
                    accept_site = false;
                }

                //cout <<"+++"<<j<<" "<<sqrt(temp_distance2)<<" "<<accept_site<<endl;

                j++;
            }

        }

        if (accept_site == true)
        {
            nucleation[site_count].xy[0] = temp_x;
            nucleation[site_count].xy[1] = temp_y;
            site_count ++;
        }

    }

    cout <<site_count<<" "<<N<<endl;
    //exit(1);


    //mapping the closest neighbors in each quadrant to each nucleation point and giving an allowance of +/- 5 deg
    //array number is storing neighbor information for quadrant number - 1
    int quadrant_number;

    for (long i=0; i<N; i++)
    {

        double **temp_neigh_dist2;

        temp_neigh_dist2 = new double *[quadrants];
        for (int j=0; j<quadrants; j++)
        {
            temp_neigh_dist2[j] = new double [nucl_neigh_count];
        }


        for (int j=0; j<quadrants; j++)
        {
            for (int k=0; k< nucl_neigh_count; k++)
            {
                temp_neigh_dist2[j][k] = L[0]*L[1];
            }
        }
        //nucl_neigh_count declared in constants.h

        //mapping the closest nucl_neigh_count neighbors in each quadrant
       // cout << "line 116" << endl;
        for (long j=0; j<N; j++)
        {

            if (i!=j)
            {   //cout << "line 121" << endl;
                int *dir_xy;
                double *temp_dist, temp_dist2;
                temp_dist = new double [Dim];
                dir_xy = new int [Dim];
                for (int k=0; k<Dim; k++)
                {   //cout << "line 127" << endl;
                    temp_dist[k] = (nucleation[j].xy[k] - nucleation[i].xy[k]); //cout << "line 145: " << k << " " << temp_dist[k] << endl;
                    if (temp_dist[k] >= L[k]/2) {dir_xy[k] = -1;}
                    else if (temp_dist[k] <= -L[k]/2) {dir_xy[k] = 1;}
                    else if (temp_dist[k] <= 0 && temp_dist[k] > -L[k]/2) {dir_xy[k] = -1;}
                    else if (temp_dist[k] >= 0 && temp_dist[k] < L[k]/2) {dir_xy[k] = 1;}
                }

                //cout << "line 135: "<< dir_xy[0] << " " << dir_xy[1] << endl;
                if (dir_xy[0] == 1 && dir_xy[1] == 1) {quadrant_number = 1;}

                if (dir_xy[0] == -1 && dir_xy[1] == 1) {quadrant_number = 2;}

                if (dir_xy[0] == -1 && dir_xy[1] == -1) {quadrant_number = 3;}

                if (dir_xy[0] == 1 && dir_xy[1] == -1) {quadrant_number = 4;}
                cout << "line 155" << endl;
                for (int k=0; k<Dim; k++)
                {
                    temp_dist[k] = abs(temp_dist[k]);
                    if(temp_dist[k] > L[k]/2) {temp_dist[k] = L[k] - temp_dist[k];}
                }
                cout << "line 161: " << nucl_neigh_count << endl;
                temp_dist2 = pow(temp_dist[0],2.0) + pow(temp_dist[1],2.0);

                bool neighbor_update;
                neighbor_update = false;

                for (int k=0; k<nucl_neigh_count; k++)
                {//cout << "line 168: " << endl;//<< temp_dist2 << " " << quadrant_number << endl;//<< ";"<< temp_neigh_dist2[quadrant_number-1][k] << endl;//<< ";" << neighbor_update << endl;
                    if (temp_dist2 < temp_neigh_dist2[quadrant_number-1][k] && neighbor_update == false)
                    {   cout << "line 158" << endl;
                        for (int l=nucl_neigh_count-2; l>k-1; l--)
                        {   cout << "line 160 " << quadrant_number << " " << i << " " << l  << " " << nucl_neigh_count << endl;
                            nucleation[i].neighbors[quadrant_number-1][l+1] = nucleation[i].neighbors[quadrant_number-1][l];cout << "line 161" << endl;
                            temp_neigh_dist2[quadrant_number-1][l+1] = temp_neigh_dist2[quadrant_number-1][l];cout << "line 162" << endl;
                        }cout << "line 163" << endl;
                        nucleation[i].neighbors[quadrant_number-1][k] = j;
                        temp_neigh_dist2[quadrant_number-1][k] = temp_dist2;
                        neighbor_update = true;cout << "line 166" << endl;
                    }
                }

                delete [] temp_dist;
            }
        }

        delete [] temp_neigh_dist2;

    }

    cout << "line 170" << endl;


    //checking if any of nucleation neighbors can be ignored
    int deletion_counter = 0;
    for (long i=0; i<N; i++)
    {
        //cout << i << endl;
        double **slope_angle, **intercept, **perpen_slope_angle;
        slope_angle = new double *[quadrants];
        intercept = new double *[quadrants];
        perpen_slope_angle = new double *[quadrants];
        for (int j=0; j<quadrants; j++)
        {
            slope_angle[j] = new double [nucl_neigh_count];
            intercept[j] = new double [nucl_neigh_count];
            perpen_slope_angle[j] = new double [nucl_neigh_count];
        }
        //calculating line equations of the neighbors to the primary atom and perpendicular slopes
        for (int j=0; j<quadrants; j++)
        {//cout << "line186" << endl;
            for (int k=0; k<nucl_neigh_count; k++)
            {//cout << "line 188" << endl;
            long atom_s = nucleation[i].neighbors[j][k]; //cout << atom_s << endl;
            if (atom_s != -1)
            {
                slope_angle[j][k] = atan2((nucleation[atom_s].xy[1] - nucleation[i].xy[1]),(nucleation[atom_s].xy[0] - nucleation[i].xy[0]));
                intercept[j][k] = nucleation[atom_s].xy[1] - (tan(slope_angle[j][k])*nucleation[atom_s].xy[0]);
                perpen_slope_angle[j][k] = slope_angle[j][k] + 90.0;
            }

            }
        }
        for (int j=0; j<quadrants; j++)
        {//cout << "line 196" << endl;
           for (int k=0; k<nucl_neigh_count; k++)
           {//cout << "line 198" << endl;
                long atom_p = nucleation[i].neighbors[j][k];
                bool delete_neighbor;
                delete_neighbor = false;
                if (atom_p != -1)
                {//cout << "line 203" << endl;
                    double dist_atom_p2;
                    dist_atom_p2 = pow(nucleation[atom_p].xy[0] - nucleation[i].xy[0], 2.0) + pow(nucleation[atom_p].xy[1] - nucleation[i].xy[1], 2.0);
                    for (int l=j-1; l<j+2; l++)
                    {//cout << "line 207" << " " << i <<  endl;
                        int actual_quadrant = l;
                        actual_quadrant = actual_quadrant + quadrants;
                        actual_quadrant = actual_quadrant%quadrants;
                        if (actual_quadrant == 0) {actual_quadrant = quadrants;}
                        //cout << actual_quadrant << " " << i << endl;
                        for (int m=0; m<nucl_neigh_count; m++)
                        {//cout << "line 214 " << actual_quadrant << " " << i << endl;
                            long atom_s = nucleation[i].neighbors[actual_quadrant-1][m]; //cout << "line 228 " << atom_s << endl;
                            if (atom_p != atom_s && atom_s != -1)
                            {  // cout << "line 230 " << atom_s << endl;
                                double dist_atom_s2;
                                dist_atom_s2 = pow(nucleation[atom_s].xy[0] - nucleation[i].xy[0], 2.0) + pow(nucleation[atom_s].xy[1] - nucleation[i].xy[1], 2.0);

                                if (dist_atom_p2 >= dist_atom_s2*4.0)
                                {//cout << "line 220" << endl;
                                    double temp_intercept;
                                    temp_intercept = nucleation[atom_s].xy[1] - (tan(perpen_slope_angle[j][k])*nucleation[atom_s].xy[0]);

                                    double intersec_x, intersec_y;
                                    intersec_x = (temp_intercept - intercept[j][k])/(tan(slope_angle[j][k])-tan(perpen_slope_angle[j][k]));
                                    intersec_y = tan(slope_angle[j][k])*intersec_x + intercept[j][k];

                                    double intersec_atom_p2, intersec_atom_s2;
                                    intersec_atom_p2 = pow(intersec_x - nucleation[atom_p].xy[0],2.0) + pow(intersec_y - nucleation[atom_p].xy[1],2.0);
                                    intersec_atom_s2 = pow(intersec_x - nucleation[atom_s].xy[0],2.0) + pow(intersec_y - nucleation[atom_s].xy[1],2.0);

                                    if (intersec_atom_p2 > intersec_atom_s2*2.25)
                                    {
                                        delete_neighbor = true;
                                    }

                                }
                            }

                        }
                    }
                }
           if (delete_neighbor == true) {nucleation[i].neighbors[j][k] = -1;}
           }
        }
       // cout << "line 261" << endl;
        for (int j=0; j<quadrants; j++)
        {
            delete [] slope_angle[j];
            delete [] intercept[j];
            delete [] perpen_slope_angle[j];
        }
        delete [] slope_angle;
        delete [] intercept;
        delete [] perpen_slope_angle;
    }
    cout << "starting angle generation" << endl;
    //generating random mis-orientation angles between 0-60 deg
    for (long i=0; i<N; i++)
    {
        //nucleation[i].angle = 0;
        cout << i << endl;
        bool angle_accept = false;
        double temp_angle;
        long counter = 0;//to check number of iterations and exit if the angle generation is stuck in infinite loop
        while (angle_accept != true)
        {
            angle_accept = true;

            temp_angle = 60 * (rand() / double(RAND_MAX));
            while (temp_angle<3.0 || temp_angle>57.0) {temp_angle = 60 * (rand() / double(RAND_MAX));}
            //************ check if you want to assign pre-determined angle*************************************************************************
            //temp_angle = 5.0;

            //cout << "line 285 " << endl;
            counter ++;
            if (counter > 1000) {cout << "angle generation failed, re-initiate prog or reduce min_angle between grains" << endl; exit(1);}

            for (int j=0; j<quadrants; j++)
            {
                for (int k=0; k<nucl_neigh_count; k++)
                {
                    long atom_s = nucleation[i].neighbors[j][k];
                    if (atom_s != -1)
                    {
                        if(nucleation[atom_s].angle != -1)
                        {
                            double angle_diff = abs(temp_angle - nucleation[atom_s].angle) ;
                            if (angle_diff >= 60.0) {angle_diff = angle_diff - 60.0;}
                            if (angle_diff >= 30.0) {angle_diff = 60.0 - angle_diff;}
                            if (angle_diff < min_angle) {angle_accept = false;}
                            //**************put this back when you want to go back to minimum angle*******************//
                            //*******remove it when you want all grains to have same angle**************//
                        }
                    }
                }


            }
        }
        nucleation[i].angle = temp_angle;
        if (i==0) {nucleation[i].angle=0.0;}
        //if (i==1) {nucleation[i].angle=59.5;cout << "line 324" << endl;}
        //else {nucleation[i].angle=30.0;}
        cout << nucleation[i].angle << endl;
    }

    int bin_width = 5;
    int num_bins = 60/bin_width;
    long * angle_bin;//bins for angle distribution
    angle_bin = new long [num_bins];// each bin is of width 5 degree

    for (int i=0; i<num_bins; i++)
    {
        angle_bin[i] = 0;
    }

    for (long i=0; i<N; i++)
    {
        double current_angle = nucleation[i].angle;
        cout << nucleation[i].angle << endl;
        int temp_bin = floor(current_angle/bin_width);

        angle_bin[temp_bin] ++;
    }


    //normalizing the bin distribution data
    for (long i=0; i<num_bins; i++)
    {
        angle_bin[i] = angle_bin[i]/N;
    }


    ofstream myfile;
    myfile.open (angledistri_file);

    site_count=0;
    for (int i=0; i<num_bins; i++)
    {
        myfile << (i+0.5)*bin_width << " " << angle_bin[i] << endl;
        site_count+=(N*angle_bin[i]);
    }
    cout <<"site_count="<<site_count<<endl;

    myfile.close();

    delete [] angle_bin;

}
