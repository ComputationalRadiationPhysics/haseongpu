#ifndef TESTDATA_TRANSPOSED_H
#define TESTDATA_TRANSPOSED_H

double host_p_in[]={
	-1,	0,	1,	-1,	1,
	-1,	0,	-1,	1,	1};

double host_n_x[] = {
	0,	-1,	1,	0,
	0.707106781186548, -0.707106781186548, 0.707106781186548, -0.707106781186548,
	0.707106781186548, 0.707106781186548, -0.707106781186548, -0.707106781186548};


double host_n_y[] = {
	-1,	0,	0,	1,
	-0.707106781186548, -0.707106781186548, 0.707106781186548, 0.707106781186548,
	0.707106781186548, -0.707106781186548, 0.707106781186548, -0.707106781186548};


int host_forbidden[] = {
	-1, -1, -1, -1,
	2,	2,	2,	2,
	1,	1,	1,	1};

int host_n_p[] = {
	0,3,2,4,
	0,3,2,4,
	2,0,4,3};

int host_neighbors[] = {
	-1,-1,-1,-1,
	1, 3, 0, 2,
	2, 0, 3, 1};

int host_t_in[] = {
	0,3,2,4,
	2,0,4,3,
	1,1,1,1};

double host_z_mesh = 1;
int host_size_t = 4;
int host_size_p = 5;
int host_mesh_z = 1;

double host_clad_abs = 5.5;
int host_clad_num = 3;

double host_N_tot = 2.76E20;
double host_beta_v[] = {0.00689445752513904,
	0.00703271177028374,
	0.00675620327999434,
	0.00689445752513904 };


double host_sigma_e = 2.4E-20;

double host_sigma_a = 1.16E-21;

int host_cell_type[] = { 1,1,1,1} ;



#endif