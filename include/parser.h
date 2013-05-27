#include <string>  /* string */
#include <iostream>
#include <fstream> /* ifstream */
#include <cstdlib> /* atof */
#include <vector> 

int parse_beta_v(std::string root, std::vector<double>* betas);
int parse_cell_type(std::string root, std::vector<unsigned>* cell_types);
int parse_forbidden(std::string root, std::vector<int>* forbidden);
int parse_neighbors(std::string root, std::vector<int>* neighbors);
int parse_n_p(std::string root, std::vector<int>* n_p);
int parse_n_x(std::string root, std::vector<double>* );
int parse_n_y(std::string root, std::vector<double>* );
int parse_p_in(std::string root, std::vector<int>* );
int parse_t_in(std::string root, std::vector<unsigned>* );
float parse_clad_abs(std::string root);
float parse_clad_num(std::string root);
float parse_n_tot(std::string root);
float parse_sigma_a(std::string root);
float parse_sigma_e(std::string root);
unsigned parse_size_p(std::string root);
unsigned parse_size_t(std::string root);
unsigned parse_mesh_z(std::string root);

int parse(std::string location, 
	  std::vector<double> * betas,
	  std::vector<double> * n_x,
	  std::vector<double> * n_y,
	  std::vector<unsigned> * cell_types,
	  std::vector<unsigned> * t_in,
	  std::vector<int> * forbidden,
	  std::vector<int> * neighbors,
	  std::vector<int> * n_p,
	  std::vector<int> * p_in,
	  float *clad_abs,
	  float *clad_num,
	  float *n_tot,
	  float *sigma_a,
	  float *sigma_e,
	  unsigned *size_p,
	  unsigned *size_t,
	  unsigned *mesh_z){
  
  if(location[location.size()-1] != '/')
    location.append("/");

  // Parse every file for its own
  if(parse_beta_v(location, betas)) return 1;
  if(parse_cell_type(location, cell_types)) return 1;
  if(parse_forbidden(location, forbidden)) return 1;
  if(parse_neighbors(location, neighbors)) return 1;
  if(parse_n_p(location, n_p)) return 1;
  if(parse_n_x(location, n_x)) return 1;
  if(parse_n_y(location, n_y)) return 1;
  if(parse_p_in(location, p_in)) return 1;
  if(parse_t_in(location, t_in)) return 1;
  
  (*clad_abs) = parse_clad_abs(location);
  (*clad_num) = parse_clad_num(location);
  (*n_tot) = parse_n_tot(location);
  (*sigma_a) = parse_sigma_a(location);
  (*sigma_e) = parse_sigma_e(location);
  (*size_p) = parse_size_p(location);
  (*size_t) = parse_size_t(location);
  (*mesh_z) = parse_mesh_z(location);

  if(!((*clad_abs) && *(clad_num)  && (*n_tot) && (*sigma_a) && *(sigma_e) && (*size_p) && (*size_t) && (*mesh_z)))
    return 1;

  return 0;
} 

/* Parse beta_v.txt */
int parse_beta_v(std::string root, std::vector<double>* betas){
  std::string line;
  std::ifstream fileStream; 
  double number = 0;
  
  fileStream.open(root.append("beta_v.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
      betas->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;
}

/* Parse cell_type.txt */
int parse_cell_type(std::string root, std::vector<unsigned>* cell_types){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("cell_type.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
      cell_types->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;

}

/* Parse clad_abs.txt */
float parse_clad_abs(std::string root){
  std::string line;
  std::ifstream fileStream; 
  float number = 0;
  
  fileStream.open(root.append("clad_abs.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return  0.0;
  }

  return number;


}

/* Parse clad_num.txt */
float parse_clad_num(std::string root){
  std::string line;
  std::ifstream fileStream; 
  float number = 0;
  
  fileStream.open(root.append("clad_num.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 0.0;
  }

  return number;


}

/* Parse forbidden.txt */
int parse_forbidden(std::string root, std::vector<int>* forbidden){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("forbidden.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
      forbidden->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;

}

/* Parse neighbors.txt */
int parse_neighbors(std::string root, std::vector<int>* neighbors){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("neighbors.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
      neighbors->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;

}

/* Parse n_p.txt */
int parse_n_p(std::string root, std::vector<int>* n_p){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("n_p.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
      n_p->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;

}

/* Parse N_tot.txt */
float parse_n_tot(std::string root){
  std::string line;
  std::ifstream fileStream; 
  float number = 0;
  
  fileStream.open(root.append("n_tot.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 0.0;
  }

  return number;


}

/* Parse n_x.txt */
int parse_n_x(std::string root, std::vector<double>* n_x){
  std::string line;
  std::ifstream fileStream; 
  double number = 0;
  
  fileStream.open(root.append("n_x.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
      n_x->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;
}

/* Parse n_y.txt */
int parse_n_y(std::string root, std::vector<double>* n_y){
  std::string line;
  std::ifstream fileStream; 
  double number = 0;
  
  fileStream.open(root.append("n_y.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
      n_y->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;
}

/* Parse p_in.txt */
int parse_p_in(std::string root, std::vector<int>* p_in){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("p_in.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
      p_in->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;

}

/* Parse sigma_a.txt */
float parse_sigma_a(std::string root){
  std::string line;
  std::ifstream fileStream; 
  float number = 0;
  
  fileStream.open(root.append("sigma_a.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 0.0;
  }

  return number;


}

/* Parse sigma_e.txt */
float parse_sigma_e(std::string root){
  std::string line;
  std::ifstream fileStream; 
  float number = 0;
  
  fileStream.open(root.append("sigma_e.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 0.0;
  }

  return number;


}

/* Parse size_p.txt */
unsigned parse_size_p(std::string root){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("size_p.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 0;
  }

  return number;


}

/* Parse size_t.txt */
unsigned parse_size_t(std::string root){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("size_t.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 0;
  }

  return number;


}


/* Parse t_in.txt */
int parse_t_in(std::string root, std::vector<unsigned>* t_in){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("t_in.txt").c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
      t_in->push_back(number);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 1;
  }

  return 0;

}

/* Parse mesh_z.txt */
unsigned parse_mesh_z(std::string root){
  std::string line;
  std::ifstream fileStream; 
  unsigned number = 0;
  
  fileStream.open(root.append("mesh_z.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atoi(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", root.c_str());
    return 0;
  }

  return number;


}
