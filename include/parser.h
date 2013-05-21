#include <string>  /* string */
#include <iostream>
#include <fstream> /* ifstream */
#include <cstdlib> /* atof */
#include <vector> 

std::vector<double>* parse_beta_v(std::string root);
std::vector<unsigned>* parse_cell_type(std::string root);
std::vector<int>* parse_forbidden(std::string root);
std::vector<int>* parse_neighbors(std::string root);
std::vector<int>* parse_n_p(std::string root);
std::vector<double>* parse_n_x(std::string root);
std::vector<double>* parse_n_y(std::string root);
std::vector<int>* parse_p_in(std::string root);
std::vector<unsigned>* parse_t_in(std::string root);
float parse_clad_abs(std::string root);
float parse_clad_num(std::string root);
float parse_n_tot(std::string root);
float parse_sigma_a(std::string root);
float parse_sigma_e(std::string root);
unsigned parse_size_p(std::string root);
unsigned parse_size_t(std::string root);
unsigned parse_mesh_z(std::string root);

int parse(std::string location){
  
  if(location[location.size()-1] != '/')
    location.append("/");

  std::vector<double> * betas = parse_beta_v(location);
  std::vector<unsigned> * cell_types = parse_cell_type(location);
  std::vector<int> * forbidden  = parse_forbidden(location);
  std::vector<int> * neighbors  = parse_neighbors(location);
  std::vector<int> * n_p = parse_n_p(location);
  std::vector<double> * n_x = parse_n_x(location);
  std::vector<double> * n_y = parse_n_y(location);
  std::vector<int> * p_in = parse_p_in(location);
  std::vector<unsigned> * t_in = parse_t_in(location);
  
  float clad_abs = parse_clad_abs(location);
  float clad_num = parse_clad_num(location);
  float n_tot = parse_n_tot(location);
  float sigma_a = parse_sigma_a(location);
  float sigma_e = parse_sigma_e(location);
  unsigned size_p = parse_size_p(location);
  unsigned size_t = parse_size_t(location);
  unsigned mesh_z = parse_mesh_z(location);
  
  return 0;
} 

/* Parse beta_v.txt */
std::vector<double>* parse_beta_v(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<double>* betas = new std::vector<double>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return betas;
}

/* Parse cell_type.txt */
std::vector<unsigned>* parse_cell_type(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<unsigned>* cell_types = new std::vector<unsigned>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return cell_types;

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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return number;


}

/* Parse forbidden.txt */
std::vector<int>* parse_forbidden(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<int>* forbidden = new std::vector<int>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return forbidden;

}

/* Parse neighbors.txt */
std::vector<int>* parse_neighbors(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<int>* neighbors = new std::vector<int>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return neighbors;

}

/* Parse n_p.txt */
std::vector<int>* parse_n_p(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<int>* n_p = new std::vector<int>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return n_p;

}

/* Parse N_tot.txt */
float parse_n_tot(std::string root){
  std::string line;
  std::ifstream fileStream; 
  float number = 0;
  
  fileStream.open(root.append("N_tot.txt").c_str());
  if(fileStream.is_open()){
    if(fileStream.good()){
      std::getline(fileStream, line);
      number = atof(line.c_str());
    }

  }
  else{
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return number;


}

/* Parse n_x.txt */
std::vector<double>* parse_n_x(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<double>* n_x = new std::vector<double>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return n_x;
}

/* Parse n_y.txt */
std::vector<double>* parse_n_y(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<double>* n_y = new std::vector<double>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return n_y;
}

/* Parse p_in.txt */
std::vector<int>* parse_p_in(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<int>* p_in = new std::vector<int>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return p_in;

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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return number;


}


/* Parse t_in.txt */
std::vector<unsigned>* parse_t_in(std::string root){
  std::string line;
  std::ifstream fileStream; 
  std::vector<unsigned>* t_in = new std::vector<unsigned>();
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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return t_in;

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
    fprintf(stderr, "Can't open file %s", root.c_str());
    return NULL;
  }

  return number;


}
