#include <logging.h>
#include <iostream>

std::ostream& dout(unsigned activation_level) {
  static nullstream dummy;

  if(!(activation_level & verbosity)){
    return dummy;
  }

  if(activation_level & V_NOLABEL){
    return std::cout;
  }

  if(activation_level & V_ERROR){
    return std::cerr << "\033[0;" << COLOR_ERROR << "m[ERROR] ";
  }

  if(activation_level & V_WARNING){
    return std::cerr << "\033[0;" << COLOR_WARN << "m[WARNING] ";
  }

  if(activation_level & V_INFO){
    //return std::cout << "\033[0;" << COLOR_INFO << "m[INFO] ";
    return std::cout << "\033[0" << "m[INFO] ";
  }

  if(activation_level & V_STAT){
    //return std::cout << "\033[0;" << COLOR_STATISTIC << "m[STATISTIC] ";
    return std::cout << "\033[0" << "m[STATISTIC] ";
  }

  if(activation_level & V_PROGRESS){
    return std::cout << "\033[0" << "m[PROGRESS] ";
  }

  if(activation_level & V_DEBUG){
    return std::cerr << "\033[0;" << COLOR_DEBUG << "m[DEBUG] ";
  }


  return std::cout;
}

