#include <cuda_utils.h>

#include <logging.h>
#include <vector>


/** 
 * @brief Queries for devices on the running mashine and collects
 *        them on the devices array. Set the first device in this 
 *        array as computation-device. On Errors the programm will
 *        be stoped by exit(). 
 * 
 * @param maxGpus max. devices which should be allocated
 * @return vector of possible devices
 */
std::vector<unsigned> getFreeDevices(unsigned maxGpus){
  cudaDeviceProp prop;
  int minMajor = MIN_COMPUTE_CAPABILITY_MAJOR;
  int minMinor = MIN_COMPUTE_CAPABILITY_MINOR;
  int count;
  std::vector<unsigned> devices;

  // Get number of devices
  CUDA_CHECK_RETURN( cudaGetDeviceCount(&count));

  // Check devices for compute capability and if device is busy
  unsigned devicesAllocated = 0;
  for(int i=0; i < count; ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, i) );
    if( (prop.major > minMajor) || (prop.major == minMajor && prop.minor >= minMinor) ){
      cudaSetDevice(i);
      int* occupy; //TODO: occupy gets allocated, but never cudaFree'd -> small memory leak!
      if(cudaMalloc((void**) &occupy, sizeof(int)) == cudaSuccess){
        devices.push_back(i);
        devicesAllocated++;
        if(devicesAllocated == maxGpus)
          break;

      }

    }

  }
  // Exit if no device was found
  if(devices.size() == 0){
    dout(V_ERROR) << "None of the free CUDA-capable devices is sufficient!" << std::endl;
    exit(1);
  }

  // Print device information
  cudaSetDevice(devices.at(0));
  dout(V_INFO) << "Found " << int(devices.size()) << " available CUDA devices with Compute Capability >= " << minMajor << "." << minMinor << "):" << std::endl;
  for(unsigned i=0; i<devices.size(); ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, devices[i]) );
    dout(V_INFO) << "[" << devices[i] << "] " << prop.name << " (Compute Capability " << prop.major << "." << prop.minor << ")" << std::endl;
  }

  return devices;

}
