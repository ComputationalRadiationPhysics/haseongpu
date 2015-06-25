#include "types.hpp"
#include <iostream>

std::ostream& operator<<(std::ostream& out, const DeviceMode value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(NO_DEVICE_MODE);
        PROCESS_VAL(GPU_DEVICE_MODE);
        PROCESS_VAL(CPU_DEVICE_MODE);
    }
#undef PROCESS_VAL
    return out << s;
}

std::ostream& operator<<(std::ostream& out, const ParallelMode value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(NO_PARALLEL_MODE);
        PROCESS_VAL(THREADED_PARALLEL_MODE);
        PROCESS_VAL(MPI_PARALLEL_MODE);
        PROCESS_VAL(GRAYBAT_PARALLEL_MODE);
    }
#undef PROCESS_VAL
    return out << s;
}
