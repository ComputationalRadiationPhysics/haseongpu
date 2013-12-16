VTK=vtk-5.8
#VTK=vtk-5.10
INCLUDES=/usr/include/$(VTK)
LIBS=/usr/lib/$(VTK)

all:
	g++ convertVTKtoVTU.cpp -I $(INCLUDES)  -L $(LIBS)  -l vtkIO -o convertVTKtoVTU
