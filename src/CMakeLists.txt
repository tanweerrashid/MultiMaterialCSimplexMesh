
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(MultiMaterialCSimplexMesh)

FIND_PACKAGE (VTK)
IF (VTK_FOUND)
         INCLUDE(${VTK_USE_FILE})
ENDIF(VTK_FOUND)



INCLUDE_DIRECTORIES(${MultiMaterialCSimplexMesh_SOURCE_DIR})

ADD_EXECUTABLE(MultiMaterialCSimplexMesh main.cpp CSimplexMeshConverter.cpp Deformation.cpp MMCSimplexSurf.cpp MultiMaterialDeformation.cpp SimplexMesh.cpp SimplexSurf.cpp StdFunctions.cpp)

TARGET_LINK_LIBRARIES (MultiMaterialCSimplexMesh ${VTK_LIBRARIES}) 



