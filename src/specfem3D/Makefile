#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  2 . 1
#               ---------------------------------------
#
#          Main authors: Dimitri Komatitsch and Jeroen Tromp
#    Princeton University, USA and University of Pau / CNRS / INRIA
# (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
#                            April 2011
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================
#
# United States Government Sponsorship Acknowledged.
#

# src/specfem3D/Makefile.  Generated from Makefile.in by configure.

# CUDA
# with configure: ./configure --with-cuda CUDA_LIB=.. CUDA_INC=.. MPI_INC=..

# default cuda libraries
# runtime library -lcudart needed, others are optional -lcuda -lcublas
#CUDA_LIBS = -lcudart   
CUDA_LIBS = 

CUDA_LIB_LOCATION = 
CUDA_LINK = $(CUDA_LIB_LOCATION) $(CUDA_LIBS)
CUDA_INC =  -I../../
MPI_INC = 

#NVCC = nvcc
NVCC = /usr/projects/hpcsoft/toss2/common/intel/2015.5.223/composer_xe_2015.5.223/bin/intel64/icc


#NVCC_FLAGS_BASE = $(CUDA_INC) $(MPI_INC) $(COND_MPI_CPPFLAGS) 
##NVCC_FLAGS = $(NVCC_FLAGS_BASE) -dc -DCUDA -gencode=arch=compute_35,code=sm_35 
#NVCC_FLAGS = $(NVCC_FLAGS_BASE) -DCUDA -DUSE_OLDER_CUDA4_GPU -gencode=arch=compute_20,code=sm_20 

##NVCCLINK_BASE = $(NVCC) $(CUDA_INC) $(MPI_INC) $(COND_MPI_CPPFLAGS) -DCUDA
##NVCCLINK = $(NVCCLINK_BASE) -dlink -gencode=arch=compute_35,code=sm_35
#NVCCLINK = $(NVCCLINK_BASE) -DUSE_OLDER_CUDA4_GPU -gencode=arch=compute_20,code=sm_20

NVCC_FLAGS = $(MPI_INC) $(COND_MPI_CPPFLAGS)
NVCCLINK = $(NVCC) $(NVCC_FLAGS)

# OpenMP
# with configure: ./configure --with-openmp FLAGS_CHECK="-openmp .." OPENMP_LIB=.. 
#OPENMP_LIBS = $(OPENMP_LIB)
OPENMP_LIBS = 
#COND_OPENMP_FFLAGS = -DOPENMP_MODE
COND_OPENMP_FFLAGS =


FC = /usr/projects/hpcsoft/toss2/common/intel/2015.5.223/composer_xe_2015.5.223/bin/intel64/ifort
FCFLAGS = #-g
MPIFC = mpif90
MPILIBS = 
FLAGS_CHECK = -O3 -check nobounds -DFORCE_VECTORIZATION -xHost -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage
FCFLAGS_f90 = 


SHARED = ../shared/
CUDAD = ../cuda/

# E : executables directory
E = ../../bin
# O : objects directory
O = ../../obj/spec
# L : libraries directory
L = ../../lib
# Output files directory
OUTPUT=../../OUTPUT_FILES

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK) $(COND_OPENMP_FFLAGS) -I$(SHARED)
MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK) $(COND_OPENMP_FFLAGS) -I$(SHARED)
FCLINK = $(MPIFCCOMPILE_CHECK)
#FCLINK = $(FCCOMPILE_CHECK)

CC = /usr/projects/hpcsoft/toss2/common/intel/2015.5.223/composer_xe_2015.5.223/bin/intel64/icc
CFLAGS = -g -O2 $(CPPFLAGS)
CPPFLAGS =  -I/usr/include/scotch $(COND_MPI_CPPFLAGS)
COND_MPI_CPPFLAGS = -DWITH_MPI
#COND_MPI_CPPFLAGS =

AR = ar
ARFLAGS = cru
RANLIB = ranlib

libspecfem_a_OBJECTS = \
	$O/assemble_MPI_scalar.shared.o \
	$O/check_mesh_resolution.shared.o \
	$O/comp_source_time_function.o \
	$O/compute_adj_source_frechet.o \
	$O/compute_arrays_source.o \
	$O/create_name_database.shared.o \
	$O/create_serial_name_database.shared.o \
	$O/define_derivation_matrices.shared.o \
	$O/detect_surface.shared.o \
	$O/exit_mpi.shared.o \
	$O/force_ftz.cc.o \
	$O/get_attenuation_model.shared.o \
	$O/get_cmt.shared.o \
	$O/get_element_face.shared.o \
	$O/get_force.shared.o \
	$O/get_jacobian_boundaries.shared.o \
	$O/get_shape3D.shared.o \
	$O/get_value_parameters.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/lagrange_poly.shared.o \
	$O/locate_receivers.o \
	$O/locate_source.o \
	$O/multiply_arrays_source.o \
	$O/netlib_specfun_erf.shared.o \
	$O/param_reader.cc.o \
	$O/prepare_assemble_MPI.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/save_header_file.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$O/write_VTK_data.shared.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)

CUDA_OBJECTS = \
	$O/check_fields_cuda.cuda.o \
	$O/compute_add_sources_acoustic_cuda.cuda.o \
	$O/compute_add_sources_viscoelastic_cuda.cuda.o \
	$O/compute_coupling_cuda.cuda.o \
	$O/compute_forces_acoustic_cuda.cuda.o \
	$O/compute_forces_viscoelastic_cuda.cuda.o \
	$O/compute_kernels_cuda.cuda.o \
	$O/compute_stacey_acoustic_cuda.cuda.o \
	$O/compute_stacey_viscoelastic_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$O/it_update_displacement_cuda.cuda.o \
	$O/noise_tomography_cuda.cuda.o \
	$O/prepare_mesh_constants_cuda.cuda.o \
	$O/save_and_compare_cpu_vs_gpu.cudacc.o \
	$O/transfer_fields_cuda.cuda.o \
	$O/write_seismograms_cuda.cuda.o \
	$(EMPTY_MACRO)

CUDA_STUBS = \
	$O/specfem3D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)

CUDA_DEVICE_OBJ = \
        $O/cuda_device_obj.o \
        $(EMPTY_MACRO)


# solver objects - no statically allocated arrays anymore
SOLVER_ARRAY_OBJECTS = \
	$O/specfem3D_par.o \
	$O/pml_par.o \
	$O/assemble_MPI_vector.o \
	$O/fault_solver_common.o \
	$O/fault_solver_dynamic.o \
	$O/fault_solver_kinematic.o \
        $O/gravity_perturbation.o \
	$O/compute_add_sources_acoustic.o \
	$O/compute_add_sources_viscoelastic.o \
	$O/compute_add_sources_poroelastic.o \
	$O/compute_boundary_kernel.o \
	$O/compute_coupling_acoustic_el.o \
	$O/compute_coupling_acoustic_po.o \
	$O/compute_coupling_viscoelastic_ac.o \
	$O/compute_coupling_viscoelastic_po.o \
	$O/compute_coupling_poroelastic_ac.o \
	$O/compute_coupling_poroelastic_el.o \
	$O/compute_forces_acoustic_calling_routine.o \
	$O/compute_forces_acoustic_noDev.o \
	$O/compute_forces_acoustic_Dev.o \
	$O/compute_forces_viscoelastic_calling_routine.o \
	$O/compute_forces_viscoelastic_Dev.o \
	$O/compute_forces_viscoelastic_noDev.o \
	$O/compute_forces_poro_fluid_part.o \
	$O/compute_forces_poroelastic_calling_routine.o \
	$O/compute_forces_poro_solid_part.o \
	$O/compute_gradient.o \
	$O/compute_interpolated_dva.o \
	$O/compute_kernels.o \
	$O/compute_stacey_acoustic.o \
	$O/compute_stacey_viscoelastic.o \
	$O/compute_stacey_poroelastic.o \
	$O/create_color_image.o \
	$O/detect_mesh_surfaces.o \
	$O/finalize_simulation.o \
	$O/initialize_simulation.o \
	$O/iterate_time.o \
	$O/make_gravity.o \
	$O/noise_tomography.o \
	$O/pml_allocate_arrays.o \
	$O/pml_output_VTKs.o \
	$O/pml_compute_accel_contribution.o \
	$O/pml_compute_memory_variables.o \
	$O/prepare_timerun.o \
	$O/program_specfem3D.o \
	$O/read_mesh_databases.o \
	$O/save_adjoint_kernels.o \
	$O/setup_GLL_points.o \
	$O/setup_movie_meshes.o \
	$O/setup_sources_receivers.o \
	$O/specfem3D.o \
	$O/write_movie_output.o \
	$O/write_output_ASCII.o \
	$O/write_output_SU.o \
	$O/write_seismograms.o \
	$(EMPTY_MACRO)


MODEL_UPD_OBJECTS = \
	$O/specfem3D_par.o \
	$O/pml_par.o \
	$O/model_update.o \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/get_value_parameters.shared.o \
	$O/get_attenuation_model.shared.o \
	$O/gll_library.shared.o \
	$O/initialize_simulation.o \
	$O/param_reader.cc.o \
	$O/read_mesh_databases.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/save_external_bin_m_up.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)


SUM_KERNELS_OBJECTS = \
	$O/sum_kernels.o \
	$O/exit_mpi.shared.o \
	$O/get_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)


# objects toggled between the parallel and serial version
COND_MPI_OBJECTS = $O/parallel.o
#COND_MPI_OBJECTS = $O/serial.o

# objects toggled between openmp and non-openmp version
#COND_OPENMP_OBJECTS = $O/older_not_maintained_compute_forces_viscoelastic_Dev_openmp.openmp.o
COND_OPENMP_OBJECTS =

LIBSPECFEM = $L/libspecfem.a

# objects for the pure Fortran version
#XSPECFEM_OBJECTS = $(SOLVER_ARRAY_OBJECTS) $(LIBSPECFEM) $(CUDA_OBJECTS)
XSPECFEM_OBJECTS = $(SOLVER_ARRAY_OBJECTS) $(LIBSPECFEM) $(CUDA_STUBS)

####
#### targets
####

# default targets for the pure Fortran version
DEFAULT = \
	specfem3D \
	combine_vol_data \
	combine_surf_data \
	convolve_source_timefunction \
	smooth_vol_data \
	sum_kernels \
	model_update \
	$(EMPTY_MACRO)

default: $(DEFAULT)

all: clean default

specfem3D: xspecfem3D


####
#### rules for executables
####


# rules for the pure Fortran version
# solver also depends on values from mesher
xspecfem3D: $(XSPECFEM_OBJECTS) $(COND_MPI_OBJECTS) $(COND_OPENMP_OBJECTS)
##	${NVCCLINK} -o $(CUDA_DEVICE_OBJ) $(CUDA_OBJECTS)
##	${FCLINK} -o ${E}/xspecfem3D $(XSPECFEM_OBJECTS) $(CUDA_DEVICE_OBJ) $(COND_MPI_OBJECTS) $(MPILIBS) $(COND_OPENMP_OBJECTS) $(OPENMP_LIBS) $(CUDA_LINK)
#	${FCLINK} -o ${E}/xspecfem3D $(XSPECFEM_OBJECTS) $(COND_MPI_OBJECTS) $(MPILIBS) $(COND_OPENMP_OBJECTS) $(OPENMP_LIBS) $(CUDA_LINK)
	${FCLINK} -o ${E}/xspecfem3D $(XSPECFEM_OBJECTS) $(COND_MPI_OBJECTS) $(MPILIBS) $(COND_OPENMP_OBJECTS) $(OPENMP_LIBS) $(CUDA_LINK)


convolve_source_timefunction: xconvolve_source_timefunction
create_movie_shakemap_AVS_DX_GMT: xcreate_movie_shakemap_AVS_DX_GMT
combine_vol_data: xcombine_vol_data
combine_surf_data: xcombine_surf_data
smooth_vol_data: xsmooth_vol_data
sum_kernels: xsum_kernels
model_update: xmodel_update

xconvolve_source_timefunction: $O/convolve_source_timefunction.shared.o
	${FCCOMPILE_CHECK} -o  ${E}/xconvolve_source_timefunction $O/convolve_source_timefunction.shared.o

xcreate_movie_shakemap_AVS_DX_GMT: $O/create_movie_shakemap_AVS_DX_GMT.shared.o $(LIBSPECFEM) $(OUTPUT)/surface_from_mesher.h
	${FCCOMPILE_CHECK} -o  ${E}/xcreate_movie_shakemap_AVS_DX_GMT $O/create_movie_shakemap_AVS_DX_GMT.shared.o $(LIBSPECFEM) -I$(OUTPUT)

xcombine_vol_data: $O/combine_vol_data.shared.o $O/write_c_binary.cc.o $O/read_parameter_file.shared.o $O/read_value_parameters.shared.o $O/get_value_parameters.shared.o $O/param_reader.cc.o
	${FCCOMPILE_CHECK} -o  ${E}/xcombine_vol_data  $O/combine_vol_data.shared.o $O/write_c_binary.cc.o $O/read_parameter_file.shared.o $O/read_value_parameters.shared.o $O/get_value_parameters.shared.o $O/param_reader.cc.o

xcombine_surf_data: $O/combine_surf_data.shared.o $O/write_c_binary.cc.o $O/param_reader.cc.o
	${FCCOMPILE_CHECK} -o  ${E}/xcombine_surf_data  $O/combine_surf_data.shared.o $O/write_c_binary.cc.o $O/param_reader.cc.o

xsmooth_vol_data: $O/smooth_vol_data.o $O/read_parameter_file.shared.o $O/read_value_parameters.shared.o $O/get_value_parameters.shared.o $O/param_reader.cc.o $O/gll_library.shared.o $O/exit_mpi.shared.o $(COND_MPI_OBJECTS)
	${FCLINK} -o  ${E}/xsmooth_vol_data  $O/smooth_vol_data.o $O/read_parameter_file.shared.o $O/read_value_parameters.shared.o $O/get_value_parameters.shared.o $O/param_reader.cc.o $O/gll_library.shared.o $O/exit_mpi.shared.o $(COND_MPI_OBJECTS) $(MPILIBS)

xsum_kernels: $(SUM_KERNELS_OBJECTS) $(COND_MPI_OBJECTS)
	${FCLINK} -o  ${E}/xsum_kernels  $(SUM_KERNELS_OBJECTS) $(COND_MPI_OBJECTS) $(MPILIBS) 

xmodel_update: $(MODEL_UPD_OBJECTS) $(COND_MPI_OBJECTS) $(CUDA_STUBS)
	${FCLINK} -o  ${E}/xmodel_update  $(MODEL_UPD_OBJECTS) $(COND_MPI_OBJECTS) $(MPILIBS) $(CUDA_STUBS)

clean:
	rm -f $O/* *.o *.gnu *.mod $(OUTPUT)/timestamp* $(OUTPUT)/starttime*txt work.pc* \
        xspecfem3D \
        xconvolve_source_timefunction \
        xcreate_movie_shakemap_AVS_DX_GMT xcombine_vol_data xcombine_surf_data \
        xsmooth_vol_data xmodel_update xsum_kernels

###
### rule for the archive library
###

$L/libspecfem.a: $(libspecfem_a_OBJECTS)
	-rm -f $L/libspecfem.a
	$(AR) $(ARFLAGS) $L/libspecfem.a $(libspecfem_a_OBJECTS)
	$(RANLIB) $L/libspecfem.a

####
#### rule to build each .o file below
####

###
### optimized flags (not dependent on values from mesher anymore)
###

$O/%.o: %.f90 $(SHARED)constants.h
	${FCCOMPILE_CHECK} -c -o $@ $<

$O/%.o: %.F90 $(SHARED)constants.h
	${FCCOMPILE_CHECK} -c -o $@ $<

$O/%.shared.o: ${SHARED}%.f90 $(SHARED)constants.h
	${FCCOMPILE_CHECK} -c -o $@ $<

$O/%.shared.o: ${SHARED}%.F90 $(SHARED)constants.h
	${FCCOMPILE_CHECK} -c -o $@ $<


###
### OpenMP compilation
###
$O/%.openmp.o: %.f90 $(SHARED)constants.h
	${FCCOMPILE_CHECK} -c -o $@ $<


###
### CUDA compilation
###
$O/%.cuda.o: ${CUDAD}%.cu ../../config.h $(CUDAD)mesh_constants_cuda.h $(CUDAD)prepare_constants_cuda.h
	$(NVCC) -c $< -o $@ $(NVCC_FLAGS)

###
### C compilation
###
force_ftz.o: ${SHARED}force_ftz.c ../../config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -I../.. -o $O/force_ftz.o ${SHARED}force_ftz.c

$O/%.cc.o: ${SHARED}%.c ../../config.h
	${CC} -c $(CFLAGS) $(MPI_INC) -o $@ ${SHARED}$< -I../../

$O/%.cudacc.o: ${CUDAD}%.c ../../config.h
	${CC} -c $(CFLAGS) $(MPI_INC) -o $@ ${CUDAD}$< -I../../


###
### MPI compilation without optimization
###

$O/parallel.o: $(SHARED)constants.h $(SHARED)parallel.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/parallel.o $(SHARED)parallel.f90

$O/serial.o: $(SHARED)constants.h $(SHARED)serial.f90
	${FCCOMPILE_CHECK} -c -o $O/serial.o $(SHARED)serial.f90

$O/smooth_vol_data.o: $(SHARED)constants.h $(SHARED)smooth_vol_data.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/smooth_vol_data.o $(SHARED)smooth_vol_data.f90




##
## kernel summation
##

$O/sum_kernels.o: $(SHARED)constants.h $(SHARED)sum_kernels.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/sum_kernels.o $(SHARED)sum_kernels.f90

##
## model update
##

$O/model_update.o: $(SHARED)constants.h model_update.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_update.o model_update.f90

$O/save_external_bin_m_up.o:  $(SHARED)constants.h save_external_bin_m_up.f90
	${FCCOMPILE_CHECK} -c -o $O/save_external_bin_m_up.o save_external_bin_m_up.f90

