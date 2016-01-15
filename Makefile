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

# Makefile.  Generated from Makefile.in by configure.

FC = /usr/projects/hpcsoft/toss2.1/common/intel/2011.11.339/composer_xe_2011.11.339/bin/intel64/ifort
FCFLAGS = #-g
MPIFC = mpif90
MPILIBS = 
FLAGS_CHECK = -axAVX,SSE4.2 -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -diag-disable 6477 -implicitnone -warn all -O3 -check nobounds -DFORCE_VECTORIZATION
FCFLAGS_f90 = 

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK)
MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK)
FCLINK = $(MPIFCCOMPILE_CHECK)
#FCLINK = $(FCCOMPILE_CHECK)

CC = /usr/projects/hpcsoft/toss2.1/common/intel/2011.11.339/composer_xe_2011.11.339/bin/intel64/icc
CFLAGS =-axAVX,SSE4.2  $(CPPFLAGS)
CPPFLAGS =  -I/usr/include/scotch $(COND_MPI_CPPFLAGS)
COND_MPI_CPPFLAGS = -DWITH_MPI
#COND_MPI_CPPFLAGS =

AR = ar
ARFLAGS = cru
RANLIB = ranlib

##.PHONY: clean default all backup bak generate_databases specfem3D meshfem3D

####
#### targets
####

# default targets for the pure Fortran version
DEFAULT = \
	xdecompose_mesh \
	meshfem3D \
	generate_databases \
	specfem3D \
	$(EMPTY_MACRO)


default: $(DEFAULT)

all: default \
	xconvolve_source_timefunction \
	xcombine_vol_data \
	xcombine_surf_data \
	xsmooth_vol_data \
	xsum_kernels \
	xmodel_update \
	xcheck_mesh_quality_CUBIT_Abaqus \
	$(EMPTY_MACRO)

required: bin lib obj

backup:
	cp -rp src DATA Makefile go_generate_databases* go_mesher* go_solver* mymachines bak

bak: backup

mesh : meshfem3D
gen : generate_databases
spec : specfem3D
dec : decompose_mesh

generate_databases: xgenerate_databases
specfem3D: xspecfem3D
meshfem3D: xmeshfem3D
decompose_mesh: xdecompose_mesh
convolve_source_timefunction: xconvolve_source_timefunction
create_movie_shakemap_AVS_DX_GMT: xcreate_movie_shakemap_AVS_DX_GMT
combine_vol_data: xcombine_vol_data
combine_surf_data: xcombine_surf_data
smooth_vol_data: xsmooth_vol_data
sum_kernels: xsum_kernels
model_update: xmodel_update
check_mesh_quality_CUBIT_Abaqus: xcheck_mesh_quality_CUBIT_Abaqus

bin:
	mkdir -p bin

lib:
	mkdir -p lib

obj:
	mkdir -p obj

reqmesh:
	(cd obj; mkdir -p mesh)

reqspec:
	(cd obj; mkdir -p spec)

reqdec:
	(cd obj; mkdir -p dec)

reqgen :
	(cd obj; mkdir -p gen)

reqche :
	(cd obj; mkdir -p che)

xmeshfem3D:  required reqmesh
	     (cd src/meshfem3D; make)

xspecfem3D:  required reqspec
	     (cd src/specfem3D; make specfem3D)

xgenerate_databases:  required reqgen
	     (cd src/generate_databases ; make generate_databases)

xdecompose_mesh:  required reqdec
	     (cd src/decompose_mesh ; make)

xcreate_movie_shakemap_AVS_DX_GMT: required
	(cd src/specfem3D ; make xcreate_movie_shakemap_AVS_DX_GMT)

xcombine_vol_data: required reqspec
	(cd src/specfem3D ; make xcombine_vol_data)

xcombine_surf_data: required
	(cd src/specfem3D ; make xcombine_surf_data)

xconvolve_source_timefunction: required
	(cd src/specfem3D ; make xconvolve_source_timefunction)

xsmooth_vol_data: required reqspec
	(cd src/specfem3D ; make xsmooth_vol_data)

xsum_kernels: required reqspec
	(cd src/specfem3D ; make xsum_kernels)

xmodel_update: required reqspec xspecfem3D
	(cd src/specfem3D ; make xmodel_update)


xcheck_mesh_quality_CUBIT_Abaqus: required reqche
	(cd src/check_mesh_quality_CUBIT_Abaqus ; make)


clean: required
	(rm -rf bin lib obj src/meshfem3D/*.mod src/decompose_mesh/*.mod src/generate_databases/*.mod src/specfem3D/*.mod ; cd src/decompose_mesh/scotch/src ; make realclean)


help:
	@echo "usage: make [executable]"
	@echo ""
	@echo "supported executables:"
	@echo "    xgenerate_databases"
	@echo "    xspecfem3D"
	@echo "    xmeshfem3D"
	@echo "    xdecompose_mesh"
	@echo "    xconvolve_source_timefunction"
	@echo "    xcreate_movie_shakemap_AVS_DX_GMT"
	@echo "    xcombine_vol_data"
	@echo "    xcombine_surf_data"
	@echo "    xsmooth_vol_data"
	@echo "    xsum_kernels"
	@echo "    xmodel_update"
	@echo "    xcheck_mesh_quality_CUBIT_Abaqus"
	@echo ""
