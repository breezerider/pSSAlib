
LINKLIBS = $(OPT_LIBS) $(GSL_LIBS) $(SBML_LDFLAGS)

AUTOMAKE_OPTIONS = subdir-objects

lib_LTLIBRARIES = libpssa.la
libpssa_la_SOURCES = \
PSSA.cpp \
datamodel/DataModels.cpp \
datamodel/PSSACR_Bins.cpp \
datamodel/SimulationInfo.cpp \
datamodel/detail/Model.cpp \
datamodel/detail/Reaction.cpp \
datamodel/detail/Species.cpp \
datamodel/detail/SpeciesReference.cpp \
datamodel/detail/XMLTypeDefinitions.cpp \
grouping/GroupingModule.cpp \
grouping/GroupingModule_DM.cpp \
grouping/GroupingModule_PDM.cpp \
grouping/GroupingModule_PSSACR.cpp \
grouping/GroupingModule_SPDM.cpp \
sampling/CompositionRejectionSampler.cpp \
sampling/inc/SamplingModule_S_PDM.inc \
sampling/SamplingModule.cpp \
sampling/SamplingModule_DM.cpp \
sampling/SamplingModule_PDM.cpp \
sampling/SamplingModule_PSSACR.cpp \
sampling/SamplingModule_SPDM.cpp \
update/UpdateModule.cpp \
update/UpdateModule_DM.cpp \
update/UpdateModule_PDM.cpp \
update/UpdateModule_PSSACR.cpp \
update/UpdateModule_SPDM.cpp \
util/MPIWrapper.cpp \
util/FileSystem.cpp

libpssa_la_CFLAGS = -DUNIX -rdynamic $(GSL_CFLAGS) $(SBML_CPPFLAGS)
libpssa_la_CXXFLAGS = -DUNIX -rdynamic $(GSL_CFLAGS) $(SBML_CPPFLAGS)
libpssa_la_LIBADD = $(LINKLIBS)
