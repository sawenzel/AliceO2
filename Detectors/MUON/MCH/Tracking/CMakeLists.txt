# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(MCHTracking
        SOURCES
        src/Cluster.cxx
        src/TrackParam.cxx
        src/Track.cxx
        src/TrackExtrap.cxx
        src/TrackFitter.cxx
        src/TrackFinderOriginal.cxx
        src/TrackFinder.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::Field O2::MCHBase)

o2_add_executable(trackfitter-workflow
        SOURCES
        src/TrackFitterWorkflow.cxx
        src/TrackSamplerSpec.cxx
        src/TrackFitterSpec.cxx
        src/TrackSinkSpec.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES O2::MCHTracking)

o2_add_executable(trackfinderoriginal-workflow
        SOURCES
        src/TrackFinderOriginalWorkflow.cxx
        src/ClusterSamplerSpec.cxx
        src/TrackFinderOriginalSpec.cxx
        src/TrackSinkSpec.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES O2::MCHTracking)

        o2_add_executable(trackfinder-workflow
        SOURCES
        src/TrackFinderWorkflow.cxx
        src/ClusterSamplerSpec.cxx
        src/VertexSamplerSpec.cxx
        src/TrackFinderSpec.cxx
        src/TrackAtVertexSpec.cxx
        src/TrackSinkSpec.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2::MCHTracking)
