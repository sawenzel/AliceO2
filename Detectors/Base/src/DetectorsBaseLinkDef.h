// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class o2::base::Detector + ;
#pragma link C++ class o2::base::Propagator + ;
#pragma link C++ class o2::base::PropagatorF + ;
#pragma link C++ class o2::base::PropagatorD + ;
#pragma link C++ class o2::base::PropagatorImpl < double> + ;
#pragma link C++ class o2::base::PropagatorImpl < float> + ;

#pragma link C++ class o2::base::GeometryManager + ;
#pragma link C++ class o2::base::GeometryManager::MatBudgetExt + ;
#pragma link C++ enum o2::base::MatbudGeomBackend;
#pragma link C++ class o2::base::MaterialManager + ;
#pragma link C++ class o2::MaterialManagerParam + ;
#pragma link C++ class o2::GeometryManagerParam + ;
#pragma link C++ class o2::base::SimFieldUtils + ;

#pragma link C++ class o2::base::Ray + ;
#pragma link C++ class o2::base::MatCell + ;
#pragma link C++ class o2::base::MatBudget + ;
#pragma link C++ class o2::base::MatLayerCyl + ;
#pragma link C++ class o2::base::MatLayerCylSet + ;
#pragma link C++ class o2::base::Aligner + ;
#pragma link C++ class o2::conf::ConfigurableParamHelper < o2::base::Aligner> + ;

#pragma link C++ class o2::GlobalParams + ;
#pragma link C++ class o2::conf::ConfigurableParamHelper < o2::GlobalParams> + ;

#pragma link C++ class o2::data::Stack + ;

#pragma link C++ class o2::base::O2Tessellated - ;
#pragma link C++ class o2::base::BVHSurfaceCurveRecord + ;
#pragma link C++ class o2::base::BVHSurfaceRecord + ;
#pragma link C++ class std::vector < o2::base::BVHSurfaceCurveRecord> + ;
#pragma link C++ class std::vector < o2::base::BVHSurfaceRecord> + ;
#pragma link C++ class o2::base::O2BVHSurfaceSolid - ;
#pragma link C++ class o2::base::O2BVHAssembly + ;
#pragma link C++ class o2::base::FlatCSGHalfspace + ;
#pragma link C++ class o2::base::FlatCSGCell + ;
#pragma link C++ class std::vector < o2::base::FlatCSGHalfspace> + ;
#pragma link C++ class std::vector < o2::base::FlatCSGCell> + ;
#pragma link C++ class o2::base::O2FlatCSG + ;
// O2FlatCSG's sub-cell boxes, its active lists and its BVH are transient by design
// (scripts/geometry/Design_FlatCSGSolid.md section 7: a stored BVH is a second thing that can
// disagree with the data it describes), so a shape read back off a geometry file arrives
// un-closed and answers every query through its _Loop twins -- correct, but with none of the
// acceleration and with nobody told. geom.C calls CloseShape() on what it loads; this makes
// every other reader, o2-sim reading geom.root above all, get the same object.
//
// A read rule rather than a hand-written Streamer. `-` on the link pragma only means "I supply
// my own Streamer", and that Streamer would have to do exactly what the rule below does while
// giving up automatic schema evolution over std::vector<FlatCSGHalfspace> and
// std::vector<FlatCSGCell> -- the wrong trade for a format that is still settling. The refusal
// path is loud rather than silent: CloseShape only refuses a cell bounding box that is missing,
// inverted or non-finite, which is a broken file and not something to carry on from quietly.
#pragma read sourceClass = "o2::base::O2FlatCSG" targetClass = "o2::base::O2FlatCSG" version = "[1-]" source = "" target = "" code = "{ newObj->CloseShape(); if (!newObj->IsClosed()) { newObj->Error(\"Streamer\", \"Shape %s was read from a file and CloseShape() refused it, so it has no sub-cell boxes and every query falls back to its _Loop twin. See the Error above: a cell bounding box is missing, inverted or non-finite.\", newObj->GetName()); } }";

#endif
