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

/// \file MatLayerCylSet.cxx
/// \brief Implementation of the wrapper for the set of cylindrical material layers

#include "DetectorsBase/MatLayerCylSet.h"
#include "CommonConstants/MathConstants.h"
#include "CommonUtils/TreeStream.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "png.h"
#include <gsl/span>
#include "TStopwatch.h"

#ifndef GPUCA_ALIGPUCODE // this part is unvisible on GPU version

#include "GPUCommonLogger.h"
#include <TFile.h>
#include "CommonUtils/TreeStreamRedirector.h"
//#define _DBG_LOC_ // for local debugging only

#endif // !GPUCA_ALIGPUCODE
#undef NDEBUG
using namespace o2::base;

using flatObject = o2::gpu::FlatObject;

#ifndef GPUCA_ALIGPUCODE // this part is unvisible on GPU version

//________________________________________________________________________________
void MatLayerCylSet::addLayer(float rmin, float rmax, float zmax, float dz, float drphi)
{
  // add new layer checking for overlaps
  assert(mConstructionMask != Constructed);
  assert(rmin < rmax && zmax > 0 && dz > 0 && drphi > 0);
  mConstructionMask = InProgress;
  int nlr = getNLayers();
  if (!nlr) {
    // book local storage
    auto sz = sizeof(MatLayerCylSetLayout);
    o2::gpu::resizeArray(mFlatBufferContainer, 0, sz);
    mFlatBufferPtr = mFlatBufferContainer;
    mFlatBufferSize = sz;
    //--------------????
    get()->mRMin = 1.e99;
    get()->mRMax = 0.;
  }

  for (int il = 0; il < nlr; il++) {
    const auto& lr = getLayer(il);
    if (lr.getRMax() > rmin && rmax > lr.getRMin()) {
      LOG(fatal) << "new layer overlaps with layer " << il;
    }
  }
  auto* oldLayers = o2::gpu::resizeArray(get()->mLayers, nlr, nlr + 1);
  // dynamyc buffers of old layers were used in new ones, detach them
  for (int i = nlr; i--;) {
    oldLayers[i].clearInternalBufferPtr();
  }
  delete[] oldLayers;
  get()->mLayers[nlr].initSegmentation(rmin, rmax, zmax, dz, drphi);
  get()->mNLayers++;
  get()->mRMin = get()->mRMin > rmin ? rmin : get()->mRMin;
  get()->mRMax = get()->mRMax < rmax ? rmax : get()->mRMax;
  get()->mZMax = get()->mZMax < zmax ? zmax : get()->mZMax;
  get()->mRMin2 = get()->mRMin * get()->mRMin;
  get()->mRMax2 = get()->mRMax * get()->mRMax;
}

//________________________________________________________________________________
void MatLayerCylSet::populateFromTGeo(int ntrPerCell)
{
  ///< populate layers, using ntrPerCell test tracks per cell
  assert(mConstructionMask == InProgress);

  int nlr = getNLayers();
  if (!nlr) {
    LOG(error) << "The LUT is not yet initialized";
    return;
  }
  if (get()->mR2Intervals) {
    LOG(error) << "The LUT is already populated";
    return;
  }
  for (int i = 0; i < nlr; i++) {
    printf("Populating with %d trials Lr  %3d ", ntrPerCell, i);
    get()->mLayers[i].print();
    get()->mLayers[i].populateFromTGeo(ntrPerCell);
  }
  finalizeStructures();
}

//________________________________________________________________________________
void MatLayerCylSet::finalizeStructures()
{
  // build layer search structures
  assert(mConstructionMask == InProgress);
  int nlr = getNLayers();
  int nR2Int = 2 * (nlr + 1);
  o2::gpu::resizeArray(get()->mR2Intervals, 0, nR2Int);
  o2::gpu::resizeArray(get()->mInterval2LrID, 0, nR2Int);
  get()->mR2Intervals[0] = get()->mRMin2;
  get()->mR2Intervals[1] = get()->mRMax2;
  get()->mInterval2LrID[0] = 0;
  auto& nRIntervals = get()->mNRIntervals;
  nRIntervals = 1;

  for (int i = 1; i < nlr; i++) {
    const auto& lr = getLayer(i);
    if (std::sqrt(lr.getRMin2()) > std::sqrt(get()->mR2Intervals[nRIntervals] + Ray::Tiny)) {
      // register gap
      get()->mInterval2LrID[nRIntervals] = -1;
      get()->mR2Intervals[++nRIntervals] = lr.getRMin2();
    }
    get()->mInterval2LrID[nRIntervals] = i;
    get()->mR2Intervals[++nRIntervals] = lr.getRMax2();
  }
  delete[] o2::gpu::resizeArray(get()->mInterval2LrID, nR2Int, nRIntervals); // rebook with precise size
  delete[] o2::gpu::resizeArray(get()->mR2Intervals, nR2Int, ++nRIntervals); // rebook with precise size
  //
}

//________________________________________________________________________________
void MatLayerCylSet::dumpToTree(const std::string& outName) const
{
  /// dump per cell info to the tree

  o2::utils::TreeStreamRedirector dump(outName.data(), "recreate");
  for (int i = 0; i < getNLayers(); i++) {
    const auto& lr = getLayer(i);
    float r = 0.5 * (lr.getRMin() + lr.getRMax());
    // per cell dump
    int nphib = lr.getNPhiBins();
    for (int ip = 0; ip < nphib; ip++) {
      float phi = 0.5 * (lr.getPhiBinMin(ip) + lr.getPhiBinMax(ip));
      float sn, cs;
      int ips = lr.phiBin2Slice(ip);
      char merge = 0; // not mergeable
      if (ip + 1 < nphib) {
        int ips1 = lr.phiBin2Slice(ip + 1);
        merge = ips == ips1 ? -1 : lr.canMergePhiSlices(ips, ips1); // -1 for already merged
      } else {
        merge = -2; // last one
      }
      o2::math_utils::sincos(phi, sn, cs);
      float x = r * cs, y = r * sn;
      for (int iz = 0; iz < lr.getNZBins(); iz++) {
        float z = 0.5 * (lr.getZBinMin(iz) + lr.getZBinMax(iz));
        auto cell = lr.getCellPhiBin(ip, iz);
        dump << "cell"
             << "ilr=" << i << "r=" << r << "phi=" << phi << "x=" << x << "y=" << y << "z=" << z << "ip=" << ip << "ips=" << ips << "iz=" << iz
             << "mrgnxt=" << merge << "val=" << cell << "\n";
      }
    }
    //
    // statistics per layer
    MatCell mean, rms;
    lr.getMeanRMS(mean, rms);
    dump << "lay"
         << "ilr=" << i << "r=" << r << "mean=" << mean << "rms=" << rms << "\n";
  }
}

//________________________________________________________________________________
void MatLayerCylSet::writeToFile(const std::string& outFName)
{
  /// store to file

  TFile outf(outFName.data(), "recreate");
  if (outf.IsZombie()) {
    return;
  }
  outf.WriteObjectAny(this, Class(), "ccdb_object");
  outf.Close();
}

GPUd() int MatLayerCylSet::searchLayerFast(float r2, int low, int high) const {
  // we can avoid the sqrt .. at the cost of more memory in the lookup
  const auto index = int(std::sqrt(r2) * InvVoxelRDelta);
  if (index >= mLayerVoxelLU.size()) {
    return -1;
  }
  auto layers = mLayerVoxelLU[index];
  if (layers.first != layers.second) {
    // this means the voxel is undecided and we revert to search
    // LOG(info) << " Searching  " << layers.first << " " << layers.second + 1; 
    return searchSegment(r2, layers.first, layers.second + 1);
  }
  // LOG(info) << " Return unique " << layers.first;
  return layers.first;
}

void MatLayerCylSet::initLayerVoxelLU() const {
  if (mLayerVoxelLU.size() > 0) {
    LOG(info) << "Layer voxel already initialized; Aborting";
    return;
  }
  LOG(info) << "Initializing Voxel layer lookup";
  const int numVoxels = int(-XMIN / VoxelRDelta);
  mLayerVoxelLU.clear();
  mLayerVoxelLU.reserve(numVoxels);
  for (int voxel = 0; voxel < numVoxels; ++voxel) {
    // check the 2 extremes of this voxel "covering"
    const auto lowerR = voxel*VoxelRDelta;
    const auto upperR = lowerR + VoxelRDelta;
    const auto lowerSegment = searchSegment(lowerR*lowerR);
    const auto upperSegment = searchSegment(upperR*upperR);
    mLayerVoxelLU.push_back(std::pair<uint16_t, u_int16_t>(lowerSegment, upperSegment));
  }
  LOG(info) << "LUT has size " << mLayerVoxelLU.size();
  mVoxelInitialized = true;
}


// writes PNG image representation of voxel data
void writePNGImage(const char* pngname, std::vector<std::vector<uint8_t>> const& imgdata) {
  // open the file
  FILE *fp = fopen(pngname, "wb");
  if (!fp) {
    // std::cerr << "Error: Could not open output.png for writing.\n";
    return;
  }
  // Initialize libpng
  png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
  if (!png) {
    // std::cerr << "Error: Could not initialize libpng.\n";
    fclose(fp);
    return;
  }

  png_infop info = png_create_info_struct(png);
  if (!info) {
    // std::cerr << "Error: Could not create PNG info struct\n";
    png_destroy_write_struct(&png, nullptr);
    fclose(fp);
    return;
  }
  png_init_io(png, fp);
  // Set image properties
  png_set_IHDR(png, info, imgdata.size(), imgdata.size(), 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  // Write image data row by row
  png_write_info(png, info);
  for (int y = 0; y < imgdata.size(); ++y) {
    png_write_row(png, &imgdata[y][0]);
  }
  // Finalize and close the PNG file
  png_write_end(png, nullptr);
  png_destroy_write_struct(&png, &info);
  fclose(fp);

}

void MatLayerCylSet::initPhiSectorVoxelLU(const char* pngname) const {
    if (mPhiSectorVoxelLU) {
      LOG(warn) << "PhiSector lookup already initialized; Aborting";
      return;
    }
    
    mPhiSectorVoxelLU = new std::vector<std::vector<uint16_t>>;
    mPhiSectorVoxelLU->resize(NVoxels1D); // size the Y-dimensions
    for (int row = 0; row < NVoxels1D; ++row) {
      (*mPhiSectorVoxelLU)[row].resize(NVoxels1D); // size the X-dimensions
    }

    std::vector<std::vector<uint8_t>> pngData;
    if (pngname) {
      pngData.resize(NVoxels1D); // size the Y-dimensions
      for (int row = 0; row < NVoxels1D; ++row) {
        pngData[row].resize(NVoxels1D); // size the X-dimensions
      }
    }

    // resize the VecGeom voxel hashmap
    mPhiSectorLUTs.resize(1, nullptr);
    mPhiSectorLUTs[0] = new HashMap2D(vecgeom::Vector2D<float>(XMIN,XMIN), vecgeom::Vector2D<float>(VDELTA, VDELTA), NVoxels1D, NVoxels1D);

    size_t longest  = 0;
    int longestindex = 0;
    std::array<int, 5> counter{0,0,0,0};

    auto PropList = [this](auto key) {
      int s = 0;
      const auto ptr = mPhiSectorLUTs[0]->getPropertiesGivenKey(key, s);
      return gsl::span<const std::pair<short, short>>(ptr, s);
    };

    // we loop over all voxel indixes
    // calculate the 4 corner coordinates of the voxel and check each of them
    // if result is unique we'll store it --> if not MAX_VALUE OF SHORT
    for (int yi = 0; yi < NVoxels1D; ++yi) { // iterate over rows
      for (int xi = 0; xi < NVoxels1D; ++xi) { // iterate within row (cash friendly)
        // LOG(info) << xi << " " << yi;
        auto points = getVoxelPoints(xi, yi);

        const auto voxelkey = mPhiSectorLUTs[0]->getKeyFromCells(xi, yi);

        const auto INVALID = std::numeric_limits<uint16_t>::max(); // IN PNG this is white
        int phiIndex = -1; // unitialized
        for (auto& p : points) {
          const auto r2 = p.first*p.first + p.second*p.second;
          if (r2 >= XMIN*XMIN) {
            break;
          }
          
          const auto interval = searchLayerFast(r2); // this is the interval - not the layer
          if (interval >= 0) {
            // getPhiID
            // we get the phi angle of that point
            float phi = o2::gpu::CAMath::ATan2(p.second, p.first);
            o2::math_utils::bringTo02Pi(phi);

            // sector
            // we first need the layer structure 
            const auto& lr = getLayer(interval);
            const auto phisector = lr.getPhiSliceID(phi);
            std::pair<short, short> p(interval, phisector);
            bool exists = false;
            for (const auto& existing : PropList(voxelkey)) {
              if (existing == p) {
                exists = true;
                break;
              }
            }
            if (!exists) {
              mPhiSectorLUTs[0]->addPropertyForKey(voxelkey, p);
            }
            if (phiIndex == -1) {
              phiIndex = phisector;
            }
            else {
              if (phisector != phiIndex) {
                phiIndex = INVALID;
              }
            }
          }
        }
        if (phiIndex == -1) {
          phiIndex = INVALID;
        }
        // the first index is the "row"
        const uint16_t final = (uint16_t)phiIndex;
        (*mPhiSectorVoxelLU)[yi][xi] = final;

        // do a check here
        auto props = PropList(voxelkey);
        if (props.size() > longest) {
          longestindex = voxelkey;
          longest = props.size();
        }
        longest = std::max(props.size(), longest);
        if (props.size() == 0 && (*mPhiSectorVoxelLU)[yi][xi]!=INVALID) {
          std::cerr << "Inconsistency type 1\n";
        }
        if (props.size() > 1 && (*mPhiSectorVoxelLU)[yi][xi]!=INVALID) {
          // here all phi sectors need to be the same
          short phi = -1;
          bool consistent = true;
          for (auto &p : PropList(voxelkey)) {
            if (phi == -1) {
              phi = p.second;
            }
            if (p.second != phi) {
              consistent = false;
            }
            if (!consistent) {
              std::cerr << "Inconsistency type 2\n";
            }
          }
        }
        if (props.size() == 1 && (*mPhiSectorVoxelLU)[yi][xi]==INVALID) {
          std::cerr << "Inconsistency type 3\n";
        }
        counter[PropList(voxelkey).size()]++;
        if (pngname) {
          pngData[yi][xi] = (final == INVALID) ? std::numeric_limits<uint8_t>::max() : final % 128;
        }
      }
    }
    std::cout << "Largest voxel " << longest << "\n";
    for (auto &p : PropList(longestindex)) {
      std::cout << "( " << p.first << " , " << p.second << " ) ";
    }
    std::cout << "\n";

    for (int i = 0; i< 5; ++i) {
      std::cout << "bin " << i << " : " << counter[i] << "\n";
    }
    if (pngname) {
      LOG(info) << "Writing voxels to PNG image";
      // trim to 8 bit gray scale values first
      writePNGImage(pngname, pngData);
    }
}

//________________________________________________________________________________
MatLayerCylSet* MatLayerCylSet::loadFromFile(const std::string& inpFName)
{
  TFile inpf(inpFName.data());
  if (inpf.IsZombie()) {
    LOG(error) << "Failed to open input file " << inpFName;
    return nullptr;
  }
  MatLayerCylSet* mb = reinterpret_cast<MatLayerCylSet*>(inpf.GetObjectChecked("ccdb_object", Class()));
  if (!mb && !(mb = reinterpret_cast<MatLayerCylSet*>(inpf.GetObjectChecked("MatBud", Class())))) { // for old objects
    LOG(error) << "Failed to load mat.LUT from " << inpFName;
    return nullptr;
  }
  auto rptr = rectifyPtrFromFile(mb);
  return rptr;
}

//________________________________________________________________________________
MatLayerCylSet* MatLayerCylSet::rectifyPtrFromFile(MatLayerCylSet* ptr)
{
  // rectify object loaded from file
  if (ptr && !ptr->get()) {
    ptr->fixPointers();
  }
  return ptr;
}

//________________________________________________________________________________
void MatLayerCylSet::optimizePhiSlices(float maxRelDiff)
{
  // merge similar (whose relative budget does not differ within maxRelDiff) phi slices
  assert(mConstructionMask == InProgress);
  for (int i = getNLayers(); i--;) {
    get()->mLayers[i].optimizePhiSlices(maxRelDiff);
  }
  // flatten();  // RS: TODO
}

//________________________________________________________________________________
void MatLayerCylSet::print(bool data) const
{
  ///< print layer data
  if (!get()) {
    printf("Not initialized yet\n");
    return;
  }
  if (mConstructionMask != Constructed) {
    LOG(warning) << "Object is not yet flattened";
  }
  for (int i = 0; i < getNLayers(); i++) {
    printf("#%3d | ", i);
    getLayer(i).print(data);
  }
  printf("%.2f < R < %.2f  %d layers with total size %.2f MB\n", getRMin(), getRMax(), getNLayers(),
         float(getFlatBufferSize()) / 1024 / 1024);
}

#endif //!GPUCA_ALIGPUCODE

#ifndef GPUCA_GPUCODE
//________________________________________________________________________________
std::size_t MatLayerCylSet::estimateFlatBufferSize() const
{
  std::size_t sz = alignSize(sizeof(MatLayerCylSetLayout), getBufferAlignmentBytes()); // hold data members

  sz = alignSize(sz + get()->mNLayers * sizeof(MatLayerCyl), MatLayerCyl::getClassAlignmentBytes());
  sz = alignSize(sz + (get()->mNRIntervals + 1) * sizeof(float), getBufferAlignmentBytes());
  sz = alignSize(sz + get()->mNRIntervals * sizeof(int), getBufferAlignmentBytes());

  for (int i = 0; i < getNLayers(); i++) {
    sz = alignSize(sz + getLayer(i).estimateFlatBufferSize(), getBufferAlignmentBytes());
  }
  return sz;
}
#endif // ! GPUCA_GPUCODE

//_________________________________________________________________________________________________
GPUd() MatBudget MatLayerCylSet::getMatBudget(float x0, float y0, float z0, float x1, float y1, float z1) const
{
  static int goodcounter = 0;
  static int badcounter = 0;
  static int cmpgood = 0;
  static int cmpbad = 0;

  // get material budget traversed on the line between point0 and point1
  MatBudget rval;
  Ray ray(x0, y0, z0, x1, y1, z1);
  short lmin, lmax; // get innermost and outermost relevant layer
  if (ray.isTooShort() || !getLayersRange(ray, lmin, lmax)) {
    rval.length = ray.getDist();
    return rval;
  }
  short lrID = lmax;
  while (lrID >= lmin) { // go from outside to inside
    const auto& lr = getLayer(lrID);
    int nphiSlices = lr.getNPhiSlices();
    int nc = ray.crossLayer(lr);  // determines how many crossings this ray has with this tubular layer 
    for (int ic = nc; ic--;) {
      float cross1, cross2;
      ray.getCrossParams(ic, cross1, cross2); // tmax,tmin of crossing the layer
      
      auto getPhiSlicesExact = [nphiSlices, &ray, &lr](auto cr1, auto cr2, int& phiID, int& phiIDLast) -> void {
        auto phi0 = ray.getPhi(cr1), phi1 = ray.getPhi(cr2), dPhi = phi0 - phi1;
        phiID = lr.getPhiSliceID(phi0);
        phiIDLast = lr.getPhiSliceID(phi1);
        // account for eventual wrapping around 0
        if (dPhi > 0.f) {
          if (dPhi > o2::constants::math::PI) { // wraps around phi=0
            phiIDLast += nphiSlices;
          }
        } else {
          if (dPhi < -o2::constants::math::PI) { // wraps around phi=0
            phiID += nphiSlices;
          }
        }
      };
      
      auto lookupOne = [&ray, &lr, this](auto t, float& phi, bool& slow) -> int {
        const auto px = ray.getPos(t, 0);
        const auto py = ray.getPos(t, 1);
        int xi, yi; // the indices into the voxel structure
        findPhiVoxel(px, py, 0, xi, yi);
        // use the lookup
        const auto phiSliceTmp = (*mPhiSectorVoxelLU)[yi][xi];
        // is it valid?
        if (phiSliceTmp != std::numeric_limits<uint16_t>::max()) {
          // yes
          goodcounter++;
          return phiSliceTmp;
        }
        badcounter++;
        // go the slow path
        std::cout << "BAD: ";
        for (const auto& p : getPhiPropList(xi, yi)) {
          std::cout << "( " << p.first << " " << p.second << ")";
        }
        std::cout << "\n";
        // ...
        slow = true;
        phi = ray.getPhi(t);
        return lr.getPhiSliceID(phi);
      };
      
      auto getPhiSlicesFast = [nphiSlices, &ray, &lr, &lookupOne](auto cr1, auto cr2, int& phiID, int& phiIDLast) {
        bool slow1 = false;
        float phi1;
        phiID = lookupOne(cr1, phi1, slow1);
      
        bool slow2 = false;
        float phi2;
        phiIDLast = lookupOne(cr2, phi2, slow2);
      };
    
      int phiID, phiIDLast;
      getPhiSlicesExact(cross1, cross2, phiID, phiIDLast);
      
      int phiIDcmp, phiIDLastcmp;
      getPhiSlicesFast(cross1, cross2, phiIDcmp, phiIDLastcmp);
    
      if ((phiID != phiIDcmp) || ( phiIDLast != phiIDLastcmp)) {
        cmpbad++;
        LOG(info) << phiID << " " << phiIDcmp << " : " << phiIDLast << " " << phiIDLastcmp << " goodsofar " << goodcounter << " badsofar " << badcounter << " " << cmpgood << " " << cmpbad;
      }
      else {
        cmpgood++;
      }

      /*
      auto phi0 = ray.getPhi(cross1), phi1 = ray.getPhi(cross2), dPhi = phi0 - phi1;
      auto phiID = lr.getPhiSliceID(phi0), phiIDLast = lr.getPhiSliceID(phi1);
      // account for eventual wrapping around 0
      if (dPhi > 0.f) {
        if (dPhi > o2::constants::math::PI) { // wraps around phi=0
          phiIDLast += nphiSlices;
        }
      } else {
        if (dPhi < -o2::constants::math::PI) { // wraps around phi=0
          phiID += nphiSlices;
        }
      }
      */
       
      
      int stepPhiID = phiID > phiIDLast ? -1 : 1;
      bool checkMorePhi = true;
      auto tStartPhi = cross1, tEndPhi = 0.f;
      do {
        // get the path in the current phi slice
        if (phiID == phiIDLast) {
          tEndPhi = cross2;
          checkMorePhi = false;
        } else { // last phi slice still not reached
          tEndPhi = ray.crossRadial(lr, (stepPhiID > 0 ? phiID + 1 : phiID) % nphiSlices);
          if (tEndPhi == Ray::InvalidT) {
            break; // ray parallel to radial line, abandon check for phi bin change
          }
        }
        auto zID = lr.getZBinID(ray.getZ(tStartPhi));
        auto zIDLast = lr.getZBinID(ray.getZ(tEndPhi));
        // check if Zbins are crossed

#ifdef _DBG_LOC_
        printf("-- Zdiff (%3d : %3d) mode: t: %+e %+e\n", zID, zIDLast, tStartPhi, tEndPhi);
#endif

        if (zID != zIDLast) {
          auto stepZID = zID < zIDLast ? 1 : -1;
          bool checkMoreZ = true;
          auto tStartZ = tStartPhi, tEndZ = 0.f;
          do {
            if (zID == zIDLast) {
              tEndZ = tEndPhi;
              checkMoreZ = false;
            } else {
              tEndZ = ray.crossZ(lr.getZBinMin(stepZID > 0 ? zID + 1 : zID));
              if (tEndZ == Ray::InvalidT) { // track normal to Z axis, abandon Zbin change test
                break;
              }
            }
            // account materials of this step
            float step = tEndZ > tStartZ ? tEndZ - tStartZ : tStartZ - tEndZ; // the real step is ray.getDist(tEnd-tStart), will rescale all later
            const auto& cell = lr.getCell(phiID % nphiSlices, zID);
            rval.meanRho += cell.meanRho * step;
            rval.meanX2X0 += cell.meanX2X0 * step;
            rval.length += step;

#ifdef _DBG_LOC_
            float pos0[3] = {ray.getPos(tStartZ, 0), ray.getPos(tStartZ, 1), ray.getPos(tStartZ, 2)};
            float pos1[3] = {ray.getPos(tEndZ, 0), ray.getPos(tEndZ, 1), ray.getPos(tEndZ, 2)};
            printf(
              "Lr#%3d / cross#%d : account %f<t<%f at phiSlice %d | Zbin: %3d (%3d) |[%+e %+e +%e]:[%+e %+e %+e] "
              "Step: %.3e StrpCor: %.3e\n",
              lrID, ic, tEndZ, tStartZ, phiID % nphiSlices, zID, zIDLast,
              pos0[0], pos0[1], pos0[2], pos1[0], pos1[1], pos1[2], step, ray.getDist(step));
#endif

            tStartZ = tEndZ;
            zID += stepZID;
          } while (checkMoreZ);
        } else {
          float step = tEndPhi > tStartPhi ? tEndPhi - tStartPhi : tStartPhi - tEndPhi; // the real step is |ray.getDist(tEnd-tStart)|, will rescale all later
          const auto& cell = lr.getCell(phiID % nphiSlices, zID);
          rval.meanRho += cell.meanRho * step;
          rval.meanX2X0 += cell.meanX2X0 * step;
          rval.length += step;

#ifdef _DBG_LOC_
          float pos0[3] = {ray.getPos(tStartPhi, 0), ray.getPos(tStartPhi, 1), ray.getPos(tStartPhi, 2)};
          float pos1[3] = {ray.getPos(tEndPhi, 0), ray.getPos(tEndPhi, 1), ray.getPos(tEndPhi, 2)};
          printf(
            "Lr#%3d / cross#%d : account %f<t<%f at phiSlice %d | Zbin: %3d ----- |[%+e %+e +%e]:[%+e %+e %+e]"
            "Step: %.3e StrpCor: %.3e\n",
            lrID, ic, tEndPhi, tStartPhi, phiID % nphiSlices, zID,
            pos0[0], pos0[1], pos0[2], pos1[0], pos1[1], pos1[2], step, ray.getDist(step));
#endif
        }
        //
        tStartPhi = tEndPhi;
        phiID += stepPhiID;

      } while (checkMorePhi);
    }
    lrID--;
  } // loop over layers

  if (rval.length != 0.f) {
    rval.meanRho /= rval.length;                                       // average
    rval.meanX2X0 *= ray.getDist();                                    // normalize
  }
  rval.length = ray.getDist();

#ifdef _DBG_LOC_
  printf("<rho> = %e, x2X0 = %e  | step = %e\n", rval.meanRho, rval.meanX2X0, rval.length);
#endif
  return rval;
}

//_________________________________________________________________________________________________
GPUd() bool MatLayerCylSet::getLayersRange(const Ray& ray, short& lmin, short& lmax) const
{
  // get range of layers corresponding to rmin/rmax
  //
  lmin = lmax = -1;
  float rmin2, rmax2;
  ray.getMinMaxR2(rmin2, rmax2);

  if (rmin2 >= getRMax2() || rmax2 <= getRMin2()) {
    return false;
  }
  int lmxInt, lmnInt;
  if (!mVoxelInitialized) {
    lmxInt = rmax2 < getRMax2() ? searchSegment(rmax2, 0) : get()->mNRIntervals - 2;
    lmnInt = rmin2 >= getRMin2() ? searchSegment(rmin2, 0, lmxInt + 1) : 0;
  }
  else {
    lmxInt = rmax2 < getRMax2() ? searchLayerFast(rmax2, 0) : get()->mNRIntervals - 2;
    lmnInt = rmin2 >= getRMin2() ? searchLayerFast(rmin2, 0, lmxInt + 1) : 0;
  }

  // cross check
  /*
  if (mVoxelInitialized) {
    int lmxIntCheck = rmax2 < getRMax2() ? searchLayerFast(rmax2, 0) : get()->mNRIntervals - 2;
    if (lmxInt != lmxIntCheck) {
      LOG(warn) << " Difference mx " << lmxInt << " " << lmxIntCheck;
    }
    int lmnIntCheck = rmin2 >= getRMin2() ? searchLayerFast(rmin2, 0, lmxIntCheck + 1) : 0;
    if (lmnInt != lmnIntCheck) {
      LOG(warn) << " Difference mn " << lmnInt 
                << " " << lmnIntCheck << " rmin2 " << rmin2 << " " << lmxInt << " " 
                << searchSegment(rmin2) << " " << getRMin2() << " " << searchSegment(rmin2, 0, lmxIntCheck + 1);
    }
  }
  */
  
  const auto* interval2LrID = get()->mInterval2LrID;
  lmax = interval2LrID[lmxInt];
  lmin = interval2LrID[lmnInt];
  // make sure lmnInt and/or lmxInt are not in the gap
  if (lmax < 0) {
    lmax = interval2LrID[lmxInt - 1]; // rmax2 is in the gap, take highest layer below rmax2
  }
  if (lmin < 0) {
    lmin = interval2LrID[lmnInt + 1]; // rmin2 is in the gap, take lowest layer above rmin2
  }
  return lmin <= lmax; // valid if both are not in the same gap
}

GPUd() int MatLayerCylSet::searchSegment(float val, int low, int high) const
{
  ///< search segment val belongs to. The val MUST be within the boundaries
  if (low < 0) {
    low = 0;
  }
  if (high < 0) {
    high = get()->mNRIntervals;
  }
  int mid = (low + high) >> 1;
  const auto* r2Intervals = get()->mR2Intervals;
  while (mid != low) {
    if (val < r2Intervals[mid]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) >> 1;
  }

  return mid;
}

#ifndef GPUCA_ALIGPUCODE // this part is unvisible on GPU version

void MatLayerCylSet::flatten()
{
  // make object flat: move all content to single internally allocated buffer
  assert(mConstructionMask == InProgress);

  int sz = estimateFlatBufferSize();
  // create new internal buffer with total size and copy data
  delete[] o2::gpu::resizeArray(mFlatBufferContainer, mFlatBufferSize, sz);
  mFlatBufferPtr = mFlatBufferContainer;
  mFlatBufferSize = sz;
  int nLr = getNLayers();

  auto offs = alignSize(sizeof(MatLayerCylSetLayout), getBufferAlignmentBytes()); // account for the alignment
  // move array of layer pointers to the flat array
  auto* oldLayers = o2::gpu::resizeArray(get()->mLayers, nLr, nLr, (MatLayerCyl*)(mFlatBufferPtr + offs));
  // dynamyc buffers of old layers were used in new ones, detach them
  for (int i = nLr; i--;) {
    oldLayers[i].clearInternalBufferPtr();
  }
  delete[] oldLayers;
  offs = alignSize(offs + nLr * sizeof(MatLayerCyl), MatLayerCyl::getClassAlignmentBytes()); // account for the alignment

  // move array of R2 boundaries to the flat array
  delete[] o2::gpu::resizeArray(get()->mR2Intervals, nLr + 1, nLr + 1, (float*)(mFlatBufferPtr + offs));
  offs = alignSize(offs + (nLr + 1) * sizeof(float), getBufferAlignmentBytes()); // account for the alignment

  // move array of R2 boundaries to the flat array
  delete[] o2::gpu::resizeArray(get()->mInterval2LrID, nLr, nLr, (int*)(mFlatBufferPtr + offs));
  offs = alignSize(offs + nLr * sizeof(int), getBufferAlignmentBytes()); // account for the alignment

  for (int il = 0; il < nLr; il++) {
    MatLayerCyl& lr = get()->mLayers[il];
    lr.flatten(mFlatBufferPtr + offs);
    offs = alignSize(offs + lr.getFlatBufferSize(), getBufferAlignmentBytes()); // account for the alignment
  }
  mConstructionMask = Constructed;
}

//______________________________________________
void MatLayerCylSet::moveBufferTo(char* newFlatBufferPtr)
{
  /// sets buffer pointer to the new address, move the buffer content there.
  flatObject::moveBufferTo(newFlatBufferPtr);
  setActualBufferAddress(mFlatBufferPtr);
}
#endif // !GPUCA_ALIGPUCODE

#ifndef GPUCA_GPUCODE
//______________________________________________
void MatLayerCylSet::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// Sets the actual location of the external flat buffer before it was created
  ///
  fixPointers(mFlatBufferPtr, futureFlatBufferPtr, false); // flag that futureFlatBufferPtr is not valid yet
  flatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

//______________________________________________
void MatLayerCylSet::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// Sets the actual location of the external flat buffer after it has been moved (i.e. to another machine)
  ///
  fixPointers(actualFlatBufferPtr);
}
//______________________________________________
void MatLayerCylSet::cloneFromObject(const MatLayerCylSet& obj, char* newFlatBufferPtr)
{
  /// Initializes from another object, copies data to newBufferPtr
  flatObject::cloneFromObject(obj, newFlatBufferPtr);
  fixPointers(mFlatBufferPtr);
}

//______________________________________________
void MatLayerCylSet::fixPointers(char* newBasePtr)
{
  // fix pointers on the internal structure of the flat buffer after retrieving it from the file
  if (newBasePtr) {
    mFlatBufferPtr = newBasePtr; // used to impose external pointer
  } else {
    mFlatBufferPtr = mFlatBufferContainer; // impose pointer after reading from file
  }
  auto offs = alignSize(sizeof(MatLayerCylSetLayout), getBufferAlignmentBytes()); // account for the alignment
  char* newPtr = mFlatBufferPtr + offs;                                           // correct pointer on MatLayerCyl*
  char* oldPtr = reinterpret_cast<char*>(get()->mLayers);                         // old pointer read from the file
  fixPointers(oldPtr, newPtr);
}

//______________________________________________
void MatLayerCylSet::fixPointers(char* oldPtr, char* newPtr, bool newPtrValid)
{
  // fix pointers on the internal structure of the flat buffer after retrieving it from the file
  auto* layPtr = get()->mLayers;
  get()->mLayers = flatObject::relocatePointer(oldPtr, newPtr, get()->mLayers);
  get()->mR2Intervals = flatObject::relocatePointer(oldPtr, newPtr, get()->mR2Intervals);
  get()->mInterval2LrID = flatObject::relocatePointer(oldPtr, newPtr, get()->mInterval2LrID);
  if (newPtrValid) {
    layPtr = get()->mLayers;
  }
  for (int i = 0; i < getNLayers(); i++) {
    layPtr[i].setFlatPointer(flatObject::relocatePointer(oldPtr, newPtr, layPtr[i].getFlatBufferPtr()));
    layPtr[i].fixPointers(oldPtr, newPtr);
  }
}
#endif // !GPUCA_GPUCODE

#ifndef GPUCA_ALIGPUCODE // this part is unvisible on GPU version

MatLayerCylSet* MatLayerCylSet::extractCopy(float rmin, float rmax, float tolerance) const
{
  Ray ray(std::max(getRMin(), rmin), 0., 0., std::min(getRMax(), rmax), 0., 0.);
  short lmin, lmax;
  if (!getLayersRange(ray, lmin, lmax)) {
    LOGP(warn, "No layers found for {} < r < {}", rmin, rmax);
    return nullptr;
  }
  LOGP(info, "Will extract layers {}:{} (out of {} layers) for {} < r < {}", lmin, lmax, getNLayers(), rmin, rmax);
  MatLayerCylSet* copy = new MatLayerCylSet();
  int lrCount = 0;
  for (int il = lmin; il <= lmax; il++) {
    const auto& lr = getLayer(il);
    float drphi = lr.getDPhi() * (lr.getRMin() + lr.getRMax()) / 2. * 0.999;
    copy->addLayer(lr.getRMin(), lr.getRMax(), lr.getZMax(), lr.getDZ(), drphi);
    auto& lrNew = copy->getLayer(lrCount);
    for (int iz = 0; iz < lrNew.getNZBins(); iz++) {
      for (int ip = 0; ip < lrNew.getNPhiBins(); ip++) {
        lrNew.getCellPhiBin(ip, iz).set(lr.getCellPhiBin(ip, iz));
      }
    }
    lrCount++;
  }

  copy->finalizeStructures();
  copy->optimizePhiSlices(tolerance);
  copy->flatten();
  return copy;
}

// builds the voxel structure; dimX; dimY are the dimensions of the VoxelStructure
vecgeom::Flat2DVoxelHashMap<short, int, false> determineVoxelSet(MatLayerCyl const& layer, int dimX, int dimY, std::vector<std::pair<float, float>> &points, std::vector<short> &flatlookup ) {
  // determines the set of all voxel cells that have overlap with Cylinder
  auto Rmax = layer.getRMax();
  auto Rmin = layer.getRMin();
  std::cout << " This layer has radius " << Rmin << " to " << Rmax << "\n";
  using HM = typename vecgeom::Flat2DVoxelHashMap<short, int, false>;
  HM voxelmap(vecgeom::Vector2D<float>(-Rmax,-Rmax), vecgeom::Vector2D<float>(Rmax, Rmax), dimX, dimY);

  flatlookup.resize(dimX * dimY);

  auto PropList = [&voxelmap](auto key) {
    int s = 0;
    const auto ptr = voxelmap.getPropertiesGivenKey(key, s);
    return gsl::span<const short>(ptr, s);
  };

  auto addProp = [&flatlookup, &voxelmap,&layer,PropList](auto x, auto y) {
    float phi = o2::gpu::CAMath::ATan2(y, x);
    o2::math_utils::bringTo02Pi(phi);
    const auto voxelkey = voxelmap.getKey(x, y);
    const auto phislice = layer.getPhiSliceID(phi);
    bool exists = false;
    flatlookup[voxelkey] = phislice;
    for (const auto& existing : PropList(voxelkey)) {
      if (existing == phislice) {
        exists = true;
          break;
        }
    }
    if (!exists) {
      voxelmap.addPropertyForKey(voxelkey, phislice);
    }
  };

  // we go through row by row y = x
  // the starting point is determined through ray - circle crossings
  // also; here we will use a smaller grid with much higher granularity
  // say 10000 x 10000
  const int NVoxels1D_ = 5000; // number of voxels in one direction
  const float xmin = -Rmax - 0.2; // 
  const float Delta = 2 * std::abs(xmin) / NVoxels1D_; // 1000 = 2*XMIN
  const float invDelta = 1./ Delta;

  for (float y = xmin; y <= -xmin; y += Delta) {
    // determine the x - range to loop over ... with ray - circular layer intersection
    // std::cout << "Probing y " << y << "\n";
    Ray ray(xmin, y, 0, -xmin, y, 0);
    float crsmax1, crsmax2, crsmin1, crsmin2;
    bool outercross = ray.crossCircleR(layer.getRMax2(), crsmax1, crsmax2);
    bool innercross = ray.crossCircleR(layer.getRMin2(), crsmin1, crsmin2);
    
    if (outercross && !innercross) {
      auto startX = ray.getPos(std::min(crsmax1, crsmax2), 0);
      auto endX = ray.getPos(std::max(crsmax1, crsmax2), 0);
      // std::cout << y << " from " << startX << " to " << endX << "\n";
      for (float x = startX - Delta; x <= endX + Delta; x += Delta) {   
        points.push_back(std::pair<float, float>(x, y));
        addProp(x, y);
      }
    }
    if (outercross && innercross) {
      auto startX1 = ray.getPos(std::min(crsmax1, crsmax2), 0);
      auto endX2 = ray.getPos(std::max(crsmax1, crsmax2), 0);
      
      auto startX2 = ray.getPos(std::max(crsmin1, crsmin2), 0);
      auto endX1 = ray.getPos(std::min(crsmin1, crsmin2), 0);
    
      // std::cout << y << " from " << startX1 << " to " << endX1 << " and " << startX2 << " to " << endX2 << "\n";
      for (float x = startX1 - Delta; x <= endX1 + Delta; x += Delta) {   
        points.push_back(std::pair<float, float>(x, y));
        addProp(x, y);
      }
      for (float x = startX2 - Delta; x <= endX2 + Delta; x += Delta) {   
        points.push_back(std::pair<float, float>(x, y));
        addProp(x, y);
      }
    }
  }
  auto S1 = voxelmap.size();
  auto S2 = voxelmap.getPropSize();
  std::cout << "Filled voxel " << voxelmap.size() << " Prop Size " << voxelmap.getPropSize() << " " << (1.*S2)/(1.*S1) << "\n";
  return voxelmap;
}


// let's assume points belong to layer l for which we have a voxmap
// Let's benchmark the old vs the lookup solution
void benchmark(MatLayerCyl const& layer, std::vector<std::pair<float, float>> const& points, vecgeom::Flat2DVoxelHashMap<short, int, false> const& voxmap, std::vector<short> const& flatlookup) {
  {
   long counter = 0;
   TStopwatch timer;
   timer.Start();
   for (int i = 0; i <= points.size()*100; ++i) {
    auto& p = points[i % points.size()];
    auto phi = std::atan2(p.second, p.first);
    o2::math_utils::bringTo02Pi(phi);
    auto phiID = layer.getPhiSliceID(phi);
    counter += phiID;
   }
   timer.Stop();
   std::cout << " Conventional time " << timer.RealTime() << " Checksum " << counter << "\n";
  }
  {
   long counter = 0;
   TStopwatch timer;
   timer.Start();
   for (int i = 0; i <= points.size()*100; ++i) {
    auto& p = points[i % points.size()];
    int s;
    auto props = voxmap.getPropertiesGivenKey(voxmap.getKey(p.first, p.second), s);
    if (s == 1) {
      counter += *props;
    }
   }
   timer.Stop();
   std::cout << " Best new time " << timer.RealTime() << " Checksum " << counter << "\n";
  }
  {
   long counter = 0;
   TStopwatch timer;
   timer.Start();
   for (int i = 0; i <= points.size()*100; ++i) {
    auto& p = points[i % points.size()];
    auto phi = flatlookup[voxmap.getKey(p.first, p.second)];
    counter += phi;
   }
   timer.Stop();
   std::cout << " Best flat lookup new time " << timer.RealTime() << " Checksum " << counter << "\n";
  }
}

void MatLayerCylSet::analyseLayers() {
  
  auto analyseOneLayer = [this](int i) {
    std::cout << " ------ " << i << " ----- \n";
    auto& layer = getLayer(i);
    const auto Rlower = layer.getRMin();
    const auto Rupper = layer.getRMax();
    auto NSlices = layer.getNPhiSlices();
    std::vector<std::pair<float, float>> points;
    if (NSlices > 1) {
      float deltaXMin = 1000000.;
      float deltaYMin = 1000000.;
      int smalledAreaIndex = 0;
      for (int slice = 0; slice < NSlices; ++slice) {
        int binlower, binupper;
        layer.getNPhiBinsInSlice(slice, binlower, binupper);
        auto phiStart = layer.getPhiBinMin(binlower);
        auto phiEnd = layer.getPhiBinMax(binupper);
      
        auto y1 = Rlower * std::sin(phiStart);
        auto x1 = Rlower * std::cos(phiStart);
        auto y2 = Rupper * std::sin(phiStart);
        auto x2 = Rupper * std::cos(phiStart);
      
        auto y3 = Rlower * std::sin(phiEnd);
        auto x3 = Rlower * std::cos(phiEnd);
        auto y4 = Rupper * std::sin(phiEnd);
        auto x4 = Rupper * std::cos(phiEnd);

        deltaYMin = std::min(deltaYMin, std::abs(y4 - y1));
        deltaXMin = std::min(deltaXMin, std::abs(x2 - x3));
      }
      auto dimX = int(Rupper / deltaXMin);
      auto dimY = int(Rupper / deltaYMin);
      std::cout << " DeltaXMin " << deltaXMin << " " << " DeltaYMin " << deltaYMin << "\n";
      std::cout << dimX << " x " << dimY << "\n";
      std::vector<short> flatlookup;
      auto voxmap = determineVoxelSet(layer, std::min(1000,dimX), std::min(dimY,1000), points, flatlookup);
      benchmark(layer, points, voxmap, flatlookup);
    }
    else {
      std::cout << " Just 1 slice. Nothing to be done\n";
    }
  };
  for (int l = 40; l < 41; ++l) {
    analyseOneLayer(l);
  }
}

#endif
