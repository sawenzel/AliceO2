/*
 * VoxelContainer.h
 *
 *  Created on: Feb 8, 2019
 *      Author: swenzel
 */

#ifndef COMMON_MATHUTILS_INCLUDE_MATHUTILS_VOXELCONTAINER_H_
#define COMMON_MATHUTILS_INCLUDE_MATHUTILS_VOXELCONTAINER_H_

#include <functional>

namespace o2 {

// a container keeping 3D voxels;
// essentially a regular grid over a finite dimensional space where
// each cell is equally sized and holds a value of type VoxelContent

template <typename VoxelContent>
class VoxelContainer {
public:
  VoxelContainer() = default;

  // returns the voxel content at that coordinate as reference
  VoxelContent& getVoxelContent(float x, float y, float z);
  // returns the voxel content at that bins (not checking boundaries) as reference
  VoxelContent& getVoxelContentByBin(int binx, int biny, int binz);

  // returns the voxel bin along some coordinate (might be negative or > dimension)
  int getBinX(float x) const { return (x - mOriginX) * mInvLx; }
  int getBinY(float y) const { return (y - mOriginY) * mInvLy; }
  int getBinZ(float z) const { return (z - mOriginZ) * mInvLz; }

  // inits the container into dimX... bins where the total space should span Lx and
  // start at some origin
  void init(int dimX, int dimY, int dimZ, float Lx, float Ly, float Lz,
            float originx, float originy, float originz);

  // inits the container into dimX... bins where the total space should span Lx and
  // assuming equal extent around the 0 origin
  void init(int dimX, int dimY, int dimZ, float Lx, float Ly, float Lz)
  {
    init(dimX, dimY, dimZ, Lx, Ly, Lz, -Lx / 2., -Ly / 2., -Lz / 2.);
  }

  // simple (linear) voxel visitor executing some hook (which is given the current
  // bin number)
  void visitAllVoxels(std::function<void(int, int, int)> hook) const;

 private:
  inline int toIndex(int binx, int biny, int binz)
  {
    return mDimX * mDimY * binz + mDimX * biny + binx;
  }

  int mDimX = 0; // number of voxels in x direction
  int mDimY = 0; // number of voxels in y direction
  int mDimZ = 0; // number of voxels in z direction
  float mLx = 0.f; // spatial size of voxel in x direction
  float mLy = 0.f; // spatial size of voxel in x direction
  float mLz = 0.f; // spatial size of voxel in x direction
  float mInvLx = 1.f; // inverse spatial size of voxel in x direction
  float mInvLy = 1.f; // inverse spatial size of voxel in x direction
  float mInvLz = 1.f; // inverse spatial size of voxel in x direction
  float mOriginX = 0.f; // location offset x
  float mOriginY = 0.f;
  float mOriginZ = 0.f;

  // the actual container in linear index space (to be seen how
  // we can make this configurable)
  std::vector<VoxelContent> mVoxels;
};

template <typename VoxelContent>
inline void VoxelContainer<VoxelContent>::init(int dimX, int dimY, int dimZ, float Lx, float Ly, float Lz,
                                               float originx, float originy, float originz)
{
  mOriginX = originx;
  mOriginY = originy;
  mOriginZ = originz;
  mDimX = dimX;
  mDimY = dimY;
  mDimZ = dimZ;

  // size per voxel
  mLx = Lx / mDimX;
  mLy = Ly / mDimY;
  mLz = Lz / mDimZ;

  // inverse size per voxel
  mInvLx = 1. / mLx;
  mInvLy = 1. / mLy;
  mInvLz = 1. / mLz;

  mVoxels.resize(mDimX*mDimY*mDimZ);
}

template <typename VoxelContent>
inline VoxelContent& VoxelContainer<VoxelContent>::getVoxelContentByBin(int binx, int biny, int binz)
{
  return mVoxels[toIndex(binx, biny, binz)];
}

template <typename VoxelContent>
inline VoxelContent& VoxelContainer<VoxelContent>::getVoxelContent(float x, float y, float z)
{
  int binx = getBinX(x);
  int biny = getBinY(y);
  int binz = getBinZ(y);
  return getVoxelContentByBin(binx, biny, binz);
}

template <typename VoxelContent>
void VoxelContainer<VoxelContent>::visitAllVoxels(std::function<void(int,int,int)> hook) const {
  for (int bz = 0; bz < mDimZ; ++bz) {
    for (int by = 0; by < mDimY; ++by) {
      for (int bx = 0; bx < mDimX; ++bx) {
        hook(bx, by, bz);
      }
    }
  }
}

} // namespace o2

#endif /* COMMON_MATHUTILS_INCLUDE_MATHUTILS_VOXELCONTAINER_H_ */
