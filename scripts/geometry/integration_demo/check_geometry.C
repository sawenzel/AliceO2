// Stage 3 of the integration demo: where did the CAD modules land, and do they overlap?
//
//   root -l -b -q 'check_geometry.C("o2sim_geometry.root", 1)'
//
// Second argument: also run TGeoManager::CheckOverlaps (slow). Its verdict on
// O2BVHSurfaceSolid is documented-unreliable (Handoff_IntegrationTest.md §8) -- a hint,
// never a verdict.
//
// The two CAD modules appear under `barrel` as the assemblies the converter emits:
// IRIS's top assembly is called "Unnamed" and Bagger's "Assembly" (both come from the CAD
// root label, which O2_CADtoTGeo.py does not rename).

static void subtreeBox(TGeoVolume* vol, const TGeoHMatrix& base, const char* label)
{
  Double_t lo[3] = {1e30, 1e30, 1e30}, hi[3] = {-1e30, -1e30, -1e30};
  long nleaf = 0;
  TGeoIterator it(vol);
  TGeoNode* nd;
  while ((nd = it.Next())) {
    if (nd->GetVolume()->IsAssembly() || !nd->GetVolume()->GetShape()) continue;
    TGeoHMatrix m = base * (*it.GetCurrentMatrix());
    auto* bb = (TGeoBBox*)nd->GetVolume()->GetShape();
    const Double_t* o = bb->GetOrigin();
    for (int i = 0; i < 8; i++) {
      Double_t l[3] = {o[0] + ((i & 1) ? 1 : -1) * bb->GetDX(),
                       o[1] + ((i & 2) ? 1 : -1) * bb->GetDY(),
                       o[2] + ((i & 4) ? 1 : -1) * bb->GetDZ()}, g[3];
      m.LocalToMaster(l, g);
      for (int k = 0; k < 3; k++) { if (g[k] < lo[k]) lo[k] = g[k]; if (g[k] > hi[k]) hi[k] = g[k]; }
    }
    nleaf++;
  }
  printf("WORLDBOX %-10s leaves=%4ld  x[%9.3f,%9.3f] y[%9.3f,%9.3f] z[%9.3f,%9.3f]\n",
         label, nleaf, lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);
}

void check_geometry(const char* geofile, int overlaps = 0, double ovlp_prec = 0.1)
{
  TGeoManager::Import(geofile);
  auto* mgr = gGeoManager;
  printf("GEOM %s\n", geofile);
  printf("COUNTS volumes=%d nodes=%d\n", mgr->GetListOfVolumes()->GetEntries(), mgr->GetNNodes());

  std::map<std::string, std::string> want = {{"Unnamed", "IRIS"}, {"Assembly", "BAGR"}};
  TGeoIterator it(mgr->GetTopVolume());
  TGeoNode* nd;
  while ((nd = it.Next())) {
    auto f = want.find(nd->GetVolume()->GetName());
    if (f == want.end()) continue;
    subtreeBox(nd->GetVolume(), *it.GetCurrentMatrix(), f->second.c_str());
    want.erase(f);
  }
  for (auto& kv : want) printf("WORLDBOX %-10s NOT FOUND\n", kv.second.c_str());

  std::map<std::string, int> cls;
  TIter nx(mgr->GetListOfVolumes()); TGeoVolume* v;
  while ((v = (TGeoVolume*)nx())) { if (v->GetShape()) cls[v->GetShape()->ClassName()]++; }
  for (auto& kv : cls) printf("SHAPE %-32s %d\n", kv.first.c_str(), kv.second);

  if (overlaps) {
    mgr->CheckOverlaps(ovlp_prec);
    auto* l = mgr->GetListOfOverlaps();
    printf("OVERLAPS prec=%g count=%d\n", ovlp_prec, l ? l->GetEntries() : 0);
    if (l) {
      for (int i = 0; i < l->GetEntries(); i++) {
        auto* ov = (TGeoOverlap*)l->At(i);
        printf("OVERLAP %.6f cm  %s\n", ov->GetOverlap(), ov->GetTitle());
      }
    }
  }
}
