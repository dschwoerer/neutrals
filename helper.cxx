#include <bout.hxx>

FieldPerp extrap_sheath(const Field3D &var) {
  // 3rd order extrapolation of a field to find the sheath value
  if (var.getLocation() != CELL_CENTRE)
    throw BoutException("extrap_sheath needs CELL_CENTRE location");
  FieldPerp result;
  result = 0.375 * sliceXZ(var, mesh->yend + 1) + 0.75 * sliceXZ(var, mesh->yend) -
           0.125 * sliceXZ(var, mesh->yend - 1);
  return result;
}

Field2D RMS_Z(const Field3D &var) {

  Field2D result = 0.;

  for (int jx = 0; jx < mesh->LocalNx; jx++) {
    for (int jy = 0; jy < mesh->LocalNy; jy++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        result(jx, jy) += var(jx, jy, jz) * var(jx, jy, jz);
      }
      result(jx, jy) /= mesh->LocalNz;
    }
  }
  return sqrt(result);
}

void spread(Field3D f) {
  for (int y = 0; y < mesh->LocalNy; ++y) {
    BoutReal fc = f(mesh->xstart, y, 0);
    for (int x = 0; x < mesh->LocalNx; ++x) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        f(x, y, z) = fc;
      }
    }
  }
}

void check_nonneg(Field3D &f) {
  BoutReal k;
  for (auto i = f.begin(); !i.done(); ++i) {
    if (f[i] <= 0) {
      k = f[i];
      throw BoutException("found neg value (%g) at %d", k, i);
    }
  }
}

void limit_at_least_smooth(Field3D &f, BoutReal limit) {
  int max = mesh->LocalNx * mesh->LocalNy * mesh->LocalNz;
  BoutReal limit2 = limit * 2;
  BoutReal *data = f(0, 0);
  for (int i = 0; i < max; ++i) {
    if (data[i] < limit2) {
      data[i] = data[i] / 2 + limit;
    }
  }
}

bool limit_at_least(Field3D &f, BoutReal l, bool dothrow) {
  bool changed = false;
  for (auto i = f.begin(); !i.done(); ++i) {
    if (f[i] < l) {
      f[i] = l;
      changed = true;
    }
  }
  return changed;
}

bool limit_at_most(Field3D &f, BoutReal l) {
  bool changed = false;
  for (auto i = f.begin(); !i.done(); ++i) {
    if (f[i] > l) {
      f[i] = l;
      changed = true;
    }
  }
  return changed;
}

void myNorm(Field3D f, BoutReal *data, int max, int counter) {
  BoutReal norm2 = 0;
  BoutReal norminf = 0;
  for (int i = 2; i < mesh->LocalNx - 2; i++) {
    for (int j = 2; j < mesh->LocalNy - 2; ++j) {
      for (int k = 0; k < mesh->LocalNz; ++k) {
        BoutReal cur = abs(f(i, j, k));
        norm2 += cur * cur;
        norminf = norminf > cur ? norminf : cur;
      }
    }
  }
  data[counter] = norm2;
  data[counter + max] = norminf;
}

bool read_equilibrium_file(BoutReal *data, const char *fname) {
  std::string datadir;
  Options::getRoot()->get("datadir", datadir, "dafuq", false);
  int max = 512;
  char *fullname = new char[max];
  int size = snprintf(fullname, max, "%s/equilibrium/%s", datadir.c_str(), fname);
  if (size > max) {
    delete[] fullname;
    fullname = new char[size + 20];
    snprintf(fullname, size + 10, "%s/equilibrium/%s", datadir.c_str(), fname);
  }

  FILE *file = fopen(fullname, "rb");
  delete[] fullname;
  if (file == NULL) {
    output.write("\tWarning: Background file `%s` not found.\n", fname);
    return false;
  } else {
    fread(data, sizeof(BoutReal), mesh->GlobalNy, file);
    output.write("\tInfo: Background read from `%s`.\n", fname);
    return true;
  }
}

void set_equilibrium_value(Field3D &f, const char *fname, bool set_default,
                           BoutReal def_value) {
  f.allocate();
  BoutReal *data = NULL;
  data = new BoutReal[mesh->GlobalNy];
  if (read_equilibrium_file(data, fname)) {
    for (int i = 0; i < mesh->LocalNx; i++) {
      for (int j = 0; j < mesh->LocalNy; j++) {
        int jglobal = mesh->YGLOBAL(j) + 2;
        for (int k = 0; k < mesh->LocalNz; k++) {
          f(i, j, k) = data[jglobal];
        }
      }
    }
  } else if (set_default) {
#pragma omp parallel
    for (auto i = f.begin(); !i.done(); ++i) {
      f[i] = def_value;
    }
  }
  delete[] data;
}

void set_equilibrium_value(Field3D &f, const char *fname, bool set_default) {
  if (set_default) {
    throw BoutException("No default value given, but is tried to set!");
  } else {
    set_equilibrium_value(f, fname, set_default, 0);
  }
}

void set_nans(Field3D &f) {
  f.allocate();
  for (int x = 0; x < mesh->LocalNx; x++) {
    for (int y = 0; y < mesh->LocalNy; y++) {
      for (int z = 0; z < mesh->LocalNz; z++) {
        f(x, y, z) = nan("");
      }
    }
  }
}
