FieldPerp extrap_sheath(const Field3D &var);
Field2D RMS_Z(const Field3D &var);
void spread(Field3D f);
Field3D integrate_dy(Field3D i);
void check_nonneg(Field3D &f);
bool limit_at_least(Field3D &f, BoutReal l, bool dothrow = false);
bool limit_at_most(Field3D &f, BoutReal l);

void limit_at_least_smooth(Field3D &f, BoutReal l);
void myNorm(Field3D f, BoutReal *data, int max, int counter);
bool read_equilibrium_file(BoutReal *data, const char *fname);
void set_equilibrium_value(Field3D &f, const char *fname, bool set_default,
                           BoutReal def_value);
void set_equilibrium_value(Field3D &f, const char *fname, bool set_default = false);
void set_nans(Field3D &f);
