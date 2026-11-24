"""Select the branch-local electrostatic wall operator in a fresh fine-grid build."""
from pathlib import Path


def select(source, target):
    text = Path(source).read_text()
    marker = '  if (unresolved) return -1;'
    if text.count(marker) != 1:
        raise ValueError('Unrecognized volume-cache refresh implementation')
    text = text.replace(marker, marker + '''
  int failed = 0;
  if (rank == 0) {
    extern int cube_project_branch_voltage(double *, const signed char *, double, double);
    extern int cube_yy_igrf_field(const double *, double, double *, double *);
    signed char *branch = malloc((size_t)count);
    if (!branch) failed=1;
    for (int n=0; !failed && n<count; n++) {
      int i=n/(CACHE_NTH*CACHE_NPH), j=(n/CACHE_NPH)%CACHE_NTH, k=n%CACHE_NPH;
      double radius=mi_state.pot_cache_r_min+i*dr, th=j*dth, ph=k*dph;
      double u[3]={sin(th)*cos(ph),sin(th)*sin(ph),cos(th)};
      double x[3]={radius*u[0],radius*u[1],radius*u[2]}, b[3], mag;
      if (cube_yy_igrf_field(x,time_s,b,&mag)) { failed=1; break; }
      branch[n]=(b[0]*u[0]+b[1]*u[1]+b[2]*u[2]<=0.) ? 1 : -1;
    }
    if (!failed) failed=cube_project_branch_voltage(mi_state.pot_cache_value,branch,
                         mi_state.pot_cache_r_min,mi_state.pot_cache_r_max);
    free(branch);
  }
  if (MPI_Bcast(&failed,1,MPI_INT,0,cache_comm)!=MPI_SUCCESS || failed) return -1;
  if (MPI_Bcast(mi_state.pot_cache_value,count,MPI_DOUBLE,0,cache_comm)!=MPI_SUCCESS) return -1;
''')
    target.write_text(text)
    return target
