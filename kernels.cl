#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define NSPEEDS         9

typedef struct
{
  float speeds[NSPEEDS];
} t_speed;

kernel void accelerate_flow(global t_speed* cells,
                            global int* obstacles,
                            int nx, int ny,
                            float density, float accel)
{
  /* compute weighting factors */
  float w1 = density * accel / 9.0;
  float w2 = density * accel / 36.0;

  /* modify the 2nd row of the grid */
  int jj = ny - 2;

  /* get column index */
  int ii = get_global_id(0);

  /* if the cell is not occupied and
  ** we don't send a negative density */
  if (!obstacles[ii + jj* nx]
      && (cells[ii + jj* nx].speeds[3] - w1) > 0.f
      && (cells[ii + jj* nx].speeds[6] - w2) > 0.f
      && (cells[ii + jj* nx].speeds[7] - w2) > 0.f)
  {
    /* increase 'east-side' densities */
    cells[ii + jj* nx].speeds[1] += w1;
    cells[ii + jj* nx].speeds[5] += w2;
    cells[ii + jj* nx].speeds[8] += w2;
    /* decrease 'west-side' densities */
    cells[ii + jj* nx].speeds[3] -= w1;
    cells[ii + jj* nx].speeds[6] -= w2;
    cells[ii + jj* nx].speeds[7] -= w2;
  }
}

kernel void propagate(global t_speed* cells,
                      global t_speed* tmp_cells,
                      global int* obstacles,
                      int nx, int ny)
{
  /* get column and row indices */
  int ii = get_global_id(0);
  int jj = get_global_id(1);

  /* determine indices of axis-direction neighbours
  ** respecting periodic boundary conditions (wrap around) */
  int y_n = (jj + 1) % ny;
  int x_e = (ii + 1) % nx;
  int y_s = (jj == 0) ? (jj + ny - 1) : (jj - 1);
  int x_w = (ii == 0) ? (ii + nx - 1) : (ii - 1);

  /* propagate densities from neighbouring cells, following
  ** appropriate directions of travel and writing into
  ** scratch space grid */
  tmp_cells[ii + jj*nx].speeds[0] = cells[ii + jj*nx].speeds[0]; /* central cell, no movement */
  tmp_cells[ii + jj*nx].speeds[1] = cells[x_w + jj*nx].speeds[1]; /* east */
  tmp_cells[ii + jj*nx].speeds[2] = cells[ii + y_s*nx].speeds[2]; /* north */
  tmp_cells[ii + jj*nx].speeds[3] = cells[x_e + jj*nx].speeds[3]; /* west */
  tmp_cells[ii + jj*nx].speeds[4] = cells[ii + y_n*nx].speeds[4]; /* south */
  tmp_cells[ii + jj*nx].speeds[5] = cells[x_w + y_s*nx].speeds[5]; /* north-east */
  tmp_cells[ii + jj*nx].speeds[6] = cells[x_e + y_s*nx].speeds[6]; /* north-west */
  tmp_cells[ii + jj*nx].speeds[7] = cells[x_e + y_n*nx].speeds[7]; /* south-west */
  tmp_cells[ii + jj*nx].speeds[8] = cells[x_w + y_n*nx].speeds[8]; /* south-east */
}
