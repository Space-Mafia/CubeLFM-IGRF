static double cubic_ramp_down(double coordinate, double start,
                              double length) {
  const double scaled = (coordinate - start) / length;
  if (scaled <= 0.0) {
    return 1.0;
  }
  if (scaled >= 1.0) {
    return 0.0;
  }
  return 1.0 - 3.0 * scaled * scaled + 2.0 * scaled * scaled * scaled;
}
