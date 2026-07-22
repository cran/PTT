#include "helpers.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

double unifRand() {
  // Use R's RNG so set.seed() controls posterior sampling. Avoid exact 0 and 1
  // because downstream code takes logarithms of the draw.
  double draw = R::runif(0.0, 1.0);
  return std::min(
    std::max(draw, std::numeric_limits<double>::min()),
    std::nextafter(1.0, 0.0)
  );
}

double log_exp_x_plus_exp_y(double x, double y) {
  if (std::isinf(x)) return x > 0 ? x : y;
  if (std::isinf(y)) return y > 0 ? y : x;
  const double maximum = std::max(x, y);
  const double minimum = std::min(x, y);
  return maximum + std::log1p(std::exp(minimum - maximum));
}

unsigned long long pow2(int k) {

  return (unsigned long long) 1 << k;
}

void print_index(INDEX_TYPE & I, int level) {
  ushort x_curr = 0;
  ushort index_prev_var = 0;

  for (int i = 0; i < level; i++) {
    x_curr += I.var[i] - index_prev_var - 1;
    Rcpp::Rcout << "X" << x_curr << "=" << ((I.var[MAXVAR] >> i) & (ushort) 1) << ",";
    index_prev_var = I.var[i];
  }
  Rcpp::Rcout << I.var[MAXVAR];
}

void print_index_2(INDEX_TYPE & I, int level, int k) {
  ushort x_curr = 0;
  ushort index_prev_var = 0;
  ushort lower = 0;
  int x_curr_count = -1;

  for (int i = 0; i < level; i++) {
    if ( I.var[i] - index_prev_var - 1 > 0 ) { // next variable
      if (x_curr > 0) Rcpp::Rcout << ";";
      Rcpp::Rcout << "X" << x_curr << "_l=" << lower << ", X" << x_curr << "_u=" << lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;

      lower = 0;
      x_curr_count = 0;
    }
    else {
      x_curr_count++;
    }

    x_curr += I.var[i] - index_prev_var - 1;
    lower |= (((I.var[MAXVAR] >> i) & (ushort) 1)) << (k-x_curr_count-1);
    index_prev_var = I.var[i];
  }

  if (level > 0) {
    if (x_curr >0) Rcpp::Rcout << ";";
    Rcpp::Rcout << "X" << x_curr << "_l=" << lower << ", X" << x_curr << "_u=" << lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;
  }

  if (level == 0) {

  }
}

unsigned int Choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);

  unsigned long long result = 1;
  for (int i = 1; i <= k; i++) {
    const unsigned long long factor = n - k + i;
    if (result > ULLONG_MAX / factor) {
      Rcpp::stop("Requested tree is too large to index safely.");
    }
    result = result * factor / i;
    if (result > UINT_MAX) {
      Rcpp::stop("Requested tree is too large to index safely.");
    }
  }
  return static_cast<unsigned int>(result);
}



INDEX_TYPE init_index(int n,int level) {
  INDEX_TYPE init;

  for (int i=0; i < level; i++) {init.var[i] = i+1;}
  for (int i=level; i <= MAXVAR; i++) {init.var[i] = 0;}

  return init;
}

INDEX_TYPE make_child_index(INDEX_TYPE& I, unsigned short part_dim, int level, ushort which) {
  INDEX_TYPE child_index = I;
  unsigned short data = part_dim+1;
  int i;
  int j;
  int x_curr;
  int child_index_var_prev;

  if (level == 0) {
    x_curr=1;
    i=0;
    child_index_var_prev = 0;
  }

  else {
    x_curr = child_index.var[0]; // current dimension
    child_index_var_prev = child_index.var[0];
    i = 1;
  }

  while (i<MAXVAR) {
    while (child_index.var[i] >0 && data >= x_curr ) {
      x_curr += child_index.var[i] - child_index_var_prev - 1;
      child_index_var_prev = child_index.var[i];
      i++;

    }

    if (child_index.var[i] == 0 && data >= x_curr) {
      child_index.var[i] = data - x_curr + 1 + child_index_var_prev;
      j=i;
      i=MAXVAR;
    } else { // this corresponds to the first i such that data < x_curr
      for (int h = level; h >= i; h--) {
        child_index.var[h] = child_index.var[h-1]+1;
      }

      child_index.var[i-1] = child_index.var[i] - (x_curr - data + 1);

      j=i-1;
      i=MAXVAR;
    }
  }
  // Update the child bits using unsigned masks. The original expression left
  // shifted a negative promoted integer, which is undefined behavior in C++.
  const unsigned int all_ones = std::numeric_limits<ushort>::max();
  const unsigned int bits = I.var[MAXVAR];
  const unsigned int high_mask = (all_ones << (j + 1)) & all_ones;
  const unsigned int low_mask = ~(all_ones << j) & all_ones;
  child_index.var[MAXVAR] = static_cast<ushort>(
    ((bits << 1) & high_mask) |
      (static_cast<unsigned int>(which) << j) |
      (bits & low_mask)
  );

  return child_index;
}


INDEX_TYPE get_next_node(INDEX_TYPE& I, int p, int level) {

  INDEX_TYPE node = I;
  int i = level-1; int j = p + level -2;
  while (i>=0 && node.var[i] == j+1) {i--;j--;}
  if (i < 0) { //reach the end of nodes
    for (int h = 0; h <= MAXVAR; h++) {
      node.var[h]=0; //invalid node
    }
  } else {
    node.var[i] += 1;
    for (j=i+1;j<level;j++) {
      node.var[j] = node.var[i]+ j-i;
    }
  }
  node.var[MAXVAR] = 0;

  return node;
}

uint convert_to_inverse_base_2(double x, int k) {
  const double clipped = std::min(
    std::max(x, 0.0),
    std::nextafter(1.0, 0.0)
  );
  uint x_base_2 = (uint) floor(clipped * pow2(k));
  uint x_base_inverse_2 = 0;

  for (int i = 0; i < k; i++) {

    x_base_inverse_2 |= ((x_base_2 >> (k - 1 -i)) & 1) << i;
  }

  return x_base_inverse_2;
}

uint convert_to_base_2(double x, int k) {
  const double clipped = std::min(
    std::max(x, 0.0),
    std::nextafter(1.0, 0.0)
  );
  uint x_base_2 = (uint) floor(clipped * pow2(k));

  return x_base_2;
}
