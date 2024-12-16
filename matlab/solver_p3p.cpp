#include "/usr/include/eigen3/Eigen/Dense"
#include "mex.h"
#include <iostream>

bool solve_cubic_single_real(double c2, double c1, double c0, double &root) {
  double a = c1 - c2 * c2 / 3.0;
  double b = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1) / 27.0 + c0;
  double c = b * b / 4.0 + a * a * a / 27.0;
  if (c != 0) {
    if (c > 0) {
      c = std::sqrt(c);
      b *= -0.5;
      root = std::cbrt(b + c) + std::cbrt(b - c) - c2 / 3.0;
      return true;
    } else {
      c = 3.0 * b / (2.0 * a) * std::sqrt(-3.0 / a);
      root =
          2.0 * std::sqrt(-a / 3.0) * std::cos(std::acos(c) / 3.0) - c2 / 3.0;
    }
  } else {
    root = -c2 / 3.0 + (a != 0 ? (3.0 * b / a) : 0);
  }
  return false;
}

bool inline root2real(double b, double c, double &r1, double &r2) {
  double THRESHOLD = -1.0e-12;
  double v = b * b - 4.0 * c;
  if (v < THRESHOLD) {
    r1 = r2 = -0.5 * b;
    return v >= 0;
  }
  if (v > THRESHOLD && v < 0.0) {
    r1 = -0.5 * b;
    r2 = -2;
    return true;
  }

  double y = std::sqrt(v);
  if (b < 0) {
    r1 = 0.5 * (-b + y);
    r2 = 0.5 * (-b - y);
  } else {
    r1 = 2.0 * c / (-b + y);
    r2 = 2.0 * c / (-b - y);
  }
  return true;
}

inline std::array<Eigen::Vector3d, 2> compute_pq(Eigen::Matrix3d C) {
  std::array<Eigen::Vector3d, 2> pq;
  Eigen::Matrix3d C_adj;

  C_adj(0, 0) = C(1, 2) * C(2, 1) - C(1, 1) * C(2, 2);
  C_adj(1, 1) = C(0, 2) * C(2, 0) - C(0, 0) * C(2, 2);
  C_adj(2, 2) = C(0, 1) * C(1, 0) - C(0, 0) * C(1, 1);
  C_adj(0, 1) = C(0, 1) * C(2, 2) - C(0, 2) * C(2, 1);
  C_adj(0, 2) = C(0, 2) * C(1, 1) - C(0, 1) * C(1, 2);
  C_adj(1, 0) = C_adj(0, 1);
  C_adj(1, 2) = C(0, 0) * C(1, 2) - C(0, 2) * C(1, 0);
  C_adj(2, 0) = C_adj(0, 2);
  C_adj(2, 1) = C_adj(1, 2);

  Eigen::Vector3d v;
  if (C_adj(0, 0) > C_adj(1, 1)) {
    if (C_adj(0, 0) > C_adj(2, 2)) {
      v = C_adj.col(0) / std::sqrt(C_adj(0, 0));
    } else {
      v = C_adj.col(2) / std::sqrt(C_adj(2, 2));
    }
  } else if (C_adj(1, 1) > C_adj(2, 2)) {
    v = C_adj.col(1) / std::sqrt(C_adj(1, 1));
  } else {
    v = C_adj.col(2) / std::sqrt(C_adj(2, 2));
  }

  C(0, 1) -= v(2);
  C(0, 2) += v(1);
  C(1, 2) -= v(0);
  C(1, 0) += v(2);
  C(2, 0) -= v(1);
  C(2, 1) += v(0);

  pq[0] = C.col(0);
  pq[1] = C.row(0);

  return pq;
}

// Performs a few newton steps on the equations
inline void refine_lambda(double &lambda1, double &lambda2, double &lambda3,
                          const double a12, const double a13, const double a23,
                          const double b12, const double b13,
                          const double b23) {

  for (int iter = 0; iter < 5; ++iter) {
    double r1 = (lambda1 * lambda1 - 2.0 * lambda1 * lambda2 * b12 +
                 lambda2 * lambda2 - a12);
    double r2 = (lambda1 * lambda1 - 2.0 * lambda1 * lambda3 * b13 +
                 lambda3 * lambda3 - a13);
    double r3 = (lambda2 * lambda2 - 2.0 * lambda2 * lambda3 * b23 +
                 lambda3 * lambda3 - a23);
    if (std::abs(r1) + std::abs(r2) + std::abs(r3) < 1e-10)
      return;
    double x11 = lambda1 - lambda2 * b12;
    double x12 = lambda2 - lambda1 * b12;
    double x21 = lambda1 - lambda3 * b13;
    double x23 = lambda3 - lambda1 * b13;
    double x32 = lambda2 - lambda3 * b23;
    double x33 = lambda3 - lambda2 * b23;
    double detJ = 0.5 / (x11 * x23 * x32 +
                         x12 * x21 * x33); // half minus inverse determinant
    // This uses the closed form of the inverse for the jacobean.
    // Due to the zero elements this actually becomes quite nice.
    lambda1 += (-x23 * x32 * r1 - x12 * x33 * r2 + x12 * x23 * r3) * detJ;
    lambda2 += (-x21 * x33 * r1 + x11 * x33 * r2 - x11 * x23 * r3) * detJ;
    lambda3 += (x21 * x32 * r1 - x11 * x32 * r2 - x12 * x21 * r3) * detJ;
  }
}

Eigen::MatrixXd solver_p3p_real(const Eigen::VectorXd &data) {

  const double *d = data.data();
  Eigen::MatrixXd x1(3, 3);
  Eigen::MatrixXd x2(3, 3);
  for (int i = 0; i < 3; i++) {
    x1.col(i) << d[i * 6], d[i * 6 + 1], d[i * 6 + 2];
    x2.col(i) << d[i * 6 + 3], d[i * 6 + 4], d[i * 6 + 5];
  }

  Eigen::Vector3d X01 = x1.col(0) - x1.col(1);
  Eigen::Vector3d X02 = x1.col(0) - x1.col(2);
  Eigen::Vector3d X12 = x1.col(1) - x1.col(2);

  double a01 = X01.squaredNorm();
  double a02 = X02.squaredNorm();
  double a12 = X12.squaredNorm();

  std::array<Eigen::Vector3d, 3> X = {x1.col(0), x1.col(1), x1.col(2)};
  std::array<Eigen::Vector3d, 3> x = {x2.col(0), x2.col(1), x2.col(2)};

  if (a01 > a02) {
    if (a01 > a12) {
      std::swap(x[0], x[2]);
      std::swap(X[0], X[2]);
      std::swap(a01, a12);
      X01 = -X12;
      X02 = -X02;
    }
  } else if (a02 > a12) {
    std::swap(x[0], x[1]);
    std::swap(X[0], X[1]);
    std::swap(a02, a12);
    X01 = -X01;
    X02 = X12;
  }

  const double a12d = 1.0 / a12;
  const double a = a01 * a12d;
  const double b = a02 * a12d;

  const double m01 = x[0].dot(x[1]);
  const double m02 = x[0].dot(x[2]);
  const double m12 = x[1].dot(x[2]);

  // Ugly parameters to simplify the calculation
  const double m12sq = -m12 * m12 + 1.0;
  const double m02sq = -1.0 + m02 * m02;
  const double m01sq = -1.0 + m01 * m01;
  const double ab = a * b;
  const double bsq = b * b;
  const double asq = a * a;
  const double m013 = -2.0 + 2.0 * m01 * m02 * m12;
  const double bsqm12sq = bsq * m12sq;
  const double asqm12sq = asq * m12sq;
  const double abm12sq = 2.0 * ab * m12sq;

  const double k3_inv = 1.0 / (bsqm12sq + b * m02sq);
  const double k2 =
      k3_inv * ((-1.0 + a) * m02sq + abm12sq + bsqm12sq + b * m013);
  const double k1 =
      k3_inv * (asqm12sq + abm12sq + a * m013 + (-1.0 + b) * m01sq);
  const double k0 = k3_inv * (asqm12sq + a * m01sq);
  double s;
  bool G = solve_cubic_single_real(k2, k1, k0, s);

  Eigen::Matrix3d C;
  C(0, 0) = -a + s * (1 - b);
  C(0, 1) = -m02 * s;
  C(0, 2) = a * m12 + b * m12 * s;
  C(1, 0) = C(0, 1);
  C(1, 1) = s + 1;
  C(1, 2) = -m01;
  C(2, 0) = C(0, 2);
  C(2, 1) = C(1, 2);
  C(2, 2) = -a - b * s + 1;

  std::array<Eigen::Vector3d, 2> pq = compute_pq(C);

  double d0, d1, d2;

  Eigen::Matrix3d XX;

  XX << X01, X02, X01.cross(X02);
  XX = XX.inverse().eval();

  Eigen::Vector3d v1, v2;
  Eigen::Matrix3d YY;

  int step = 0;
  Eigen::MatrixXd sols(3, 16);
  for (int i = 0; i < 2; ++i) {
    // [p0 p1 p2] * [1; x; y] = 0, or [p0 p1 p2] * [d2; d0; d1] = 0
    double p0 = pq[i](0);
    double p1 = pq[i](1);
    double p2 = pq[i](2);
    // here we run into trouble if p0 is zero,
    // so depending on which is larger, we solve for either d0 or d1
    // The case p0 = p1 = 0 is degenerate and can be ignored
    bool switch_12 = std::abs(p0) <= std::abs(p1);

    if (switch_12) {
      // eliminate d0
      double w0 = -p0 / p1;
      double w1 = -p2 / p1;
      double ca = 1.0 / (w1 * w1 - b);
      double cb = 2.0 * (b * m12 - m02 * w1 + w0 * w1) * ca;
      double cc = (w0 * w0 - 2 * m02 * w0 - b + 1.0) * ca;
      double taus[2];
      if (!root2real(cb, cc, taus[0], taus[1]))
        continue;
      for (double tau : taus) {
        if (tau <= 0)
          continue;
        // positive only
        d2 = std::sqrt(a12 / (tau * (tau - 2.0 * m12) + 1.0));
        d1 = tau * d2;
        d0 = (w0 * d2 + w1 * d1);
        if (d0 < 0)
          continue;

        refine_lambda(d0, d1, d2, a01, a02, a12, m01, m02, m12);
        v1 = d0 * x[0] - d1 * x[1];
        v2 = d0 * x[0] - d2 * x[2];
        YY << v1, v2, v1.cross(v2);
        Eigen::Matrix3d R = YY * XX;
        sols.block<3, 4>(0, step) << R, d0 * x[0] - R * X[0];
        step = step + 4;
      }
    } else {
      double w0 = -p1 / p0;
      double w1 = -p2 / p0;
      double ca = 1.0 / (-a * w1 * w1 + 2 * a * m12 * w1 - a + 1);
      double cb = 2 * (a * m12 * w0 - m01 - a * w0 * w1) * ca;
      double cc = (1 - a * w0 * w0) * ca;

      double taus[2];
      if (!root2real(cb, cc, taus[0], taus[1]))
        continue;
      for (double tau : taus) {
        if (tau <= 0)
          continue;
        d0 = std::sqrt(a01 / (tau * (tau - 2.0 * m01) + 1.0));
        d1 = tau * d0;
        d2 = w0 * d0 + w1 * d1;

        if (d2 < 0)
          continue;

        refine_lambda(d0, d1, d2, a01, a02, a12, m01, m02, m12);
        v1 = d0 * x[0] - d1 * x[1];
        v2 = d0 * x[0] - d2 * x[2];
        YY << v1, v2, v1.cross(v2);
        Eigen::Matrix3d R = YY * XX;
        sols.block<3, 4>(0, step) << R, d0 * x[0] - R * X[0];
        step = step + 4;
      }
    }

    if (step > 0 && G)
      break;
  }

  sols.conservativeResize(3, step);

  return sols;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (mxGetNumberOfElements(prhs[0]) != 18) {
    mexErrMsgIdAndTxt("incorrectSize", "Input size must be double 18.");
  }
  const Eigen::VectorXd data =
      Eigen::Map<const Eigen::VectorXd>(mxGetPr(prhs[0]), 18);
  Eigen::MatrixXd sols = solver_p3p_real(data);
  plhs[0] = mxCreateDoubleMatrix(sols.rows(), sols.cols(), mxREAL);
  double *zr = mxGetPr(plhs[0]);
  for (Eigen::Index i = 0; i < sols.size(); i++) {
    zr[i] = sols(i);
  }
}