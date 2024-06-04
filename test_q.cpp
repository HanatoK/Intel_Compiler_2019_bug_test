#include <stdexcept>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cmath>

#if defined(COLVARS_LAMMPS)
#include "math_eigen_impl.h"
#elif defined(EIGEN3)
#include "eigen3/Eigen/Dense"
#else
#include "nr_jacobi.h"
#endif

void split_string(const std::string& data, const std::string& delim, std::vector<std::string>& dest) {
    size_t index = 0, new_index = 0;
    std::string tmpstr;
    while (index != data.length()) {
        new_index = data.find(delim, index);
        if (new_index != std::string::npos) tmpstr = data.substr(index, new_index - index);
        else tmpstr = data.substr(index, data.length());
        if (!tmpstr.empty()) {
            dest.push_back(tmpstr);
        }
        if (new_index == std::string::npos) break;
        index = new_index + 1;
    }
}

struct Coordinate {
  double x;
  double y;
  double z;
};

using AtomGroup = std::vector<Coordinate>;

AtomGroup load_xyz(const std::string& filename) {
  std::ifstream ifs(filename.c_str());
  std::string line;
  AtomGroup result;
  std::vector<std::string> fields;
  while (std::getline(ifs, line)) {
    fields.clear();
    split_string(line, " ", fields);
    if (fields.size() == 3) {
      Coordinate coor;
      coor.x = std::stod(fields[0]);
      coor.y = std::stod(fields[1]);
      coor.z = std::stod(fields[2]);
      result.push_back(coor);
    }
  }
  return result;
}

void bring_to_center(AtomGroup& ag) {
  Coordinate center{0, 0, 0};
  for (auto it = ag.begin(); it != ag.end(); ++it) {
    center.x += it->x;
    center.y += it->y;
    center.z += it->z;
  }
  center.x /= ag.size();
  center.y /= ag.size();
  center.z /= ag.size();
  for (auto it = ag.begin(); it != ag.end(); ++it) {
    it->x -= center.x;
    it->y -= center.y;
    it->z -= center.z;
  }
}

void print_matrix_transpose(const double mat[4][4]) {
  std::printf("%15.7e %15.7e %15.7e %15.7e\n", mat[0][0], mat[1][0], mat[2][0], mat[3][0]);
  std::printf("%15.7e %15.7e %15.7e %15.7e\n", mat[0][1], mat[1][1], mat[2][1], mat[3][1]);
  std::printf("%15.7e %15.7e %15.7e %15.7e\n", mat[0][2], mat[1][2], mat[2][2], mat[3][2]);
  std::printf("%15.7e %15.7e %15.7e %15.7e\n", mat[0][3], mat[1][3], mat[2][3], mat[3][3]);
}

struct eigen_result {
  double eigenvalues[4];
  double eigenvectors[4][4];
};

eigen_result optimal_rotation(
  const AtomGroup& ag, const AtomGroup& ag_ref) {
  double mat_R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (size_t i = 0; i < ag.size(); ++i) {
    mat_R[0][0] += ag[i].x * ag_ref[i].x;
    mat_R[0][1] += ag[i].x * ag_ref[i].y;
    mat_R[0][2] += ag[i].x * ag_ref[i].z;
    mat_R[1][0] += ag[i].y * ag_ref[i].x;
    mat_R[1][1] += ag[i].y * ag_ref[i].y;
    mat_R[1][2] += ag[i].y * ag_ref[i].z;
    mat_R[2][0] += ag[i].z * ag_ref[i].x;
    mat_R[2][1] += ag[i].z * ag_ref[i].y;
    mat_R[2][2] += ag[i].z * ag_ref[i].z;
  }
  double F[4][4];
  F[0][0] = mat_R[0][0] + mat_R[1][1] + mat_R[2][2]  ;
  F[0][1] = mat_R[1][2] - mat_R[2][1]                ;
  F[0][2] = mat_R[2][0] - mat_R[0][2]                ;
  F[0][3] = mat_R[0][1] - mat_R[1][0]                ;
  F[1][0] = F[0][1]                                  ;
  F[1][1] = mat_R[0][0] - mat_R[1][1] - mat_R[2][2]  ;
  F[1][2] = mat_R[0][1] + mat_R[1][0]                ;
  F[1][3] = mat_R[0][2] + mat_R[2][0]                ;
  F[2][0] = F[0][2]                                  ;
  F[2][1] = F[1][2]                                  ;
  F[2][2] = -mat_R[0][0] + mat_R[1][1] - mat_R[2][2] ;
  F[2][3] = mat_R[1][2] + mat_R[2][1]                ;
  F[3][0] = F[0][3]                                  ;
  F[3][1] = F[1][3]                                  ;
  F[3][2] = F[2][3]                                  ;
  F[3][3] = -mat_R[0][0] - mat_R[1][1] + mat_R[2][2] ;
  double w[4] = {0, 0, 0, 0};
  double v[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
#if defined(COLVARS_LAMMPS)
  MathEigen::Jacobi<double,
                    double[4],
                    double[4][4]> ecalc(4);
  ecalc.Diagonalize(F, w, v);
#elif defined (EIGEN3)
  Eigen::Matrix4d F_eigen;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      F_eigen(i, j) = F[i][j];
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigensolver;
  eigensolver.compute(F_eigen);
  for (size_t i = 0; i < 4; ++i) {
    w[3-i] = eigensolver.eigenvalues().col(0)[i];
    for (size_t j = 0; j < 4; ++j) {
      v[3-i][j] = eigensolver.eigenvectors().col(i)[j];
    }
  }
#else
  int nrot = 0;
  NR_Jacobi::jacobi(F, w, v, &nrot);
  NR_Jacobi::eigsrt(w, v);
  NR_Jacobi::transpose(v);
#endif
  eigen_result result;
  std::memcpy(&(result.eigenvalues), w, 4*sizeof(double));
  std::memcpy(&(result.eigenvectors), v, 4*4*sizeof(double));
  return result;
}

void test1() {
  AtomGroup protein_ref = load_xyz("refpos/protein.xyz");
  bring_to_center(protein_ref);
  AtomGroup protein_eq = load_xyz("traj/protein.xyz");
  bring_to_center(protein_eq);
  const eigen_result align = optimal_rotation(protein_eq, protein_ref);
  const auto max_eigenvector = align.eigenvectors[0];
  std::printf("Eigenvalues: %15.7e %15.7e %15.7e %15.7e\n",
              align.eigenvalues[0],
              align.eigenvalues[1],
              align.eigenvalues[2],
              align.eigenvalues[3]);
  std::printf("Eigenvectors:\n");
  print_matrix_transpose(align.eigenvectors);
  std::printf("Leading eigenvector: %15.7e %15.7e %15.7e %15.7e\n",
              max_eigenvector[0],
              max_eigenvector[1],
              max_eigenvector[2],
              max_eigenvector[3]);
}

void get_rot_matrix(const double q[4], double rotation_matrix[3][3]) {
  rotation_matrix[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  rotation_matrix[0][1] = 2.0 * (q[1] * q[2] - q[0] * q[3])                    ;
  rotation_matrix[0][2] = 2.0 * (q[1] * q[3] + q[0] * q[2])                    ;
  rotation_matrix[1][0] = 2.0 * (q[1] * q[2] + q[0] * q[3])                    ;
  rotation_matrix[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  rotation_matrix[1][2] = 2.0 * (q[2] * q[3] - q[0] * q[1])                    ;
  rotation_matrix[2][0] = 2.0 * (q[1] * q[3] - q[0] * q[2])                    ;
  rotation_matrix[2][1] = 2.0 * (q[2] * q[3] + q[0] * q[1])                    ;
  rotation_matrix[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}

void rotate(AtomGroup& ag, const double rotation_matrix[3][3]) {
  for (size_t i = 0; i < ag.size(); ++i) {
    double new_x = rotation_matrix[0][0] * ag[i].x + rotation_matrix[0][1] * ag[i].y + rotation_matrix[0][2] * ag[i].z;
    double new_y = rotation_matrix[1][0] * ag[i].x + rotation_matrix[1][1] * ag[i].y + rotation_matrix[1][2] * ag[i].z;
    double new_z = rotation_matrix[2][0] * ag[i].x + rotation_matrix[2][1] * ag[i].y + rotation_matrix[2][2] * ag[i].z;
    ag[i] = Coordinate{new_x, new_y, new_z};
  }
}

void test2() {
  AtomGroup protein_ref = load_xyz("refpos/protein.xyz");
  bring_to_center(protein_ref);
  AtomGroup protein_eq = load_xyz("traj/protein.xyz");
  bring_to_center(protein_eq);
  const eigen_result align = optimal_rotation(protein_eq, protein_ref);
  const auto max_eigenvector = align.eigenvectors[0];
  double rotation_matrix[3][3];
  get_rot_matrix(max_eigenvector, rotation_matrix);
  // Load the ligand
  AtomGroup ligand_ref = load_xyz("refpos/ligand.xyz");
  bring_to_center(ligand_ref);
  AtomGroup ligand_eq = load_xyz("traj/ligand.xyz");
  bring_to_center(ligand_eq);
  rotate(ligand_eq, rotation_matrix);
  const eigen_result rot_lig = optimal_rotation(ligand_eq, ligand_ref);
  const auto v = rot_lig.eigenvectors[0];
  std::printf("Eigenvalues: %15.7e %15.7e %15.7e %15.7e\n",
              rot_lig.eigenvalues[0],
              rot_lig.eigenvalues[1],
              rot_lig.eigenvalues[2],
              rot_lig.eigenvalues[3]);
  std::printf("Eigenvectors:\n");
  print_matrix_transpose(rot_lig.eigenvectors);
  std::printf("Leading eigenvector: %15.7e %15.7e %15.7e %15.7e\n",
              v[0], v[1], v[2], v[3]);
  // Calculate Euler angles
  const double theta = std::asin(2.0 * (v[0] * v[2] - v[3] * v[1])) * 180.0 / M_PI;
  const double phi = std::atan2(2.0 * (v[0] * v[1] + v[2] * v[3]), 1.0 - 2.0 * (v[1] * v[1] + v[2] * v[2])) * 180.0 / M_PI;
  const double psi = std::atan2(2.0 * (v[0] * v[3] + v[1] * v[2]), 1.0 - 2.0 * (v[2] * v[2] + v[3] * v[3])) * 180.0 / M_PI;
  std::printf("theta = %15.7e phi = %15.7e psi = %15.7e\n", theta, phi, psi);
}

int main() {
  std::cout << "\n======================= RUN TEST1 =======================\n";
  test1();
  std::cout << "\n======================= RUN TEST2 =======================\n";
  test2();
  return 0;
}
