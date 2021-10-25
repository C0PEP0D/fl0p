// Std includes
#include <iostream> // cout, endl
#include <vector>
#include <memory> // shared_ptr
#include <map>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "m0sh/regular.h"
#include "m0sh/structured_sub.h"
#include "fl0p/stationary.h"

const unsigned int DIM = 3;

using TypeScalar = double;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
using TypeMatrix = Eigen::Matrix<TypeScalar, DIM, DIM>;
template<typename... Args>
using TypeRef = Eigen::Ref<Args...>;

template<typename ...Args>
using TypeContainer = std::vector<Args...>;
using TypeMesh = m0sh::Regular<TypeVector, TypeRef, TypeContainer>;
using TypeMeshSub = m0sh::StructuredSub<TypeVector, TypeRef, TypeContainer>;
using TypeFlow = fl0w::fl0p::Stationary<TypeVector, TypeMatrix, TypeRef, TypeMesh, TypeContainer, TypeMeshSub, TypeContainer>;

void print(const TypeFlow& flow, const TypeVector& x, const TypeScalar& t) {
    std::cout << std::endl;
    std::cout << "flow.getVelocity(" << x.transpose() << ", " << t << ") = \n" << flow.getVelocity(x, t).transpose() << std::endl;
    std::cout << "flow.getJacobian(" << x.transpose() << ", " << t << ") = \n" << flow.getJacobian(x, t) << std::endl;
    std::cout << "flow.getVorticity(" << x.transpose() << ", " << t << ") = \n" << flow.getVorticity(x, t).transpose() << std::endl;
    std::cout << "flow.getAcceleration(" << x.transpose() << ", " << t << ") = \n" << flow.getAcceleration(x, t).transpose() << std::endl;
    std::cout << std::endl;
}

double f(const double x, const double y) {
    return x;
}

int main () { 
    // parameters
    std::size_t n = 8;
    // mesh
    std::vector<std::size_t> dimensions({n, n, n});
    std::vector<double> lengths({1.0, 1.0, 1.0});
    TypeVector origin({-0.5, -0.5, -0.5});
    // data
    std::vector<std::vector<float>> velocity(3, std::vector<float>(std::pow(n, DIM)));
    for(unsigned int i = 0; i < n; i++) {
        for(unsigned int j = 0; j < n; j++) {
            for(unsigned int k = 0; k < n; k++) {
                velocity[0][i + n*j + n*n*k] = (i + j + k)/3.0f;
                velocity[1][i + n*j + n*n*k] = (i + j + k)/3.0f;
                velocity[2][i + n*j + n*n*k] = (i + j + k)/3.0f;
            }
        }
    }
    // flow
    TypeFlow flow(std::make_shared<TypeMesh>(dimensions, lengths, origin), velocity, 4);
    // print
    print(flow, TypeVector({0.0, 0.0, 0.0}), 0.0);
}
