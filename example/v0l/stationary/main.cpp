// Std includes
#include <iostream> // cout, endl
#include <vector>
#include <memory> // shared_ptr
#include <map>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "v0l/bin/file_data.h"
#include "m0sh/uniform.h"
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
using TypeMesh = m0sh::Uniform<TypeVector, TypeRef, TypeContainer>;
using TypeMeshSub = m0sh::StructuredSub<TypeVector, TypeRef, TypeContainer>;
using TypeFlow = fl0w::fl0p::Stationary<TypeVector, TypeMatrix, TypeRef, TypeMesh, TypeContainer, TypeMeshSub, v0l::FileData>;

void print(const TypeFlow& flow, const TypeVector& x, const TypeScalar& t) {
    std::cout << std::endl;
    std::cout << "flow.getVelocity(" << x.transpose() << ", " << t << ") = \n" << flow.getVelocity(x, t).transpose() << std::endl;
    std::cout << std::endl;
}

int main () { 
    v0l::FileData<float> vx("../data/v.vtk", 0);
    // mesh data
    // // create lengths
    TypeVector origin;
    std::vector<double> lengths;
    for(std::size_t i = 0; i < vx.meta.dimensions.size(); i++) {
        lengths.push_back(vx.meta.dimensions[i] * vx.meta.spacing[i]);
        origin[i] = vx.meta.origin[i];
    }
    // // other scalars
    v0l::FileData<float> vy("../data/v.vtk", 1);
    v0l::FileData<float> vz("../data/v.vtk", 2);
    // // flow
    std::cout << "building flow..." << std::endl;
    TypeFlow flow(std::make_shared<TypeMesh>(vx.meta.dimensions, lengths, origin, std::vector<bool>(DIM, true)), std::vector<v0l::FileData<float>>({vx, vy, vz}), 1);
    std::cout << "flow built !" << std::endl;

    // Print
    print(flow, TypeVector({-0.5, -0.5, -0.5}), 0.0);
    print(flow, TypeVector({0.0, 0.0, 0.0}), 0.0);
    print(flow, TypeVector({0.5, 0.5, 0.5}), 0.0);
}
