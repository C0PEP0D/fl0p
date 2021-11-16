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
#include "fl0p/unstationary.h"

const unsigned int DIM = 3;

using TypeScalar = double;
using TypeVectorScalar = Eigen::Matrix<TypeScalar, 1, 1>;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
using TypeMatrix = Eigen::Matrix<TypeScalar, DIM, DIM>;
template<typename... Args>
using TypeRef = Eigen::Ref<Args...>;

template<typename ...Args>
using TypeContainer = std::vector<Args...>;
using TypeMesh = m0sh::Uniform<TypeVector, TypeRef, TypeContainer>;
using TypeMeshSub = m0sh::StructuredSub<TypeVector, TypeRef, TypeContainer>;
using TypeTimeMesh = m0sh::Uniform<TypeVectorScalar, TypeRef, TypeContainer>;
using TypeTimeMeshSub = m0sh::StructuredSub<TypeVectorScalar, TypeRef, TypeContainer>;
using TypeFlow = fl0w::fl0p::Unstationary<TypeVector, TypeMatrix, TypeRef, TypeMesh, TypeContainer, TypeMeshSub, TypeTimeMesh, TypeTimeMeshSub, TypeVectorScalar, v0l::FileData>;

void print(const TypeFlow& flow, const TypeVector& x, const TypeScalar& t) {
    std::cout << std::endl;
    std::cout << "flow.getVelocity(" << x.transpose() << ", " << t << ") = \n" << flow.getVelocity(x, t).transpose() << std::endl;
    //std::cout << "flow.getJacobian(" << x.transpose() << ", " << t << ") = \n" << flow.getJacobian(x, t) << std::endl;
    //std::cout << "flow.getVorticity(" << x.transpose() << ", " << t << ") = \n" << flow.getVorticity(x, t).transpose() << std::endl;
    //std::cout << "flow.getAcceleration(" << x.transpose() << ", " << t << ") = \n" << flow.getAcceleration(x, t).transpose() << std::endl;
    std::cout << std::endl;
}

int main () {
    std::vector<std::vector<v0l::FileData<float>>> velocity(2);
    for(std::size_t i = 0; i < velocity.size(); i++) {
        std::string fileName = std::string("../data/v050") + std::to_string(i) + ".vtk";
        velocity[i].emplace_back(fileName, 0);
        velocity[i].emplace_back(fileName, 1);
        velocity[i].emplace_back(fileName, 2);
    }
    // mesh data
    // // create lengths
    TypeVector origin;
    std::vector<double> lengths;
    for(std::size_t i = 0; i < velocity[0][0].meta.dimensions.size(); i++) {
        lengths.push_back(velocity[0][0].meta.dimensions[i] * velocity[0][0].meta.spacing[i]);
        origin[i] = velocity[0][0].meta.origin[i];
    }
    // // flow
    std::cout << "building flow..." << std::endl;
    TypeFlow flow(std::make_shared<TypeMesh>(velocity[0][0].meta.dimensions, lengths, origin, TypeContainer<bool>(DIM, true)), velocity, 4, std::make_shared<TypeTimeMesh>(std::vector<std::size_t>(1, 1), std::vector<double>(1, 0.02), TypeVectorScalar(0.0), std::vector<bool>(1, false)), 1);
    std::cout << "flow built !" << std::endl;
    // print
    print(flow, TypeVector({-0.5, -0.5, -0.5}), 0.0);
    print(flow, TypeVector({0.0, 0.0, 0.0}), 0.0);
    print(flow, TypeVector({0.5, 0.5, 0.5}), 0.0);
    print(flow, TypeVector({-0.5, -0.5, -0.5}), 0.1);
    print(flow, TypeVector({0.0, 0.0, 0.0}), 0.1);
    print(flow, TypeVector({0.5, 0.5, 0.5}), 0.1);
}
