#ifndef FL0P_STATIONARY_H
#define FL0P_STATIONARY_H
#pragma once

// include std
// // i/o
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <cstdlib>
#include <filesystem>
// // Data types
#include <tuple>
#include <memory> // shared_ptr
// include modules
#include "fl0w/flow.h"
#include "p0l/interpolation.h"

namespace fl0w {

namespace fl0p {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, typename TypeMesh, template<typename...> typename TypeContainer, typename TypeSubMesh, template<typename...> typename TypeDataContainer>
class Stationary : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        Stationary(const std::shared_ptr<TypeMesh>& p_sMesh, const TypeContainer<TypeDataContainer<float>>& p_velocity, const TypeContainer<TypeDataContainer<float>>& p_gradients, const std::size_t& p_order) : Flow<TypeVector, TypeMatrix, TypeRef>::Flow(), sMesh(p_sMesh), velocity(p_velocity), gradients(p_gradients), order(p_order) {
        }
        Stationary(const std::shared_ptr<TypeMesh>& p_sMesh, const TypeContainer<TypeDataContainer<float>>& p_velocity, const std::size_t& p_order) : Flow<TypeVector, TypeMatrix, TypeRef>::Flow(), sMesh(p_sMesh), velocity(p_velocity), order(p_order) {
            //computeJacobian();
        }
    public:
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override {
            TypeVector u;
            for(std::size_t i = 0; i < u.size(); i++) {
                u[i] = p0l::lagrangeMeshCell<TypeMesh, TypeDataContainer, float, TypeVector, TypeRef, TypeSubMesh>(sMesh, velocity[i], x, order+1);
            }
            return u;
        }
        TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const override {
            TypeMatrix grads;
            for(std::size_t i = 0; i < grads.reshaped().size(); i++) {
                grads.reshaped()[i] = p0l::lagrangeMeshCell<TypeMesh, TypeDataContainer, float, TypeVector, TypeRef, TypeSubMesh>(sMesh, gradients[i], x, order+1);
            }
            return grads;
        }
        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override {
            return TypeVector::Zero();
        }
    public:
        //void computeJacobian() {
        //    for(std::size_t k = 0; k < velocity.size(); k++) { 
        //        TypeContainer<int> ijk = sMesh->ijk(k);
        //        for(std::size_t j = 0; j < ijk.size(); j++) {
        //            TypeContainer<int> ijkp = ijk;
        //            ijkp[j] += 1;
        //            TypeContainer<int> ijkm = ijk;
        //            ijkm[j] -= 1;
        //            gradients[k].col(j) = (velocity[sMesh->index(ijkp)] - velocity[sMesh->index(ijkm)]) / (sMesh->x(ijkp)(j) - sMesh->x(ijkm)(j));
        //        }
        //    }
        //}
    public:
        TypeContainer<TypeDataContainer<float>> velocity;
        TypeContainer<TypeDataContainer<float>> gradients;
        std::shared_ptr<TypeMesh> sMesh;
        std::size_t order;
};

}

}

#endif
