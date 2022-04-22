#ifndef FL0P_UNSTATIONARY_H
#define FL0P_UNSTATIONARY_H
#pragma once

// include std
#include <tuple>
#include <memory> // shared_ptr
#include <cmath> // round
// include modules
#include "fl0p/stationary.h"

namespace fl0w {

namespace fl0p {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, typename TypeMesh, template<typename...> class TypeContainer, typename TypeSubMesh, typename TypeTimeMesh, typename TypeTimeSubMesh, typename TypeVectorScalar, template<typename...> typename TypeDataContainer>
class Unstationary : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        Unstationary() : Flow<TypeVector, TypeMatrix, TypeRef>::Flow() {
        }
        Unstationary(const std::shared_ptr<TypeMesh>& sMesh, const TypeContainer<TypeContainer<TypeDataContainer<float>>>& velocity, const TypeContainer<TypeContainer<TypeDataContainer<float>>>& gradients, const std::size_t& spaceOrder, const std::shared_ptr<TypeTimeMesh>& p_sTimeMesh, const std::size_t& p_timeOrder) : Flow<TypeVector, TypeMatrix, TypeRef>::Flow() {
            build(sMesh, velocity, gradients, spaceOrder, p_sTimeMesh, p_timeOrder);
        }
        Unstationary(const std::shared_ptr<TypeMesh>& sMesh, const TypeContainer<TypeContainer<TypeDataContainer<float>>>& velocity, const std::size_t& spaceOrder, const std::shared_ptr<TypeTimeMesh>& p_sTimeMesh, const std::size_t& p_timeOrder) : Flow<TypeVector, TypeMatrix, TypeRef>::Flow() {
            build(sMesh, velocity, spaceOrder, p_sTimeMesh, p_timeOrder);
        }
    public:
        void build(const std::shared_ptr<TypeMesh>& sMesh, const TypeContainer<TypeContainer<TypeDataContainer<float>>>& velocity, const TypeContainer<TypeContainer<TypeDataContainer<float>>>& gradients, const std::size_t& spaceOrder, const std::shared_ptr<TypeTimeMesh>& p_sTimeMesh, const std::size_t& p_timeOrder) {
            sTimeMesh = p_sTimeMesh;
            timeOrder = p_timeOrder;
            // data
            data.clear();
            for(std::size_t i = 0; i < velocity.size(); i++) {
                data.emplace_back(sMesh, velocity[i], gradients[i], spaceOrder);
            }
        }
        void build(const std::shared_ptr<TypeMesh>& sMesh, const TypeContainer<TypeContainer<TypeDataContainer<float>>>& velocity, const std::size_t& spaceOrder, const std::shared_ptr<TypeTimeMesh>& p_sTimeMesh, const std::size_t& p_timeOrder) {
            sTimeMesh = p_sTimeMesh;
            timeOrder = p_timeOrder;
            // data
            data.clear();
            for(std::size_t i = 0; i < velocity.size(); i++) {
                data.emplace_back(sMesh, velocity[i], spaceOrder);
            }
        }
    public:
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override {
            int tStartIndex = sTimeMesh->indexPoint(TypeVectorScalar(t)) - timeOrder/2;
            if(tStartIndex < 0) {
                tStartIndex = 0;
            } else if(tStartIndex > sTimeMesh->nPoints[0] - (timeOrder + 1)) {
                tStartIndex = sTimeMesh->nPoints[0] - (timeOrder + 1);
            }
            TypeContainer<TypeVector> velocity(timeOrder + 1);
            for(std::size_t i = 0; i < timeOrder + 1; i++) {
                velocity[i] = data[tStartIndex + i].getVelocity(x, t);
            }
            return p0l::lagrangeMeshPoint<TypeTimeSubMesh, TypeContainer, TypeVector, TypeVectorScalar, TypeRef>(TypeTimeSubMesh(TypeContainer<std::size_t>(1, timeOrder+1), TypeContainer<int>(1, tStartIndex), sTimeMesh), velocity, TypeVectorScalar(t));
        }

        TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const override {
            int tStartIndex = sTimeMesh->indexPoint(TypeVectorScalar(t)) - timeOrder/2;
            if(tStartIndex < 0) {
                tStartIndex = 0;
            } else if(tStartIndex > sTimeMesh->nPoints[0] - (timeOrder + 1)) {
                tStartIndex = sTimeMesh->nPoints[0] - (timeOrder + 1);
            }
            TypeContainer<TypeMatrix> gradients(timeOrder + 1);
            for(std::size_t i = 0; i < timeOrder + 1; i++) {
                gradients[i] = data[tStartIndex + i].getVelocityGradients(x, t);
            }
            return p0l::lagrangeMeshPoint<TypeTimeSubMesh, TypeContainer, TypeMatrix, TypeVectorScalar, TypeRef>(TypeTimeSubMesh(TypeContainer<std::size_t>(1, timeOrder+1), TypeContainer<int>(1, tStartIndex), sTimeMesh), gradients, TypeVectorScalar(t));
        }

        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override {
            return TypeVector::Zero(); // TODO
        }
    public:
        TypeContainer<Stationary<TypeVector, TypeMatrix, TypeRef, TypeMesh, TypeContainer, TypeSubMesh, TypeDataContainer>> data;
        std::shared_ptr<TypeTimeMesh> sTimeMesh;
        std::size_t timeOrder;
        //std::shared_ptr<TypeMesh> sMesh;
        //std::size_t spaceOrder;
};

}

}

#endif
