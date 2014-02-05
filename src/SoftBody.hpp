#ifndef __SoftBody_h__
#define __SoftBody_h__

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <fstream>
#include <vector>
#include "Primitives.hpp"
#include "Coordinates.hpp"
#include "Vector3.hpp"
#include "Plane.hpp"

//==============================================================================
// Class SoftBody
//==============================================================================
class SoftBody
{
public:
    //==========================================================================
    // Public Types
    //==========================================================================
    enum Surface
    {
        SURFACE_SPHERE,
        SURFACE_TORUS,
        SURFACE_FROM_FILE
    };

    //==========================================================================
    // Public Methods
    //==========================================================================
    SoftBody();
    ~SoftBody();

    //--------------------------------------------------------------------------
    // GenerateSphere
    //--------------------------------------------------------------------------
    void GenerateSphere(const float radius,
                        const uint32_t slices,
                        const uint32_t stacks,
                        const float mass,
                        const float springConstant);

    //--------------------------------------------------------------------------
    // GenerateTorus
    //--------------------------------------------------------------------------
    void GenerateTorus(const float innerRadius,
                       const float tubeRadius,
                       const uint32_t thetaSteps,
                       const uint32_t phiSteps,
                       const float mass,
                       const float springConstant);

    //--------------------------------------------------------------------------
    // GetSurface
    //--------------------------------------------------------------------------
    Surface GetSurface();

    //--------------------------------------------------------------------------
    // GetInternalPressure
    //--------------------------------------------------------------------------
    float GetInternalPressure();

    //--------------------------------------------------------------------------
    // SetInternalPressure
    //--------------------------------------------------------------------------
    void SetInternalPressure(const float pressure);

    //--------------------------------------------------------------------------
    // GetTargetPressure
    //--------------------------------------------------------------------------
    float GetTargetPressure();

    //--------------------------------------------------------------------------
    // SetTargetPressure
    //--------------------------------------------------------------------------
    void SetTargetPressure(const float pressure);

    //--------------------------------------------------------------------------
    // GetInternalTemperature
    //--------------------------------------------------------------------------
    float GetInternalTemperature();

    //--------------------------------------------------------------------------
    // SetInternalTemperature
    //--------------------------------------------------------------------------
    void SetInternalTemperature(const float temperature);

    //--------------------------------------------------------------------------
    // SetSpringDamping
    //--------------------------------------------------------------------------
    void SetSpringDamping(const float damping);

    //--------------------------------------------------------------------------
    // SetSpringConstant
    //--------------------------------------------------------------------------
    void SetSpringConstant(const float constant);

    //--------------------------------------------------------------------------
    // SetGravity
    //--------------------------------------------------------------------------
    void SetGravity(const Vector3f& f);

    //--------------------------------------------------------------------------
    // Render
    //--------------------------------------------------------------------------
    void Render(bool renderNormals,
                bool renderVelocites,
                bool renderForces,
                bool renderWireframe);

    //--------------------------------------------------------------------------
    // AccumulateForces
    //--------------------------------------------------------------------------
    void AccumulateForces();

    //--------------------------------------------------------------------------
    // SolveMotion
    //--------------------------------------------------------------------------
    void SolveMotion(const float dt);

    //--------------------------------------------------------------------------
    // AddClippingPlane
    //--------------------------------------------------------------------------
    void AddClipPlane(const Planef& plane);

    //--------------------------------------------------------------------------
    // ClearClipPlanes
    //--------------------------------------------------------------------------
    void ClearClipPlanes();

    //--------------------------------------------------------------------------
    // FixNorthPole
    //--------------------------------------------------------------------------
    void FixNorthPole(const Cartesian3f& position);

    //--------------------------------------------------------------------------
    // IsNorthPoleFixed
    //--------------------------------------------------------------------------
    bool IsNorthPoleFixed();

    //--------------------------------------------------------------------------
    // ReleaseNorthPole
    //--------------------------------------------------------------------------
    void ReleaseNorthPole();

    //--------------------------------------------------------------------------
    // FixSouthPole
    //--------------------------------------------------------------------------
    void FixSouthPole(const Cartesian3f& position);

    //--------------------------------------------------------------------------
    // IsSouthPoleFixed
    //--------------------------------------------------------------------------
    bool IsSouthPoleFixed();

    //--------------------------------------------------------------------------
    // ReleaseSouthPole
    //--------------------------------------------------------------------------
    void ReleaseSouthPole();

private:
    //==========================================================================
    // Private Methods
    //==========================================================================
    SoftBody(const SoftBody&) {}
    const SoftBody& operator=(const SoftBody&);
    Vector3f RungeKutta4(const Vector3f& x,
                         const Vector3f& v,
                         const float dt);
    void CalculateVolume();
    void ClearDataStructures();

    //==========================================================================
    // Private Data
    //==========================================================================
    uint32_t            _numVertices;
    PointMassf*         _vertices;

    uint32_t            _numEdges;
    SpringEdgef*        _edges;

    uint32_t            _numFaces;
    Facef*              _faces;

    float               _damping;

    Surface             _surface;

    //
    // Sphere attributes
    //
    uint32_t            _stacks;
    uint32_t            _slices;

    bool                _fixNorthPole;
    bool                _fixSouthPole;

    //
    // Torus attributes
    //
    uint32_t            _thetaSlices;
    uint32_t            _phiSlices;


    //
    // Force attributes
    //
    float               _pressure;
    float               _targetPressure;
    float               _temperature;
    float               _volume;
    float               _volumeInverse;
    Vector3f            _buoyancy;
    Vector3f            _gravity;

    Vector3f            _mins;
    Vector3f            _maxes;

    std::vector<Planef> _clipPlanes;
};

//==============================================================================
// Inline Implementation
//==============================================================================
#include "SoftBody.inl"


#endif
