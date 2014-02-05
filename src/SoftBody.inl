//------------------------------------------------------------------------------
// SoftBody::GetSurface
//------------------------------------------------------------------------------
inline SoftBody::Surface
SoftBody::GetSurface()
{
    return _surface;
}

//------------------------------------------------------------------------------
// SoftBody::GetInternalPressure
//------------------------------------------------------------------------------
inline float
SoftBody::GetInternalPressure()
{
    return _pressure;
}

//------------------------------------------------------------------------------
// SoftBody::SetInternalPressure
//------------------------------------------------------------------------------
inline void
SoftBody::SetInternalPressure(const float pressure)
{
    _pressure = pressure;
}

//------------------------------------------------------------------------------
// SoftBody::GetTargetPressure
//------------------------------------------------------------------------------
inline float
SoftBody::GetTargetPressure()
{
    return _targetPressure;
}

//------------------------------------------------------------------------------
// SetTargetPressure
//------------------------------------------------------------------------------
inline void
SoftBody::SetTargetPressure(const float targetPressure)
{
    _targetPressure = targetPressure;
}

//------------------------------------------------------------------------------
// SoftBody::GetInternalTemperature
//------------------------------------------------------------------------------
inline float
SoftBody::GetInternalTemperature()
{
    return _temperature;
}

//------------------------------------------------------------------------------
// SoftBody::SetInternalTemperature
//------------------------------------------------------------------------------
inline void
SoftBody::SetInternalTemperature(const float temperature)
{
    _temperature = temperature;
    _buoyancy = Vector3f(0.0f, temperature*0.1f, 0.0f);
}

//------------------------------------------------------------------------------
// SoftBody::SetSpringDamping
//------------------------------------------------------------------------------
inline void
SoftBody::SetSpringDamping(const float damping)
{
    _damping = damping;
}

//------------------------------------------------------------------------------
// SoftBody::SetSpringConstant
//------------------------------------------------------------------------------
inline void
SoftBody::SetSpringConstant(const float constant)
{
    for (uint32_t i = 0; i < _numEdges; ++i) {
        _edges[i].springConstant = constant;
    }
}

//------------------------------------------------------------------------------
// SoftBody::SetGravity
//------------------------------------------------------------------------------
inline void
SoftBody::SetGravity(const Vector3f& gravity)
{
    _gravity = gravity;
}

//------------------------------------------------------------------------------
// SoftBody::AddClipPlane
//------------------------------------------------------------------------------
inline void
SoftBody::AddClipPlane(const Planef& plane)
{
    _clipPlanes.push_back(plane);
}

//------------------------------------------------------------------------------
// SoftBody::ClearClipPlanes
//------------------------------------------------------------------------------
inline void
SoftBody::ClearClipPlanes()
{
    _clipPlanes.clear();
}

//------------------------------------------------------------------------------
// SoftBody::FixNorthPole
//------------------------------------------------------------------------------
inline void
SoftBody::FixNorthPole(const Cartesian3f& position)
{
    //
    // Only fix the pole if the surface is a sphere.
    // (Otherwise this method has no meaning.)
    //
    if (_surface == SURFACE_SPHERE) {
        _vertices[0].position = position;
        _vertices[0].velocity = Vector3f(0.0f, 0.0f, 0.0f);
        _fixNorthPole = true;
    }
}

//------------------------------------------------------------------------------
// SoftBody::IsNorthPoleFixed
//------------------------------------------------------------------------------
inline bool
SoftBody::IsNorthPoleFixed()
{
    return _fixNorthPole;
}

//------------------------------------------------------------------------------
// SoftBody::ReleaseNorthPole
//------------------------------------------------------------------------------
inline void
SoftBody::ReleaseNorthPole()
{
    _fixNorthPole = false;
}

//------------------------------------------------------------------------------
// SoftBody::FixSouthPole
//------------------------------------------------------------------------------
inline void
SoftBody::FixSouthPole(const Cartesian3f& position)
{
    //
    // Only fix the pole if the surface is a sphere.
    // (Otherwise this method has no meaning.)
    //
    if (_surface == SURFACE_SPHERE) {
        _vertices[_numVertices - 1].position = position;
        _vertices[_numVertices - 1].velocity = Vector3f(0.0f, 0.0f, 0.0f);
        _fixSouthPole = true;
    }
}

//------------------------------------------------------------------------------
// SoftBody::IsSouthPoleFixed
//------------------------------------------------------------------------------
inline bool
SoftBody::IsSouthPoleFixed()
{
    return _fixSouthPole;
}

//------------------------------------------------------------------------------
// SoftBody::ReleaseSouthPole
//------------------------------------------------------------------------------
inline void
SoftBody::ReleaseSouthPole()
{
    _fixSouthPole = false;
}
