#include <iostream>
#include "SoftBody.h++"


//------------------------------------------------------------------------------
// SoftBody::SoftBody
//------------------------------------------------------------------------------
SoftBody::SoftBody()
: _numVertices(0),
  _numEdges(0),
  _numFaces(0),
  _stacks(0),
  _slices(0),
  _fixNorthPole(false),
  _fixSouthPole(false),
  _pressure(0.0f),
  _temperature(0.0f),
  _volume(0.0f),
  _volumeInverse(0.0f)
{
}

//------------------------------------------------------------------------------
// SoftBody::~SoftBody
//------------------------------------------------------------------------------
SoftBody::~SoftBody()
{
    if (_vertices) {
        delete [] _vertices;
        _vertices = 0;
    }
    if (_edges) {
        delete [] _edges;
        _edges = 0;
    }
    if (_faces) {
        delete [] _faces;
        _faces = 0;
    }
}

//------------------------------------------------------------------------------
// SoftBody::GenerateSphere
//------------------------------------------------------------------------------
void
SoftBody::GenerateSphere(const float radius,
                         const uint32_t slices,
                         const uint32_t stacks,
                         const float mass,
                         const float springConstant)
{
    ClearDataStructures();
    std::cout << "Generating a soft sphere ... \n";

    _slices = slices;
    _stacks = stacks;

    Cartesian3f p;

    //
    // Generate the vertices.
    //
    _numVertices = slices * (stacks - 1) + 2;
    _vertices = new PointMassf[_numVertices];
    assert(_vertices);

    std::cout << "\t" << _numVertices << " point masses anticipated.\n";

    uint32_t index = 0;

    //
    // Generate the North pole first.
    //
    _vertices[index] = PointMassf(mass,
                                  Cartesian3f(0, radius, 0),
                                  Cartesian3f(0, radius, 0),
                                  Vector3f(0, 0, 0),
                                  Vector3f(0, 0, 0));

    //
    // Generate the vertices on lines of latitude
    //
    index++;
    float phiStep = PI / stacks;
    float thetaStep = 2.0f * PI / slices;
    for (float phi = phiStep; phi < PI; phi += phiStep) {
        for (float theta = 0.0f; theta < 2.0f * PI; theta += thetaStep) {
            p.x = radius * FloatUtilf::Cos(theta) * FloatUtilf::Sin(phi);
            p.z = radius * FloatUtilf::Sin(theta) * FloatUtilf::Sin(phi);
            p.y = radius * FloatUtilf::Cos(phi);
            _vertices[index] = PointMassf(mass,
                                          p,
                                          p,
                                          Vector3f(0, 0, 0),
                                          Vector3f(0, 0, 0));
            index++;
        }
    }

    //
    // Generate the South pole last.
    //
    _vertices[index] = PointMassf(mass,
                                  Cartesian3f(0, -radius, 0),
                                  Cartesian3f(0, -radius, 0),
                                  Vector3f(0, 0, 0),
                                  Vector3f(0, 0, 0));

    std::cout << "\t" << (index + 1) << " point masses created.\n";

    //
    // Generate the springs.
    //
    _numEdges = 3 * slices * (stacks - 1);
    std::cout << "\t" << _numEdges << " springs anticipated.\n";
    _edges = new SpringEdgef[_numEdges];
    assert(_edges);

    uint32_t i(0), j(0);
    uint32_t edgeIndex = 0;
    SpringEdgef edge;

    //
    // First add all the edges from the North pole to the first
    // ring of latitude.
    //
    for (i = 1; i <= slices; ++i) {
        edge.vertices[0] = 0;
        edge.vertices[1] = i;
        edge.springConstant = springConstant;
        edge.restLength = Length(_vertices[edge.vertices[0]].position
                                 - _vertices[edge.vertices[1]].position);
        _edges[edgeIndex] = edge;
        ++edgeIndex;
    }

    //
    // Now add the first line of latitude below the North pole.
    //
    for (i = 1; i <= slices; ++i) {
        if (i == slices) {
            edge.vertices[0] = i;
            edge.vertices[1] = 1;
        } else {
            edge.vertices[0] = i;
            edge.vertices[1] = i + 1;
        }
            edge.springConstant = springConstant;
            edge.restLength = Length(_vertices[edge.vertices[0]].position
                                    - _vertices[edge.vertices[1]].position);
            _edges[edgeIndex] = edge;
            ++edgeIndex;
    }

    //
    // Now add each slice in the middle.
    //
    uint32_t vertexIndices[4] = { 0 };
    for (j = 0; j < stacks - 2; ++j) {
        for (i = 0; i < slices; ++i) {
            if (i == (slices - 1)) {
                vertexIndices[0] = (j + 1) * slices + i + 1;
                vertexIndices[1] = j * slices + i + 1;
                vertexIndices[2] = (j + 1) * slices + 1;
                vertexIndices[3] = j * slices + 1;
            } else {
                vertexIndices[0] = (j + 1) * slices + i + 1;
                vertexIndices[1] = j * slices + i + 1;
                vertexIndices[2] = (j + 1) * slices + i + 2;
                vertexIndices[3] = j * slices + i + 2;
            }

            //
            // There are three edges generated per each vertex.
            // There are *five* total edges incident on each vertex
            // with the exception of the poles.
            //
            edge.vertices[0] = vertexIndices[0];
            edge.vertices[1] = vertexIndices[1];
            edge.springConstant = springConstant;
            edge.restLength = Length(_vertices[edge.vertices[0]].position
                                    - _vertices[edge.vertices[1]].position);
            _edges[edgeIndex] = edge;
            ++edgeIndex;

            edge.vertices[0] = vertexIndices[0];
            edge.vertices[1] = vertexIndices[2];
            edge.springConstant = springConstant;
            edge.restLength = Length(_vertices[edge.vertices[0]].position
                                    - _vertices[edge.vertices[1]].position);
            _edges[edgeIndex] = edge;
            ++edgeIndex;

            edge.vertices[0] = vertexIndices[0];
            edge.vertices[1] = vertexIndices[3];
            edge.springConstant = springConstant;
            edge.restLength = Length(_vertices[edge.vertices[0]].position
                                    - _vertices[edge.vertices[1]].position);
            _edges[edgeIndex] = edge;
            ++edgeIndex;
        }
    }

    //
    // Create the edges emanating from the last line of latitude and
    // proceeding to the South pole.
    //
    for (i = 0; i < slices; ++i) {
        edge.vertices[0] = _numVertices - 1;
        edge.vertices[1] = _numVertices - i - 2;
        edge.springConstant = springConstant;
        edge.restLength = Length(_vertices[_numVertices - 1].position
                                 -_vertices[_numVertices - i - 2].position);
        _edges[edgeIndex] = edge;
        ++edgeIndex;
    }

    std::cout << "\t" << edgeIndex << " springs created.\n";

    //
    // Create the faces.
    //
    _numFaces = 2 * slices + (stacks - 2) * slices * 2;
    std::cout << "\t" << _numFaces << " faces anticipated.\n";
    _faces = new Facef[_numFaces];
    assert(_faces);

    uint32_t faceIndex = 0;

    //
    // Add the fan of triangles around the North pole.
    //
    for (i = 1; i <= slices; ++i) {
        if (i == slices) {
            _faces[faceIndex].vertices[0] = 0;
            _faces[faceIndex].vertices[1] = i;
            _faces[faceIndex].vertices[2] = 1;
            faceIndex++;
        } else {
            _faces[faceIndex].vertices[0] = 0;
            _faces[faceIndex].vertices[1] = i;
            _faces[faceIndex].vertices[2] = i + 1;
            faceIndex++;
        }
    }

    //
    // Add all the faces in the middle slices.
    //
    for (j = 0; j < stacks - 2; ++j) {
        for (i = 0; i < slices; ++i) {
            if (i == (slices - 1)) {
                vertexIndices[0] = (j + 1) * slices + i + 1;
                vertexIndices[1] = j * slices + i + 1;
                vertexIndices[2] = (j + 1) * slices + 1;
                vertexIndices[3] = j * slices + 1;
            } else {
                vertexIndices[0] = (j + 1) * slices + i + 1;
                vertexIndices[1] = j * slices + i + 1;
                vertexIndices[2] = (j + 1) * slices + i + 2;
                vertexIndices[3] = j * slices + i + 2;
            }

            //
            // Four vertices gives us two triangles.
            //
            _faces[faceIndex].vertices[0] = vertexIndices[1];
            _faces[faceIndex].vertices[1] = vertexIndices[0];
            _faces[faceIndex].vertices[2] = vertexIndices[3];
            faceIndex++;

            _faces[faceIndex].vertices[0] = vertexIndices[3];
            _faces[faceIndex].vertices[1] = vertexIndices[0];
            _faces[faceIndex].vertices[2] = vertexIndices[2];
            faceIndex++;
        }
    }

    //
    // Create the faces on the fan of triangles along the South pole.
    //
    for (i = 1; i <= slices; ++i) {
        if (i == slices) {
            _faces[faceIndex].vertices[0] = _numVertices - 1;
            _faces[faceIndex].vertices[1] = _numVertices - slices - 1;
            _faces[faceIndex].vertices[2] = _numVertices - 2;
            faceIndex++;
        } else {
            _faces[faceIndex].vertices[0] = _numVertices - 1;
            _faces[faceIndex].vertices[1] = _numVertices - 1 - i;
            _faces[faceIndex].vertices[2] = _numVertices - 2 - i;
            faceIndex++;
        }
    }

    std::cout << "\t" << faceIndex << " faces created.\n";

    //
    // Normalize the normal vectors.
    //
    for (i = 0; i < _numVertices; ++i) {
        Normalize(_vertices[i].normal);
    }

    _surface = SURFACE_SPHERE;
    std::cout << "\nDone." << std::endl;
}

//------------------------------------------------------------------------------
// SoftBody::GenerateTorus
//------------------------------------------------------------------------------
void
SoftBody::GenerateTorus(const float innerRadius,
                        const float tubeRadius,
                        const uint32_t thetaSlices,
                        const uint32_t phiSlices,
                        const float mass,
                        const float springConstant)
{
    ClearDataStructures();
    std::cout << "Generating a soft torus ... \n";

    _thetaSlices = thetaSlices;
    _phiSlices = phiSlices;

    Cartesian3f p;
    Cartesian3f c;
    Vector3f n;

    //
    // Generate the vertices.
    //
    _numVertices = thetaSlices * phiSlices;
    _vertices = new PointMassf[_numVertices];
    assert(_vertices);

    std::cout << "\t" << _numVertices << " point masses anticipated.\n";

    uint32_t index = 0;
    float phiStep = 2.0f * PI / phiSlices;
    float thetaStep = 2.0f * PI / thetaSlices;
    for (float theta = 0.0f; theta < 2.0f * PI; theta += thetaStep) {
        for (float phi = 0.0f; phi < 2.0f * PI; phi += phiStep) {
            //
            // Compute the position of the vertex.
            //
            float tubeCosPhi = tubeRadius * FloatUtilf::Cos(phi);
            p.x = (innerRadius + tubeCosPhi) * FloatUtilf::Cos(theta);
            p.z = (innerRadius + tubeCosPhi) * FloatUtilf::Sin(theta);
            p.y = tubeRadius * FloatUtilf::Sin(phi);

            //
            // Compute the normal.
            //
            c.x = innerRadius * FloatUtilf::Cos(theta);
            c.y = innerRadius * FloatUtilf::Sin(theta);
            c.z = 0.0f;
            n = p - c;
            _vertices[index] = PointMassf(mass,
                                          p,
                                          n,
                                          Vector3f(0, 0, 0),
                                          Vector3f(0, 0, 0));
            ++index;
        }
    }
    std::cout << "\t" << index << " point masses created.\n";

    //
    // Generate the springs.
    //
    _numEdges = 3 * thetaSlices * phiSlices;
    std::cout << "\t" << _numEdges << " springs anticipated.\n";
    _edges = new SpringEdgef[_numEdges];
    assert(_edges);

    uint32_t    i(0), j(0); 
    uint32_t    edgeIndex = 0;
    SpringEdgef edge;
    uint32_t    vertexIndices[4] = { 0 };
    for (i = 0; i < thetaSlices; ++i) {
        for (j = 0; j < phiSlices; ++j) {
            if (i == thetaSlices - 1) {
                vertexIndices[0] = i * phiSlices + j;
                vertexIndices[2] = j;
                if (j == phiSlices - 1) {
                    vertexIndices[1] = i * phiSlices;
                    vertexIndices[3] = 0;
                } else {
                    vertexIndices[1] = i * phiSlices + (j + 1);
                    vertexIndices[3] = (j + 1);
                }
            } else {
                vertexIndices[0] = i * phiSlices + j;
                vertexIndices[2] = (i + 1) * phiSlices + j;
                if (j == phiSlices - 1) {
                    vertexIndices[1] = i * phiSlices;
                    vertexIndices[3] = (i + 1) * phiSlices;
                } else {
                    vertexIndices[1] = i * phiSlices + (j + 1);
                    vertexIndices[3] = (i + 1) * phiSlices + (j + 1);
                }
            }

            //
            // There are three edges generated per each vertex.
            // There are *five* total edges incident on each vertex.
            //
            edge.vertices[0] = vertexIndices[0];
            edge.vertices[1] = vertexIndices[1];
            edge.springConstant = springConstant;
            edge.restLength = Length(_vertices[edge.vertices[0]].position
                                    - _vertices[edge.vertices[1]].position);
            _edges[edgeIndex] = edge;
            ++edgeIndex;

            edge.vertices[0] = vertexIndices[0];
            edge.vertices[1] = vertexIndices[3];
            edge.springConstant = springConstant;
            edge.restLength = Length(_vertices[edge.vertices[0]].position
                                    - _vertices[edge.vertices[1]].position);
            _edges[edgeIndex] = edge;
            ++edgeIndex;

            edge.vertices[0] = vertexIndices[0];
            edge.vertices[1] = vertexIndices[2];
            edge.springConstant = springConstant;
            edge.restLength = Length(_vertices[edge.vertices[0]].position
                                    - _vertices[edge.vertices[1]].position);
            _edges[edgeIndex] = edge;
            ++edgeIndex;
        }
    }

    std::cout << "\t" << edgeIndex << " springs created.\n";

    //
    // Generate the faces.
    //
    _numFaces = 2 * thetaSlices * phiSlices;
    std::cout << "\t" << _numFaces << " faces anticipated.\n";
    _faces = new Facef[_numFaces];
    assert(_faces);

    uint32_t faceIndex = 0;
    //
    // Add all the faces in the middle slices.
    //
    for (i = 0; i < thetaSlices; ++i) {
        for (j = 0; j < phiSlices; ++j) {
            if (i == thetaSlices - 1) {
                vertexIndices[0] = i * phiSlices + j;
                vertexIndices[2] = j;
                if (j == phiSlices - 1) {
                    vertexIndices[1] = i * phiSlices;
                    vertexIndices[3] = 0;
                } else {
                    vertexIndices[1] = i * phiSlices + (j + 1);
                    vertexIndices[3] = (j + 1);
                }
            } else {
                vertexIndices[0] = i * phiSlices + j;
                vertexIndices[2] = (i + 1) * phiSlices + j;
                if (j == phiSlices - 1) {
                    vertexIndices[1] = i * phiSlices;
                    vertexIndices[3] = (i + 1) * phiSlices;
                } else {
                    vertexIndices[1] = i * phiSlices + (j + 1);
                    vertexIndices[3] = (i + 1) * phiSlices + (j + 1);
                }
            }

            //
            // Four vertices gives us two triangles.
            //
            _faces[faceIndex].vertices[0] = vertexIndices[1];
            _faces[faceIndex].vertices[1] = vertexIndices[0];
            _faces[faceIndex].vertices[2] = vertexIndices[3];
            faceIndex++;

            _faces[faceIndex].vertices[0] = vertexIndices[3];
            _faces[faceIndex].vertices[1] = vertexIndices[0];
            _faces[faceIndex].vertices[2] = vertexIndices[2];
            faceIndex++;
        }
    }

    std::cout << "\t" << faceIndex << " faces created.\n";

    //
    // Normalize the normal vectors.
    //
    for (i = 0; i < _numVertices; ++i) {
        Normalize(_vertices[i].normal);
    }

    _surface = SURFACE_TORUS;
    std::cout << "\nDone." << std::endl;
}

#define RENDER_POINTS   0
#define RENDER_EDGES    0
#define RENDER_FACES    1

//------------------------------------------------------------------------------
// SoftBody::Render
//------------------------------------------------------------------------------
void
SoftBody::Render(bool renderNormals,
                 bool renderVelocities,
                 bool renderForces,
                 bool renderWireframe)
{
#if RENDER_POINTS
    glBegin(GL_POINTS);
        glColor3f(1.0f, 0.0f, 0.0f);
        for (uint32_t i = 0; i < _numVertices; ++i) {
            glVertex3fv(_vertices[i].position.Ptr());
        }
    glEnd();
#endif

#if RENDER_EDGES
    glBegin(GL_LINES);
        glColor3f(1.0f, 0.0f, 0.0f);
        for (uint32_t i = 0; i < _numEdges; ++i) {
            glVertex3fv(_vertices[_edges[i].vertices[0]].position.Ptr());
            glVertex3fv(_vertices[_edges[i].vertices[1]].position.Ptr());
        }
    glEnd();
#endif

#if RENDER_FACES
    if (renderWireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    glBegin(GL_TRIANGLES);
        for (uint32_t i = 0; i < _numFaces; ++i) {
            glVertex3fv(_vertices[_faces[i].vertices[0]].position.Ptr());
            glNormal3fv(_vertices[_faces[i].vertices[0]].normal.Ptr());
            glVertex3fv(_vertices[_faces[i].vertices[1]].position.Ptr());
            glNormal3fv(_vertices[_faces[i].vertices[1]].normal.Ptr());
            glVertex3fv(_vertices[_faces[i].vertices[2]].position.Ptr());
            glNormal3fv(_vertices[_faces[i].vertices[2]].normal.Ptr());
        }
    glEnd();
#endif

    glDisable(GL_LIGHTING);
    if (renderNormals) {
        glBegin(GL_LINES);
        glColor3f(1.0f, 1.0f, 0.0f);
        for (uint32_t i = 0; i < _numVertices; ++i) {
            glVertex3fv(_vertices[i].position.Ptr());
            Normalize(_vertices[i].normal);
            Vector3f end = _vertices[i].normal + _vertices[i].position;
            glVertex3fv(end.Ptr());
        }
        glEnd();
    }
    if (renderVelocities) {
        glBegin(GL_LINES);
        glColor3f(1.0f, 1.0f, 0.0f);
        for (uint32_t i = 0; i < _numVertices; ++i) {
            glVertex3fv(_vertices[i].position.Ptr());
            Vector3f end = _vertices[i].velocity + _vertices[i].position;
            glVertex3fv(end.Ptr());
        }
        glEnd();
    }
    if (renderForces) {
        glBegin(GL_LINES);
        glColor3f(1.0f, 1.0f, 0.0f);
        for (uint32_t i = 0; i < _numVertices; ++i) {
            glVertex3fv(_vertices[i].position.Ptr());
            Vector3f end = _vertices[i].force + _vertices[i].position;
            glVertex3fv(end.Ptr());
        }
        glEnd();
    }
    glEnable(GL_LIGHTING);
}

//------------------------------------------------------------------------------
// SoftBody::AccumulateForces
//------------------------------------------------------------------------------
void
SoftBody::AccumulateForces()
{
    //
    // Inflate the mesh.
    //
    if (_pressure < _targetPressure) {
        _pressure += _targetPressure/900.0f;
    }

    //
    // Deflate it if needed.
    //
    if (_pressure > _targetPressure) {
        _pressure = _targetPressure;
    }

    //
    // Accumulate external force and buoyancy force.
    //
    uint32_t i(0);
    for (i = 0; i < _numVertices; ++i) {
        _vertices[i].force = Vector3f(0.0f, 0.0f, 0.0f);
        _vertices[i].force = _gravity
                            * _vertices[i].mass
                            * (_pressure - _targetPressure >= 0.0f);
        _vertices[i].force += _buoyancy;
    }

    //
    // Accumulate spring force.
    //
    Vector3f p[2];
    Vector3f v[2];
    for (i = 0; i < _numEdges; ++i) {
        p[0] = _vertices[_edges[i].vertices[0]].position;
        p[1] = _vertices[_edges[i].vertices[1]].position;
        v[0] = _vertices[_edges[i].vertices[0]].velocity;
        v[1] = _vertices[_edges[i].vertices[1]].velocity;

        float springLength = Length(p[0] - p[1]);
        if (springLength > 0.0f) {
            Vector3f x12 = p[0] - p[1];
            Vector3f v12 = v[0] - v[1];
            float forceMagnitude(0.0f);

            forceMagnitude = (springLength - _edges[i].restLength)
                            * _edges[i].springConstant
                            + (_damping / springLength)
                            * Dot(v12, x12);

            Vector3f force;
            force = x12 * (forceMagnitude / springLength);
            _vertices[_edges[i].vertices[0]].force -= force;
            _vertices[_edges[i].vertices[1]].force += force;
        }
    }

    //
    // Calculate the volume.
    //
    CalculateVolume();

    //
    // Accumulate forces from the Ideal Gas Law.
    //
    for (i = 0; i < _numFaces; ++i) {
        //
        // Calculate the pressure acting on this face.
        //
        float pressure = _pressure * _volumeInverse * _temperature;

        //
        // Calculate the area of the face.
        //
        Vector3f v1 = _vertices[_faces[i].vertices[1]].position
                    - _vertices[_faces[i].vertices[0]].position;
        Vector3f v2 = _vertices[_faces[i].vertices[2]].position
                    - _vertices[_faces[i].vertices[0]].position;
        Vector3f parallelogram;
        Cross(parallelogram, v1, v2);
        float area = 0.5f * Length(parallelogram);

        //
        // Add the pressure due to the force accumulator
        // of each point mass.
        //
        _vertices[_faces[i].vertices[0]].force += 
            _vertices[_faces[i].vertices[0]].normal * pressure * area * 60.0f;

        _vertices[_faces[i].vertices[1]].force += 
            _vertices[_faces[i].vertices[1]].normal * pressure * area * 60.0f;

        _vertices[_faces[i].vertices[2]].force += 
            _vertices[_faces[i].vertices[2]].normal * pressure * area * 60.0f;
    }
}

#define EULER           1
#define RUNGE_KUTTA4    0

//------------------------------------------------------------------------------
// SoftBody::SolveMotion
//------------------------------------------------------------------------------
void
SoftBody::SolveMotion(const float dt)
{
    for (uint32_t i = 0; i < _numVertices; ++i) {

        //
        // Euler integrate acceleration and velocity to get the
        // new velocity and position, respectively. Probably should
        // upgrade this integrator.
        //
        _vertices[i].force /= _vertices[i].mass;

        //
        // Do a forward euler step.
        //
#if EULER
        _vertices[i].velocity += _vertices[i].force * dt;
        _vertices[i].position += _vertices[i].velocity * dt;

        //
        // Do a Runge-Kutta4 step.
        //
#endif
#if RUNGE_KUTTA4
        _vertices[i].velocity += RungeKutta4(_vertices[i].velocity,
                                             _vertices[i].force,
                                             0.0001f);
        _vertices[i].position += RungeKutta4(_vertices[i].position,
                                             _vertices[i].velocity,
                                             0.0001f);
#endif

        //
        // Collide soft body with the clipping planes. If a vertex
        // crosses into the negative side of a plane, reflect its
        // velocity about the normal.
        //
        for (uint32_t j = 0; j < _clipPlanes.size(); ++j) {
            PlaneSide planeSide;
            planeSide = _clipPlanes[j].PointOnPlaneSide(_vertices[i].position);
            if (planeSide == PLANE_SIDE_NEGATIVE ||
                planeSide == PLANE_SIDE_ON)
            {
                //
                // Reflect the velocity about the normal using the law
                // of reflection.
                //
                Vector3f oldVelocity = -_vertices[i].velocity;
                float nDotV = Dot(_clipPlanes[j].normal, oldVelocity);

                if (nDotV < 0.0f) {
                    _vertices[i].velocity = _clipPlanes[j].normal
                                            * 2.0f * nDotV - oldVelocity;
                }

                //
                // Hard-code the resting condition.
                //
                if (_vertices[i].position.x < -3.8f) {
                    _vertices[i].position.x = -3.8f + 0.01f;
                }
                if (_vertices[i].position.x > 3.8f) {
                    _vertices[i].position.x = 3.8f - 0.01f;
                }
                if (_vertices[i].position.y < -3.8f) {
                    _vertices[i].position.y = -3.8f + 0.01f;
                }
                if (_vertices[i].position.y > 3.8f) {
                    _vertices[i].position.y = 3.8f - 0.01f;
                }
                if (_vertices[i].position.z < -3.8f) {
                    _vertices[i].position.z = -3.8f + 0.01f;
                }
                if (_vertices[i].position.z > 3.8f) {
                    _vertices[i].position.z = 3.8f - 0.01f;
                }
            }

            //
            // No mass may move faster than 3 units per time step
            // in any direction. This further promotes stability.
            //
            if (Length(_vertices[i].velocity) > 3.0f) {
                Normalize(_vertices[i].velocity);
                _vertices[i].velocity *= 3.0f;
            }
        }

        //
        // If either poles are fixed, they are not allowed to move.
        //
        if (_surface == SURFACE_SPHERE) {
            if (_fixNorthPole) {
                _vertices[0].velocity = Vector3f(0, 0, 0);
            }
            if (_fixSouthPole) {
                _vertices[_numVertices - 1].velocity = Vector3f(0, 0, 0);
            }
        }
    }


    //
    // Zero the current per-vertex normal vectors.
    //
    for (uint32_t i = 0; i < _numVertices; ++i) {
        _vertices[i].normal = Vector3f(0.0f, 0.0f, 0.0f);
    }

    //
    // Calculate the corrected normal vectors.
    //
    for (uint32_t i = 0; i < _numFaces; ++i) {
        Vector3f v1, v2;

        //
        // Calculate the normal of the face by taking the cross
        // product of vectors formed from two sides of the triangle.
        //
        v1 = _vertices[_faces[i].vertices[1]].position
            -_vertices[_faces[i].vertices[0]].position;
        v2 = _vertices[_faces[i].vertices[2]].position
            -_vertices[_faces[i].vertices[0]].position;

        Cross(_faces[i].normal, v2, v1);
        Normalize(_faces[i].normal);

        _vertices[_faces[i].vertices[0]].normal += _faces[i].normal;
        _vertices[_faces[i].vertices[1]].normal += _faces[i].normal;
        _vertices[_faces[i].vertices[2]].normal += _faces[i].normal;
    }

    for (uint32_t i = 0; i < _numVertices; ++i) {   
        Normalize(_vertices[i].normal);
    }
}

//------------------------------------------------------------------------------
// SoftBody::RungeKutta4
//------------------------------------------------------------------------------
Vector3f SoftBody::RungeKutta4(const Vector3f& x,
                               const Vector3f& v,
                               const float dt)
{
    const float half = 0.5f;
    const float third = float(1)/float(3);
    const float sixth = half * third;

    Vector3f k1 = v * dt;
    Vector3f k2 = (v + Vector3f(dt, dt, dt) * 0.5f) * dt;
    Vector3f k3 = (v + Vector3f(dt, dt, dt)) * dt;

    return x + k1 * sixth + k2 * third + k3 * sixth;
}

#define VOLUME_BY_BOUNDS    0
#define VOLUME_BY_GAUSS     1

//------------------------------------------------------------------------------
// SoftBody::CalculateVolume
//------------------------------------------------------------------------------
void
SoftBody::CalculateVolume()
{
    //
    // Calculate the volume.
    //
    _volume = 0.0f;

    uint32_t i = 0;

#if VOLUME_BY_BOUNDS
    //
    // Reset the bounds.
    //
    _mins = Vector3f(FloatUtilf::INFINITY,
                     FloatUtilf::INFINITY,
                     FloatUtilf::INFINITY);
    _maxes = Vector3f(-FloatUtilf::INFINITY,
                      -FloatUtilf::INFINITY,
                      -FloatUtilf::INFINITY);

    for (i = 0; i < _numVertices; ++i) {
        //
        // Set the new minimums and maximums if appropriate.
        //
        if (_vertices[i].position.x < _mins.x) {
            _mins.x = _vertices[i].position.x;
        }
        if (_vertices[i].position.y < _mins.y) {
            _mins.y = _vertices[i].position.y;
        }
        if (_vertices[i].position.x < _mins.z) {
            _mins.z = _vertices[i].position.z;
        }

        if (_vertices[i].position.x > _maxes.x) {
            _maxes.x = _vertices[i].position.x;
        }
        if (_vertices[i].position.y > _maxes.y) {
            _maxes.y = _vertices[i].position.y;
        }
        if (_vertices[i].position.z > _maxes.z) {
            _maxes.z = _vertices[i].position.z;
        }
    }

    _volume = (_maxes.x - _mins.x)
                * (_maxes.y - _mins.y)
                * (_maxes.z - _mins.z);
#endif

#if VOLUME_BY_GAUSS 
    for (i = 0; i < _numFaces; ++i) {

        //
        // Compute the average x coordinate.
        //
        float minX;
        minX = FloatUtilf::Min(_vertices[_faces[i].vertices[0]].position.x,
                               _vertices[_faces[i].vertices[1]].position.y,
                               _vertices[_faces[i].vertices[2]].position.z);

        float maxX;
        maxX = FloatUtilf::Max(_vertices[_faces[i].vertices[0]].position.x,
                               _vertices[_faces[i].vertices[1]].position.y,
                               _vertices[_faces[i].vertices[2]].position.z);

        float averageX = 0.5f * (maxX - minX);

        //
        // Compute the average normal.
        //
        Vector3f normal;
        normal =  _vertices[_faces[i].vertices[0]].normal
                + _vertices[_faces[i].vertices[1]].normal
                + _vertices[_faces[i].vertices[2]].normal;
        normal /= 3.0f;

        //
        // Compute the area.
        //
        Vector3f v1 = _vertices[_faces[i].vertices[1]].position
                    - _vertices[_faces[i].vertices[0]].position;
        Vector3f v2 = _vertices[_faces[i].vertices[2]].position
                    - _vertices[_faces[i].vertices[0]].position;
        Vector3f parallelogram;
        Cross(parallelogram, v1, v2);

        //
        // Add to the volume the product of the average x coordinate,
        // the x coordinate of the normal vector and the area of the
        // face.
        //
        _volume += FloatUtilf::Abs(averageX)
                 * FloatUtilf::Abs(normal.x)
                 * 0.5f * Length(parallelogram);
    }
#endif

    _volumeInverse = 1.0f / _volume;
}

//------------------------------------------------------------------------------
// SoftBody::ClearDataStructures
//------------------------------------------------------------------------------
void
SoftBody::ClearDataStructures()
{
    if (_vertices) {
        delete [] _vertices;
        _vertices = 0;
    }
    if (_edges) {
        delete [] _edges;
        _edges = 0;
    }
    if (_faces) {
        delete [] _faces;
        _faces = 0;
    }
    _numVertices = _numEdges = _numFaces = 0;
}
