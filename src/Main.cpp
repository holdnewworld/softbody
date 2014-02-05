//==============================================================================
// Tim Thirion's Pressure Soft Body Simulation
// CS 530 - Purdue University
// Term Project
// Fall 2004
//==============================================================================
#include <iostream>
#include <sstream>
#include "Common.hpp"
#include "SoftBody.hpp"
#include "Clock.hpp"

//==============================================================================
// Global Variables
//==============================================================================
static uint32_t g_windowWidth           = 600;
static uint32_t g_windowHeight          = 600;
static bool     g_showFPS               = true;
static bool     g_renderNormals         = false;
static bool     g_renderVelocities      = false;
static bool     g_renderForces          = false;
static bool     g_renderWireframe       = false;
static float    g_sphereRadius          = 0.75f;
static uint32_t g_sphereSlices          = 32;
static uint32_t g_sphereStacks          = 32;
static float    g_torusInnerRadius      = 2.0f;
static float    g_torusTubeRadius       = 0.5f;
static uint32_t g_torusThetaSlices      = 32;
static uint32_t g_torusPhiSlices        = 32;
static float    g_pointMass             = 1.0f;
static float    g_springConstant        = 75.5f;
static float    g_springDamping         = 10.0f;
static Vector3f g_gravity               = Vector3f(0.0f, -1.0f, 0.0f);
//static float  g_internalPressure      = 0.0f;
static float    g_targetPressure        = 1.0f;
static float    g_internalTemperature   = 1.0f;
static Vector3f g_lightPosition         = Vector3f(0.0f, 0.0f, 0.0f);
static int      g_mousePos[2]           = { 0 };
static bool     g_mouseDown[2]          = { false };
static float    g_objectPosition[3]     = { 0.0f, 0.0f, -9.0f };
static float    g_objectRotation[3]     = { 0.0f };
static float    g_halfCubeSide          = 3.8f;
static SoftBody g_softBody;
static Clock    g_clock;

//==============================================================================
// Global Functions
//==============================================================================

//------------------------------------------------------------------------------
// OnOff
//------------------------------------------------------------------------------
inline std::string OnOff(const bool power)
{
    if (power) {
        return std::string("on");
    }
    return std::string("off");
}

//------------------------------------------------------------------------------
// PrintString
//------------------------------------------------------------------------------
void PrintString(const std::string s, float x, float y, float z)
{
    for (uint32_t i = 0; i < s.length(); ++i) {
        glColor3f(0.0f, 0.0f, 0.0f);
        glRasterPos3f(x, y, z);
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, s[i]);
        x -= (glutBitmapWidth(GLUT_BITMAP_8_BY_13, s[i]) * 0.0125f);
    }
}

//------------------------------------------------------------------------------
// PrintMenu
//------------------------------------------------------------------------------
void PrintMenu()
{
    std::cout << "Pressure Model Soft Body Simulation\n"
        << "Tim Thirion\n"
        << "CS 530 - Purdue University\n"
        << "Fall 2004\n"
        << "-------------------------------------------------\n\n"
        << "Menu\n"
        << "\tm...Print the menu.\n"
        << "\tr...Reset the simulation.\n"
        << "\t1...Change soft body to sphere.\n"
        << "\t2...Change soft body to torus.\n"
        << "\ti...Invert the direction of gravity.\n"
        << "\td...Decrease the spring damping.\n"
        << "\tD...Increase the spring damping.\n"
        << "\tk...Decrease the spring constant.\n"
        << "\tK...Increase the spring constant.\n"
        << "\tp...Decrease the target pressure.\n"
        << "\tP...Increase the target pressure.\n"
        << "\tt...Decrease the internal temperature.\n"
        << "\tT...Increase the internal temperature.\n"
        << "\tn...Toggle per-vertex surface normal rendering.\n"
        << "\tv...Toggle the rendering of the per-vertex velocities.\n"
        << "\tf...Toggle the rendering of the per-vertex forces.\n"
        << "\tw...Render wireframe.\n"
        << "\tc...Constrain the North pole of the soft sphere.\n"
        << "\tC...Constrain the South pole of the soft sphere.\n"
        << "\tq...Toggle measurement of frames per second.\n";
    std::cout << std::endl;
}

//------------------------------------------------------------------------------
// Initialize
//------------------------------------------------------------------------------
void Initialize()
{
    //
    // Print menu.
    //
    PrintMenu();

    //
    // Initialize OpenGL.
    //
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    //
    // Initialize the mesh.
    //
    g_softBody.GenerateSphere(g_sphereRadius,
            g_sphereSlices,
            g_sphereStacks,
            g_pointMass,
            g_springConstant);

    //
    // Set the spring damping (critically damped).
    //
    g_softBody.SetSpringDamping(g_springDamping);

    //
    // Add the force of gravity.
    //
    g_softBody.SetGravity(g_gravity);
    g_softBody.SetInternalPressure(0.0f);
    g_softBody.SetTargetPressure(45.0f);
    g_softBody.SetInternalTemperature(1.0f);

    //
    // Add the clipping planes. (Define a cube.)
    //

    // Bottom
    g_softBody.AddClipPlane(Planef(Vector3f(0.0f, 1.0f, 0.0f),
                Cartesian3f(0.0f, -g_halfCubeSide, 0.0f)));

    // Top
    g_softBody.AddClipPlane(Planef(Vector3f(0.0f, -1.0f, 0.0f),
                Cartesian3f(0.0f, g_halfCubeSide, 0.0f)));

    // Left
    g_softBody.AddClipPlane(Planef(Vector3f(1.0f, 0.0f, 0.0f),
                Cartesian3f(-g_halfCubeSide, 0.0f, 0.0f)));

    // Right
    g_softBody.AddClipPlane(Planef(Vector3f(-1.0f, 0.0f, 0.0f),
                Cartesian3f(g_halfCubeSide, 0.0f, 0.0f)));

    // Back
    g_softBody.AddClipPlane(Planef(Vector3f(0.0f, 0.0f, 1.0f),
                Cartesian3f(0.0f, 0.0f, -g_halfCubeSide)));

    // Front
    g_softBody.AddClipPlane(Planef(Vector3f(0.0f, 0.0f, -1.0f),
                Cartesian3f(0.0f, 0.0f, g_halfCubeSide)));

    // Cockeyed
    //  g_softBody.AddClipPlane(Planef(Vector3f(1.0f, 1.0f, 1.0f),
    //                                 Cartesian3f(0.0f, -3.5f, 0)));

    //
    // Add some lighting.
    //
    GLfloat diffuse[3] = { 0.73f, 0.0f, 0.0f };
    GLfloat specular[3] = { 0.25f, 0.25f, 0.85f };
    GLfloat shininess = 99.0f;

    GLfloat lightPosition[3] = { 0.0f, 0.0f, 10.0f };
    GLfloat lightAmbient[3] = { 0.0f, 0.0f, 0.0f };
    GLfloat lightColor[3] = { 1.0f, 1.0f, 1.0f };

    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);

    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
}

//------------------------------------------------------------------------------
// Reshape
//------------------------------------------------------------------------------
void Reshape(int width, int height)
{
    g_windowWidth = uint32_t(width);
    g_windowHeight = uint32_t(height);

    glViewport(0, 0, g_windowWidth, g_windowHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-1.0f, 1.0f, -1.0f, 1.0f, 1.0f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

//------------------------------------------------------------------------------
// Display
//------------------------------------------------------------------------------
void Display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    glTranslatef(g_objectPosition[0],
            g_objectPosition[1],
            g_objectPosition[2]);
    glRotatef(g_objectRotation[0], 1.0f, 0.0f, 0.0f);
    glRotatef(g_objectRotation[1], 0.0f, 1.0f, 0.0f);
    glRotatef(g_objectRotation[2], 0.0f, 0.0f, 1.0f);

    //
    // Render the floor.
    //
    GLfloat diffuse[3] = { 0.3f, 0.3f, 0.3f };
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);

    glBegin(GL_QUADS);
    glVertex3f(-g_halfCubeSide, -g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(0.0f, 1.0f, 0.0f);
    glVertex3f(-g_halfCubeSide, -g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(0.0f, 1.0f, 0.0f);
    glVertex3f( g_halfCubeSide, -g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(0.0f, 1.0f, 0.0f);
    glVertex3f( g_halfCubeSide, -g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(0.0f, 1.0f, 0.0f);
    glEnd();

    //
    // Render the body.
    //
    diffuse[0] = 0.73f;
    diffuse[1] = diffuse[2] = 0.0f;
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    g_softBody.Render(g_renderNormals,
            g_renderVelocities,
            g_renderForces,
            g_renderWireframe);

    //
    // Render the frame of the cube.
    //
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_QUADS);
    // Top
    glVertex3f(-g_halfCubeSide, g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(0.0f, -1.0f, 0.0f);
    glVertex3f(-g_halfCubeSide, g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(0.0f, -1.0f, 0.0f);
    glVertex3f( g_halfCubeSide, g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(0.0f, -1.0f, 0.0f);
    glVertex3f( g_halfCubeSide, g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(0.0f, -1.0f, 0.0f);

    // Left
    glVertex3f(-g_halfCubeSide, -g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-g_halfCubeSide, -g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-g_halfCubeSide,  g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-g_halfCubeSide,  g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(1.0f, 0.0f, 0.0f);

    // Right
    glVertex3f(g_halfCubeSide, -g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(-1.0f, 0.0f, 0.0f);
    glVertex3f(g_halfCubeSide, -g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(-1.0f, 0.0f, 0.0f);
    glVertex3f(g_halfCubeSide,  g_halfCubeSide, -g_halfCubeSide);
    glNormal3f(-1.0f, 0.0f, 0.0f);
    glVertex3f(g_halfCubeSide,  g_halfCubeSide,  g_halfCubeSide);
    glNormal3f(-1.0f, 0.0f, 0.0f);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


    glPopMatrix();
    if (g_showFPS) {
        glDisable(GL_LIGHTING);
        glPushMatrix();
        glTranslatef(0.0f, 0.0f, -2.0f);
        glScalef(-0.5f, 1.0f, 1.0f);
        std::stringstream os;
        os << "FPS: " << g_clock.GetFPS();
        PrintString(os.str(), 3.8f, 1.9f, 0.0f);
        glPopMatrix();
        glEnable(GL_LIGHTING);
    }
    glutSwapBuffers();
}

//------------------------------------------------------------------------------
// Keyboard
//------------------------------------------------------------------------------
void Keyboard(unsigned char key, int, int)
{
    switch (key)
    {
        //
        // Print the menu with 'm'.
        //
        case 'm':
            {
                PrintMenu();
            } break;

            //
            // Reset the simulation with '0'.
            //
        case '0':
            {
                if (g_softBody.GetSurface() == SoftBody::SURFACE_SPHERE) {
                    g_springConstant = 75.0f;
                    g_pointMass = 1.0f;
                    g_targetPressure = 1.0f;
                    g_softBody.GenerateSphere(g_sphereRadius,
                            g_sphereSlices,
                            g_sphereStacks,
                            g_pointMass,
                            g_springConstant);
                }
                if (g_softBody.GetSurface() == SoftBody::SURFACE_TORUS) {
                    g_springConstant = 2000.0f;
                    g_pointMass = 1.0f;
                    g_targetPressure = 200.0f;
                    g_softBody.GenerateTorus(g_torusInnerRadius,
                            g_torusTubeRadius,
                            g_torusThetaSlices,
                            g_torusPhiSlices,
                            g_pointMass,
                            g_springConstant);
                }

                std::cout << "Simulation restarted."
                    << std::endl;
            } break;

            //
            // Generate a soft sphere with '1'.
            //
        case '1':
            {
                g_springConstant = 75.0f;
                g_pointMass = 1.0f;
                g_targetPressure = 1.0f;
                g_softBody.GenerateSphere(g_sphereRadius,
                        g_sphereSlices,
                        g_sphereStacks,
                        g_pointMass,
                        g_springConstant);
            } break;

            //
            // Generate a soft torus with '2'.
            //
        case '2':
            {
                g_springConstant = 2000.0f;
                g_pointMass = 1.0f;
                g_targetPressure = 200.0f;
                g_softBody.GenerateTorus(g_torusInnerRadius,
                        g_torusTubeRadius,
                        g_torusThetaSlices,
                        g_torusPhiSlices,
                        g_pointMass,
                        g_springConstant);
            } break;

            //
            // Negate the direction of gravity with 'i'.
            //
        case 'i':
            {
                g_gravity = -g_gravity;
                g_softBody.SetGravity(g_gravity);
                std::cout << "Gravity inverted."
                    << std::endl;
            } break;

            //
            // Toggle rendering of normal vectors with 'n'.
            //
        case 'n':
            {
                g_renderNormals = !g_renderNormals;
                std::cout << "Rendering of normals toggled to "
                    << OnOff(g_renderNormals).c_str()
                    << "."
                    << std::endl;
            } break;

            //
            // Toggle the rendering of per-vertex velocities with 'v'.
            //
        case 'v':
            {
                g_renderVelocities = !g_renderVelocities;
                std::cout << "Rendering of velocity vectors toggled to "
                    << OnOff(g_renderVelocities).c_str()
                    << "."
                    << std::endl;
            } break;

            //
            // Toggle rendering of force vectors with 'f'.
            //
        case 'f':
            {
                g_renderForces = !g_renderForces;
                std::cout << "Rendering of forces toggled to "
                    << OnOff(g_renderForces).c_str()
                    << "."
                    << std::endl;
            } break;

            //
            // Toggle rendering as wireframe with 'w'.
            //
        case 'w':
            {
                g_renderWireframe = !g_renderWireframe;
                std::cout << "Rendering wireframe toggled to "
                    << OnOff(g_renderWireframe).c_str()
                    << "."
                    << std::endl;
            };

            //
            // Decrease the damping parameter with 'd'.
            //
        case 'd':
            {
                g_springDamping -= 0.25f;
                g_softBody.SetSpringDamping(g_springDamping);
                std::cout << "Spring damping decreased to "
                    << g_springDamping
                    << std::endl;
            } break;

            //
            // Increase the spring damping parameter with 'D'.
            //
        case 'D':
            {
                g_springDamping += 0.25f;
                g_softBody.SetSpringDamping(g_springDamping);
                std::cout << "Spring damping increased to "
                    << g_springDamping
                    << std::endl;
            } break;

            //
            // Decrease the spring constant with 'k'.
            //
        case 'k':
            {
                g_springConstant -= 0.25f;
                g_softBody.SetSpringConstant(g_springConstant);
                std::cout << "Spring constant decreased to "
                    << g_springConstant
                    << std::endl;
            } break;

            //
            // Increase the spring constant with 'K'.
            //
        case 'K':
            {
                g_springConstant += 0.25f;
                g_softBody.SetSpringConstant(g_springConstant);
                std::cout << "Sprint constant increased to "
                    << g_springConstant
                    << std::endl;
            } break;

            //
            // Decrease the target pressure with 'p'.
            //
        case 'p':
            {
                g_targetPressure -= 0.75f;
                g_softBody.SetTargetPressure(g_targetPressure);
                std::cout << "Target pressure decreased to "
                    << g_targetPressure
                    << std::endl;
            } break;

            //
            // Increase the target pressure with 'P'.
            //
        case 'P':
            {
                g_targetPressure += 0.75f;
                g_softBody.SetTargetPressure(g_targetPressure);
                std::cout << "Target pressure increased to "
                    << g_targetPressure
                    << std::endl;
            } break;

            //
            // Decrease the internal temperature with 't'.
            //
        case 't':
            {
                g_internalTemperature -= 0.25f;
                g_softBody.SetInternalTemperature(g_internalTemperature);
                std::cout << "Internal temperature decreased to "
                    << g_internalTemperature
                    << std::endl;
            } break;

            //
            // Increase the internal temperature with 'T'.
            //
        case 'T':
            {
                g_internalTemperature += 0.25f;
                g_softBody.SetInternalTemperature(g_internalTemperature);
                std::cout << "Internal temperature increased to "
                    << g_internalTemperature
                    << std::endl;
            } break;

            //
            // Constrain the North pole with 'c'.
            //
        case 'c':
            {
                SoftBody::Surface surface = g_softBody.GetSurface();
                if (surface == SoftBody::SURFACE_SPHERE) {
                    if (g_softBody.IsNorthPoleFixed()) {
                        g_softBody.ReleaseNorthPole();
                        std::cout << "North pole released."
                            << std::endl;
                    } else {
                        g_softBody.FixNorthPole(Cartesian3f(0.0f,
                                    g_halfCubeSide,
                                    0.0f));
                        std::cout << "North pole constrained to point "
                            << "(0.0, " << g_halfCubeSide << ", 0.0)"
                            << std::endl;
                    }
                }
            } break;

            //
            // Constrain the South pole with 'C'.
            //
        case 'C':
            {
                SoftBody::Surface surface = g_softBody.GetSurface();
                if (surface == SoftBody::SURFACE_SPHERE) {
                    if (g_softBody.IsSouthPoleFixed()) {
                        g_softBody.ReleaseSouthPole();
                        std::cout << "South pole released."
                            << std::endl;
                    } else {
                        g_softBody.FixSouthPole(Cartesian3f(0.0f,
                                    -g_halfCubeSide,
                                    0.0f));
                        std::cout << "South pole constrained to point "
                            << "(0.0, " << -g_halfCubeSide << ", 0.0)"
                            << std::endl;
                    }
                }
            } break;

            //
            // Toggle rendering FPS counter with 'q'.
            //
        case 'q':
            {
                g_showFPS = !g_showFPS;
                std::cout << "FPS counter toggled "
                    << OnOff(g_renderNormals).c_str()
                    << "."
                    << std::endl;
            } break;
    }
}

//------------------------------------------------------------------------------
// Motion
//------------------------------------------------------------------------------
void Motion(int x, int y)
{
    int difference[2];

    difference[0] = g_mousePos[0] - x;
    difference[1] = g_mousePos[1] - y;

    g_mousePos[0] = x;
    g_mousePos[1] = y;

    if (g_mouseDown[0])
    {
        g_objectRotation[0] -= difference[1];
        g_objectRotation[1] -= difference[0];
    }

    if (g_mouseDown[1])
        g_objectPosition[2] += 0.25f * difference[1];
    glutPostRedisplay();
}

//------------------------------------------------------------------------------
// Mouse
//------------------------------------------------------------------------------
void Mouse(int button, int down, int x, int y)
{
    switch (button)
    {
        case GLUT_LEFT_BUTTON:
            if (down == GLUT_DOWN)
                g_mouseDown[0] = true;
            else
                g_mouseDown[0] = false;
            break;
        case GLUT_RIGHT_BUTTON:
            if (down == GLUT_DOWN)
                g_mouseDown[1] = true;
            else
                g_mouseDown[1] = false;
            break;
    }

    g_mousePos[0] = x;
    g_mousePos[1] = y;
}

//------------------------------------------------------------------------------
// Idle
//------------------------------------------------------------------------------
void Idle()
{
    g_softBody.AccumulateForces();
    g_clock.Tick();
    g_softBody.SolveMotion(0.01f);
    glutPostRedisplay();
}

//------------------------------------------------------------------------------
// Shutdown
//------------------------------------------------------------------------------
void Shutdown()
{
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(g_windowWidth, g_windowHeight);
    glutCreateWindow("Soft Body");
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutIdleFunc(Idle);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse);
    glutMotionFunc(Motion);
    Initialize();
    glutMainLoop();
    Shutdown();
    return 0;
}
